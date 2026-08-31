/**************************************************************************************/
/**                                                                                \n**/
/**            l  i  t  t  e  r  s  o  m  _  n  o  m  e  t  h  a  n  e  .  c       \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**               Vertical soil carbon distribution                                \n**/
/**                                                                                \n**/
/**               Litter and soil decomposition                                    \n**/
/**                                                                                \n**/
/**     Calculate daily litter decomposition using equation                        \n**/
/**       (1) dc/dt = -kc     where c=pool size, t=time,                           \n**/
/**           k=decomposition rate                                                 \n**/
/**     from (1),                                                                  \n**/
/**       (2) c = c0*exp(-kt) where c0=initial pool size                           \n**/
/**     from (2), decomposition in any month given by                              \n**/
/**       (3) delta_c = c0 - c0*exp(-k)                                            \n**/
/**     from (4)                                                                   \n**/
/**       (4) delta_c = c0*(1.0-exp(-k))                                           \n**/
/**     If LINEAR_DECAY is defined linear version of equations are used:           \n**/
/**       (3) delta_c = - c0*k                                                     \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"
#include "crop.h"
#include "agriculture.h"

#define CN_ratio_fast 8
#define CN_ratio_slow 12

static Real f_ph(Real ph)
{
  return 0.56 + atan(M_PI*0.45*(-5 + ph)) / M_PI;
} /* of 'f_ph' */

/* Between 57% and 98% of the exp() calls below are handed the argument they
   were handed last time: the responses do not change inside the litter loop
   and the PFT table holds only a handful of distinct decay constants.
   Returning a remembered result is exact -- the same argument gives the same
   bits -- and there are a hundred million of these calls in a three-year run.
   Not thread-safe; LPJmL divides the grid over MPI ranks and has no threads. */
static Real decay_memo(Real *lastx,Real *lasty,Real x)
{
  if(*lastx==x)
    return *lasty;
  *lastx=x;
  return *lasty=1.0-exp(-x);
}

static Real mx_agtop_leaf=NAN,my_agtop_leaf;
static Real mx_agtop_wood=NAN,my_agtop_wood;
static Real mx_agsub_leaf=NAN,my_agsub_leaf;
static Real mx_agsub_wood=NAN,my_agsub_wood;
static Real mx_slow=NAN,my_slow;
static Real mx_fast=NAN,my_fast;

/* q10_wood takes only a handful of distinct values over the whole PFT table --
   five in the standard parameter set, and 1.0 for every crop -- while the two
   exponents are fixed for the entire litter loop.  pow() was therefore being
   asked the same question again and again.  Eight slots covers the table with
   room to spare; anything beyond that just calls pow() as before. */
#define NQ10CACHE 8

Stocks littersom_nomethane(Stand *stand,                /**< pointer to stand data */
                           const Real gtemp_soil[NSOILLAYER], /**< respiration coefficents */
                           Real cellfrac_agr,           /**< stand fraction of agricultural cells (0..1) */
                           int npft,                    /**< number of natural PFTs */
                           int ncft,                    /**< number of crop PFTs */
                           const Config *config         /**< nitrogen enabled? */
                          ) /** \return decomposed carbon/nitrogen (g/m2) */
{
  Real response[NSOILLAYER];
  Real response_agtop_leaves,response_agtop_wood,response_agsub_leaves,response_agsub_wood,response_bg_litter,w_agtop;
  Real decay_slow,decay_fast;         /* decay fractions of the two soil pools */
  Real fmoist_bg,fmoist_agtop;        /* moisture response, below ground and above-ground top */
  Real q10exp_agsub,q10exp_agtop;     /* temperature exponents of the q10 response */
  Real cshift_c_fast,cshift_c_slow;   /* carbon shifted into the fast and slow pools */
  const struct Pftpar *litterpft;     /* PFT of the litter item being decomposed */
  int pftid;
  Real decay_litter;
  Pool flux_soil[LASTLAYER];
  Real decom,soil_cflux;
  Stocks decom_litter;
  Stocks decom_sum,flux;
  Real moist[NSOILLAYER];
  Real N_sum;
  Real n_immo;
  int i,p,l;
  Soil *soil;
  Real yedoma_flux;
  Real F_NO3=0;                /* soil nitrification rate gN *m-2*d-1*/
  Real F_N2O=0;                /* soil nitrification rate gN *m-2*d-1*/
  Real F_Nmineral;  /* net mineralization flux gN *m-2*d-1*/
  Real fac_wfps, fac_temp;
  Real q10base[NQ10CACHE],q10agsub[NQ10CACHE],q10agtop[NQ10CACHE];
  Real q10w,p_agsub,p_agtop;
  Real decay_bg;
  Real cshift_fast_c,cshift_slow_c,cshift_fast_n,cshift_slow_n;
  Real nmin_base,cn_excess,nimmo_fast,nimmo_slow;
  Bool has_c,has_n;
  int nq10=0,q;
#ifdef SAFE
  String line;
#endif
  Pft *pft;
  Pftcrop *crop;
  Irrigation *data;
  soil=&stand->soil;
  flux.nitrogen=0;
  foreachsoillayer(l)
    response[l]=0.0;
  decom_litter.carbon=decom_litter.nitrogen=soil_cflux=yedoma_flux=decom_sum.carbon=decom_sum.nitrogen=0.0;

  foreachsoillayer(l)
  {
    if(gtemp_soil[l]>0)
    {
      if (soil->wsats[l]-soil->ice_depth[l]-soil->ice_fw[l]-(soil->wpwps[l]*soil->ice_pwp[l])>epsilon)
        moist[l]=(soil->w[l]*soil->whcs[l]+(soil->wpwps[l]*(1-soil->ice_pwp[l]))+soil->w_fw[l])
                 /(soil->wsats[l]-soil->ice_depth[l]-soil->ice_fw[l]-(soil->wpwps[l]*soil->ice_pwp[l]));
      else
        moist[l]=epsilon;
      /* moist[l] must be in the inteval [0,1] */
      if (moist[l]<epsilon)
        moist[l]=epsilon;
      else if (moist[l]>1)
        moist[l]=1;

      response[l]=gtemp_soil[l]*(INTERCEPT+MOIST_3*(moist[l]*moist[l]*moist[l])+MOIST_2*(moist[l]*moist[l])+MOIST*moist[l]);
      if (response[l]<epsilon)
        response[l]=0.0;
      if (response[l]>1)
        response[l]=1.0;


      if(l<LASTLAYER)
      {
        if(stand->type->landusetype==NATURAL)
          getoutputindex(&stand->cell->output,RESPONSE_LAYER_NV,l,config)+=response[l];
        if(isagriculture(stand))
          getoutputindex(&stand->cell->output,RESPONSE_LAYER_AGR,l,config)+=response[l]*stand->frac/cellfrac_agr;
#ifdef SAFE
        if(soil->NO3[l]<-epsilon)
          fail(NEGATIVE_SOIL_NO3_ERR,TRUE,TRUE,"littersom: Negative soil NO3=%g in layer %d in cell (%s) before update",soil->NO3[l],l,sprintcoord(line,&stand->cell->coord));
        if(soil->NH4[l]<-epsilon)
          fail(NEGATIVE_SOIL_NH4_ERR,TRUE,TRUE,"littersom: Negative soil NH4=%g in layer %d in cell (%s) before update",soil->NH4[l],l,sprintcoord(line,&stand->cell->coord));
#endif
#ifdef LINEAR_DECAY
        flux_soil[l].slow.carbon=max(0,soil->pool[l].slow.carbon*param.k_soil10.slow*response[l]);
        flux_soil[l].fast.carbon=max(0,soil->pool[l].fast.carbon*param.k_soil10.fast*response[l]);
        flux_soil[l].slow.nitrogen=max(0,soil->pool[l].slow.nitrogen*param.k_soil10.slow*response[l]);
        flux_soil[l].fast.nitrogen=max(0,soil->pool[l].fast.nitrogen*param.k_soil10.fast*response[l]);
#else
        /* the same two decay fractions serve carbon and nitrogen; exp() is not
           free, and the compiler cannot merge the calls itself while libm may
           set errno */
        decay_slow=decay_memo(&mx_slow,&my_slow,param.k_soil10.slow*response[l]);
        decay_fast=decay_memo(&mx_fast,&my_fast,param.k_soil10.fast*response[l]);
        flux_soil[l].slow.carbon=max(0,soil->pool[l].slow.carbon*decay_slow);
        flux_soil[l].fast.carbon=max(0,soil->pool[l].fast.carbon*decay_fast);
        flux_soil[l].slow.nitrogen=max(0,soil->pool[l].slow.nitrogen*decay_slow);
        flux_soil[l].fast.nitrogen=max(0,soil->pool[l].fast.nitrogen*decay_fast);
#endif

/* TODO nitrogen limitation of decomposition including variable decay rates Brovkin,
      fn_som=flux_soil[l].fast.carbon/CN_ratio_fast+flux_soil[l].slow.carbon/CN_ratio_slow;

        if (fn_som>(flux_soil[l].fast.nitrogen+flux_soil[l].slow.nitrogen+soil->NH4[l]))
        {
          fnlim=max(0,(soil->NH4[l]+flux_soil[l].fast.nitrogen+flux_soil[l].slow.nitrogen)/fn_som);
        }*/
        if(soil->pool[l].slow.carbon>epsilon)
          soil->decay_rate[l].slow+=flux_soil[l].slow.carbon/soil->pool[l].slow.carbon;
        if(soil->pool[l].fast.carbon>epsilon)
          soil->decay_rate[l].fast+=flux_soil[l].fast.carbon/soil->pool[l].fast.carbon;

        soil->pool[l].slow.carbon-=flux_soil[l].slow.carbon;
        soil->pool[l].fast.carbon-=flux_soil[l].fast.carbon;
        soil->pool[l].slow.nitrogen-=flux_soil[l].slow.nitrogen;
        soil->pool[l].fast.nitrogen-=flux_soil[l].fast.nitrogen;
        soil_cflux+=flux_soil[l].slow.carbon+flux_soil[l].fast.carbon;
        F_Nmineral=flux_soil[l].slow.nitrogen+flux_soil[l].fast.nitrogen;
        soil->NH4[l]+=F_Nmineral;
#ifdef SAFE
        if(soil->NH4[l]<-epsilon)
          fail(NEGATIVE_SOIL_NH4_ERR,TRUE,TRUE,"Negative soil NH4=%g in layer %d in cell (%s) at mineralization",soil->NH4[l],l,sprintcoord(line,&stand->cell->coord));
#endif
        getoutput(&stand->cell->output,N_MINERALIZATION,config)+=F_Nmineral*stand->frac;
        if(isagriculture(stand))
          getoutput(&stand->cell->output,NMINERALIZATION_AGR,config)+=F_Nmineral*stand->frac;
        soil->k_mean[l].fast+=(param.k_soil10.fast*response[l]);
        soil->k_mean[l].slow+=(param.k_soil10.slow*response[l]);
#ifdef MICRO_HEATING
        soil->decomC[l]=flux_soil[l].slow.carbon+flux_soil[l].fast.carbon;
#endif
      }
      else
      {
        if (soil->YEDOMA>0.0 && response[l]>0.0)
        {
          yedoma_flux=soil->YEDOMA*(1.0-exp(-(K10_YEDOMA*response[l])));
          soil->YEDOMA-=yedoma_flux;
          soil_cflux+=yedoma_flux;
#ifdef MICRO_HEATING
          soil->decomC[l]+=yedoma_flux;
#endif
        }
      }
    }
  } /* end foreachsoillayer */

  /*
   *   Calculate daily decomposition rates (k, /month) as a function of
   *   temperature and moisture
   *
   */

  if(gtemp_soil[0]>0)
  {
    w_agtop=soil->litter.agtop_wcap>epsilon ? soil->litter.agtop_moist/soil->litter.agtop_wcap : moist[0];
    fmoist_bg=INTERCEPT+MOIST_3*(moist[0]*moist[0]*moist[0])+MOIST_2*(moist[0]*moist[0])+MOIST*moist[0];
    fmoist_agtop=INTERCEPT+MOIST_3*(w_agtop*w_agtop*w_agtop)+MOIST_2*(w_agtop*w_agtop)+MOIST*w_agtop;
    q10exp_agsub=(soil->temp[0]-10)/10.0;
    q10exp_agtop=(soil->litter.agtop_temp-10)/10.0;
    response_agsub_leaves=response[0];
    response_bg_litter=response[0];
    response_agtop_leaves=temp_response(soil->litter.agtop_temp,soil->amean_temp[0])*fmoist_agtop;
    /* neither factor changes inside the loop, so this is computed once */
#ifdef LINEAR_DECAY
    decay_bg=param.k_litter10*response_bg_litter;
#else
    decay_bg=1.0-exp(-(param.k_litter10*response_bg_litter));
#endif
    for(p=0;p<soil->litter.n;p++)
    {
      /* one pointer chase for the whole item, and one evaluation of the two
         c_shift terms that the pool update and three output blocks all share */
      litterpft=soil->litter.item[p].pft;
      pftid=litterpft->id;
      /* the same base and exponent return the same bits, so reusing an
         answer already computed in this loop is exact */
      q10w=litterpft->k_litter10.q10_wood;
      for(q=0;q<nq10;q++)
        if(q10base[q]==q10w)
          break;
      if(q<nq10)
      {
        p_agsub=q10agsub[q];
        p_agtop=q10agtop[q];
      }
      else
      {
        p_agsub=pow(q10w,q10exp_agsub);
        p_agtop=pow(q10w,q10exp_agtop);
        if(nq10<NQ10CACHE)
        {
          q10base[nq10]=q10w;
          q10agsub[nq10]=p_agsub;
          q10agtop[nq10]=p_agtop;
          nq10++;
        }
      }
      response_agsub_wood=p_agsub*fmoist_bg;
      response_agtop_wood=p_agtop*fmoist_agtop;

      decom_sum.carbon=decom_sum.nitrogen=0;
      /* agtop leaves */
#ifdef LINEAR_DECAY
      decay_litter=litterpft->k_litter10.leaf*response_agtop_leaves;
#else
      decay_litter=decay_memo(&mx_agtop_leaf,&my_agtop_leaf,litterpft->k_litter10.leaf*response_agtop_leaves);
#endif
      decom=soil->litter.item[p].agtop.leaf.carbon*decay_litter;
      soil->litter.item[p].agtop.leaf.carbon-=decom;
      decom_sum.carbon+=decom;
      decom=soil->litter.item[p].agtop.leaf.nitrogen*decay_litter;
      soil->litter.item[p].agtop.leaf.nitrogen-=decom;
      decom_sum.nitrogen+=decom;

      /* agtop wood */
#ifdef LINEAR_DECAY
      decay_litter=litterpft->k_litter10.wood*response_agtop_wood;
#else
      decay_litter=decay_memo(&mx_agtop_wood,&my_agtop_wood,litterpft->k_litter10.wood*response_agtop_wood);
#endif
      for(i=0;i<NFUELCLASS;i++)
      {
        decom=soil->litter.item[p].agtop.wood[i].carbon*decay_litter;
        soil->litter.item[p].agtop.wood[i].carbon-=decom;
        decom_sum.carbon+=decom;
        decom=soil->litter.item[p].agtop.wood[i].nitrogen*decay_litter;
        soil->litter.item[p].agtop.wood[i].nitrogen-=decom;
        decom_sum.nitrogen+=decom;
      }

      /* agsub leaves */
#ifdef LINEAR_DECAY
      decay_litter=litterpft->k_litter10.leaf*response_agsub_leaves;
#else
      decay_litter=decay_memo(&mx_agsub_leaf,&my_agsub_leaf,litterpft->k_litter10.leaf*response_agsub_leaves);
#endif
      decom=soil->litter.item[p].agsub.leaf.carbon*decay_litter;
      soil->litter.item[p].agsub.leaf.carbon-=decom;
      decom_sum.carbon+=decom;
      decom=soil->litter.item[p].agsub.leaf.nitrogen*decay_litter;
      soil->litter.item[p].agsub.leaf.nitrogen-=decom;
      decom_sum.nitrogen+=decom;

      /* agsub wood */
#ifdef LINEAR_DECAY
      decay_litter=litterpft->k_litter10.wood*response_agsub_wood;
#else
      decay_litter=decay_memo(&mx_agsub_wood,&my_agsub_wood,litterpft->k_litter10.wood*response_agsub_wood);
#endif
      for(i=0;i<NFUELCLASS;i++)
      {
        decom=soil->litter.item[p].agsub.wood[i].carbon*decay_litter;
        soil->litter.item[p].agsub.wood[i].carbon-=decom;
        decom_sum.carbon+=decom;
        decom=soil->litter.item[p].agsub.wood[i].nitrogen*decay_litter;
        soil->litter.item[p].agsub.wood[i].nitrogen-=decom;
        decom_sum.nitrogen+=decom;
      }

      /* bg litter */
      decay_litter=decay_bg;
      decom=soil->litter.item[p].bg.carbon*decay_litter;
      soil->litter.item[p].bg.carbon-=decom;
      decom_sum.carbon+=decom;
      decom=soil->litter.item[p].bg.nitrogen*decay_litter;
      soil->litter.item[p].bg.nitrogen-=decom;
      decom_sum.nitrogen+=decom;
      decom_litter.carbon+=decom_sum.carbon;
      decom_litter.nitrogen+=decom_sum.nitrogen;
      soil->decomp_litter_pft[pftid].carbon+=(1-param.atmfrac)*decom_sum.carbon;
      soil->decomp_litter_pft[pftid].nitrogen+=(1-param.atmfrac)*decom_sum.nitrogen;
      /* None of this changes with the layer: decom_sum is finished above and
         the fractions are parameters.  Each product keeps the grouping the
         expression had, so every result below is bit-for-bit what it was --
         the leading factors of a left-associated product are exactly what may
         be lifted out of it.  cn_ratio is refused at read time unless it is
         positive, so hoisting the division cannot divide by zero where the
         original would not have. */
      cshift_fast_c=param.fastfrac*(1-param.atmfrac)*decom_sum.carbon;
      cshift_slow_c=(1-param.fastfrac)*(1-param.atmfrac)*decom_sum.carbon;
      cshift_fast_n=param.fastfrac*(1-param.atmfrac)*decom_sum.nitrogen;
      cshift_slow_n=(1-param.fastfrac)*(1-param.atmfrac)*decom_sum.nitrogen;
      nmin_base=decom_sum.nitrogen*param.atmfrac;
      cn_excess=decom_sum.carbon/soil->par->cn_ratio-decom_sum.nitrogen;
      nimmo_fast=param.fastfrac*(1-param.atmfrac)*cn_excess;
      nimmo_slow=(1-param.fastfrac)*(1-param.atmfrac)*cn_excess;
      has_c=decom_sum.carbon>0;
      has_n=decom_sum.nitrogen>0;
      forrootsoillayer(l)
      {
        cshift_c_fast=cshift_fast_c*soil->c_shift[l][pftid].fast;
        cshift_c_slow=cshift_slow_c*soil->c_shift[l][pftid].slow;
        soil->pool[l].fast.carbon+=cshift_c_fast;
        soil->pool[l].slow.carbon+=cshift_c_slow;
        if(has_c && stand->type->landusetype==NATURAL)
        {
          getoutputindex(&stand->cell->output,CSHIFT_FAST_NV,l,config)+=cshift_c_fast;
          getoutputindex(&stand->cell->output,CSHIFT_SLOW_NV,l,config)+=cshift_c_slow;
        }
        /* WHEP: agricultural and managed-grassland counterparts. MUST be kept in
           step with the identical block in littersom.c: daily_littersom() calls
           THIS function whenever methane is off, which is the WHEP production
           configuration, so an addition made only in littersom.c is silently
           absent from every WHEP run while still working in any methane-enabled
           test. Stand-fraction weighted, unlike _nv, because these stand types
           occupy a fraction of the cell. */
        if(has_c && isagriculture(stand))
        {
          getoutputindex(&stand->cell->output,CSHIFT_FAST_AGR,l,config)+=cshift_c_fast*stand->frac;
          getoutputindex(&stand->cell->output,CSHIFT_SLOW_AGR,l,config)+=cshift_c_slow*stand->frac;
        }
        if(has_c && getlandusetype(stand)==GRASSLAND)
        {
          getoutputindex(&stand->cell->output,CSHIFT_FAST_MGRASS,l,config)+=cshift_c_fast*stand->frac;
          getoutputindex(&stand->cell->output,CSHIFT_SLOW_MGRASS,l,config)+=cshift_c_slow*stand->frac;
        }
        if(has_n)
        {
          soil->pool[l].slow.nitrogen+=cshift_slow_n*soil->c_shift[l][pftid].slow;
          soil->pool[l].fast.nitrogen+=cshift_fast_n*soil->c_shift[l][pftid].fast;
          /* NO3 and N2O from mineralization of organic matter */
          F_Nmineral=nmin_base*(param.fastfrac*soil->c_shift[l][pftid].fast+(1-param.fastfrac)*soil->c_shift[l][pftid].slow);
          soil->NH4[l]+=F_Nmineral;
#ifdef SAFE
          if(soil->NH4[l]<-epsilon)
            fail(NEGATIVE_SOIL_NH4_ERR,TRUE,TRUE,"Negative soil NH4=%g in layer %d in cell (%s) at mineralization",soil->NH4[l],l,sprintcoord(line,&stand->cell->coord));
#endif
          getoutput(&stand->cell->output,N_MINERALIZATION,config)+=F_Nmineral*stand->frac;
          if(isagriculture(stand))
            getoutput(&stand->cell->output,NMINERALIZATION_AGR,config)+=F_Nmineral*stand->frac;
        }
        N_sum=soil->NH4[l]+soil->NO3[l];
        if(N_sum>epsilon) /* immobilization of N */
        {
          n_immo=nimmo_fast*soil->c_shift[l][pftid].fast*N_sum/soildepth[l]*1e3/(k_N+N_sum/soildepth[l]*1e3);
          if(n_immo >0)
          {
            if(n_immo>N_sum)
              n_immo=N_sum;
            soil->pool[l].fast.nitrogen+=n_immo;
            soil->decomp_litter_pft[pftid].nitrogen+=n_immo;
            getoutput(&stand->cell->output,N_IMMO,config)+=n_immo*stand->frac;
            if(isagriculture(stand))
              getoutput(&stand->cell->output,NIMMOBILIZATION_AGR,config)+=n_immo*stand->frac;
            soil->NH4[l]-=n_immo*soil->NH4[l]/N_sum;
            soil->NO3[l]-=n_immo*soil->NO3[l]/N_sum;
#ifdef SAFE
            if(soil->NO3[l]<-epsilon)
              fail(NEGATIVE_SOIL_NO3_ERR,TRUE,TRUE,"Negative soil NO3=%g in layer %d in cell (%s) at immobilization in littersom()",soil->NO3[l],l,sprintcoord(line,&stand->cell->coord));
            if(soil->NH4[l]<-epsilon)
              fail(NEGATIVE_SOIL_NH4_ERR,TRUE,TRUE,"Negative soil NH4=%g in layer %d in cell (%s) at immobilization in littersom()",soil->NH4[l],l,sprintcoord(line,&stand->cell->coord));
#endif
          }
        }
        N_sum=soil->NH4[l]+soil->NO3[l];
        if(N_sum>epsilon)
        {
          n_immo=nimmo_slow*soil->c_shift[l][pftid].slow*N_sum/soildepth[l]*1e3/(k_N+N_sum/soildepth[l]*1e3);
          if(n_immo >0)
          {
            if(n_immo>N_sum)
              n_immo=N_sum;
            soil->pool[l].slow.nitrogen+=n_immo;
            soil->decomp_litter_pft[pftid].nitrogen+=n_immo;
            getoutput(&stand->cell->output,N_IMMO,config)+=n_immo*stand->frac;
            if(isagriculture(stand))
              getoutput(&stand->cell->output,NIMMOBILIZATION_AGR,config)+=n_immo*stand->frac;
            soil->NH4[l]-=n_immo*soil->NH4[l]/N_sum;
            soil->NO3[l]-=n_immo*soil->NO3[l]/N_sum;
#ifdef SAFE
            if(soil->NO3[l]<-epsilon)
              fail(NEGATIVE_SOIL_NO3_ERR,TRUE,TRUE,"Negative soil NO3=%g in layer %d in cell (%s) at immobilization in littersom()",soil->NO3[l],l,sprintcoord(line,&stand->cell->coord));
            if(soil->NH4[l]<-epsilon)
              fail(NEGATIVE_SOIL_NH4_ERR,TRUE,TRUE,"Negative soil NH4=%g in layer %d in cell (%s) at immobilization in littersom()",soil->NH4[l],l,sprintcoord(line,&stand->cell->coord));
#endif
          }
        }
      } /* of forrootlayer */
      /*sum for equilsom-routine*/
    }   /*end soil->litter.n*/
    /*sum for equilsom-routine*/
    soil->decomp_litter_mean.carbon+=decom_litter.carbon;
    soil->decomp_litter_mean.nitrogen+=decom_litter.nitrogen;
  } /* end of gtemp_soil[0]>0 */

  /* NO3 and N2O from nitrification */
  forrootsoillayer(l)
  {
    fac_wfps = f_wfps(soil,l);
    fac_temp = f_temp(soil->temp[l]);
    F_NO3=param.k_max*soil->NH4[l]*fac_temp*fac_wfps*f_ph(stand->cell->soilph);
    if(F_NO3>soil->NH4[l])
      F_NO3=soil->NH4[l];
    F_N2O=param.k_2*F_NO3;
    soil->NO3[l]+=F_NO3*(1-param.k_2);
#ifdef SAFE
    if(soil->NO3[l]<-epsilon)
      fail(NEGATIVE_SOIL_NO3_ERR,TRUE,TRUE,"littersom: Negative soil NO3=%g in layer %d in cell (%s)",
           soil->NO3[l],l,sprintcoord(line,&stand->cell->coord));
#endif

    soil->NH4[l]-=F_NO3;
#ifdef SAFE
    if(soil->NH4[l]<-epsilon)
      fail(NEGATIVE_SOIL_NH4_ERR,TRUE,TRUE,"Negative soil NH4=%g in layer %d in cell (%s)",
           soil->NH4[l],l,sprintcoord(line,&stand->cell->coord));
#endif
    flux.nitrogen += F_N2O;
    /* F_N2O is given back for output */
    if(stand->type->landusetype==AGRICULTURE)
    {
      foreachpft(pft,p,&stand->pftlist)
      {
        crop=pft->data;
        if(crop->sh!=NULL)
        {
          crop->sh->n2o_nitsum+=F_N2O;
          crop->sh->c_emissum+=decom_litter.carbon*param.atmfrac+soil_cflux;
        }
        else
        {
          data=stand->data;
          getoutputindex(&stand->cell->output,CFT_N2O_NIT,pft->par->id-npft+data->irrigation*ncft,config)+=F_N2O;
          getoutputindex(&stand->cell->output,CFT_C_EMIS,pft->par->id-npft+data->irrigation*ncft,config)+=decom_litter.carbon *param.atmfrac+soil_cflux;
        }
      }
    }
  }
#ifdef MICRO_HEATING
  soil->litter.decomC=decom_litter.carbon*param.atmfrac; /*only for mircobiological heating*/
#endif
  soil->count++;
  getoutput(&stand->cell->output,RH_LITTER,config)+=decom_litter.carbon*param.atmfrac*stand->frac;
  if(isagriculture(stand))
    getoutput(&stand->cell->output,RH_AGR,config) += (decom_litter.carbon*param.atmfrac+soil_cflux)*stand->frac;
  flux.carbon=decom_litter.carbon*param.atmfrac+soil_cflux;
  return flux;
} /* of 'littersom_nomethane' */

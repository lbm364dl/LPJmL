/**************************************************************************************/
/**                                                                                \n**/
/**                  t  r  a  c  e  _  d  a  y  .  c                               \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Checksums of a few global quantities, for finding where two runs of the    \n**/
/**     same configuration stop agreeing.                                          \n**/
/**                                                                                \n**/
/**     Enabled with -DTRACE_DAILY. Each quantity is printed in hexadecimal float  \n**/
/**     so the comparison is exact rather than to some number of digits. The       \n**/
/**     daily checksums are taken either side of the river routing and the annual  \n**/
/**     ones either side of the annual update, so a divergence can be attributed   \n**/
/**     to a part of the year rather than only to the year.                        \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"

#ifdef TRACE_DAILY

#define NTRACE 6

static void emit(const char *what,int year,int day,const char *tag,
                 Real sum[],const Config *config)
{
  Real all[NTRACE];
  int i;
#ifdef USE_MPI
  MPI_Reduce(sum,all,NTRACE,MPI_DOUBLE,MPI_SUM,0,config->comm);
#else
  for(i=0;i<NTRACE;i++)
    all[i]=sum[i];
#endif
  if(isroot(*config))
  {
    printf("%s %d %3d %-8s",what,year,day,tag);
    for(i=0;i<NTRACE;i++)
      printf(" %a",all[i]);
    putchar('\n');
    fflush(stdout);
  }
} /* of 'emit' */

void trace_day(const Cell grid[],     /**< cell array */
               int day,               /**< day of year */
               int year,              /**< simulation year */
               const char *tag,       /**< where in the day this was taken */
               const Config *config   /**< LPJmL configuration */
              )
{
  Real sum[NTRACE];
  int cell,i;
  for(i=0;i<NTRACE;i++)
    sum[i]=0;
  for(cell=0;cell<config->ngridcell;cell++)
    if(!grid[cell].skip)
    {
      sum[0]+=grid[cell].discharge.drunoff;      /* local runoff, before routing */
      sum[1]+=grid[cell].discharge.dmass_river;  /* water in transit */
      sum[2]+=grid[cell].discharge.dmass_lake;
      sum[3]+=grid[cell].discharge.dfout;        /* routed outflow */
      sum[4]+=grid[cell].balance.aprec;          /* climate, hence the RNG */
      sum[5]+=grid[cell].balance.anpp;
    }
  emit("TRACE",year,day,tag,sum,config);
  if(day==1 && !strcmp(tag,"cells"))
    /* one line per cell, so two runs can be diffed down to the cell that
       first disagrees */
    for(cell=0;cell<config->ngridcell;cell++)
      if(!grid[cell].skip)
        printf("CELL %6d %a %a %a\n",cell+config->startgrid,
               grid[cell].discharge.drunoff,grid[cell].balance.anpp,
               grid[cell].balance.aevap);
} /* of 'trace_day' */

void trace_annual(const Cell grid[],   /**< cell array */
                  int year,            /**< simulation year */
                  const char *tag,     /**< where in the year this was taken */
                  const Config *config /**< LPJmL configuration */
                 )
{
  Real sum[NTRACE];
  const Stand *stand;
  const Pft *pft;
  int cell,i,s,p;
  for(i=0;i<NTRACE;i++)
    sum[i]=0;
  for(cell=0;cell<config->ngridcell;cell++)
    if(!grid[cell].skip)
    {
      foreachstand(stand,s,grid[cell].standlist)
      {
        sum[0]+=standstocks(stand).carbon*stand->frac;
        sum[1]+=standstocks(stand).nitrogen*stand->frac;
        sum[2]+=stand->frac;
        foreachpft(pft,p,&stand->pftlist)
          sum[3]+=pft->nind;
        sum[4]+=getnpft(&stand->pftlist);
      }
      /* soil water, ice and heat, which the stock sums above do not cover */
      foreachstand(stand,s,grid[cell].standlist)
      {
        for(i=0;i<NSOILLAYER;i++)
        {
          sum[5]+=stand->soil.w[i]+stand->soil.w_fw[i];
          sum[5]+=stand->soil.ice_depth[i]+stand->soil.ice_fw[i];
          sum[5]+=stand->soil.ice_pwp[i]+stand->soil.freeze_depth[i];
          sum[5]+=stand->soil.temp[i]+stand->soil.k_dry[i];
          sum[5]+=stand->soil.perc_energy[i];
        }
        for(i=0;i<NHEATGRIDP;i++)
          sum[5]+=stand->soil.enth[i];
      }
    }
  emit("TRACEY",year,0,tag,sum,config);
} /* of 'trace_annual' */

#endif

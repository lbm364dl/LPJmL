/**************************************************************************************/
/**                                                                                \n**/
/**                   p  e  t  p  a  r  .  c                                       \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Function calculates daily photosynthetically active                        \n**/
/**     radiation flux, daylength and daily potential evapotranspiration           \n**/
/**     given temperature, longwave radiation, shortwave radiation, and            \n**/
/**     latitude                                                                   \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include <stdio.h>
#include <math.h>
#include "types.h"
#include "date.h"
#include "hash.h"
#include "bstruct.h"
#include "numeric.h"
#include "units.h"

#define dayseconds 86400 /* seconds per day */

#define gamma_t(temp)  (65.05+temp*0.064)
#define lambda(temp) (2.495e6-temp*2380)

void petpar_cell(Petpar *p,      /**< [out] the parts that do not vary by stand */
                 Real lat,        /**< latitude (deg) */
                 int day,         /**< day (1..365) */
                 Real temp,       /**< temperature (deg C) */
                 Real lw,         /**< longwave net/downward flux (W m-2) */
                 Real swdown,     /**< shortwave downward flux (W m-2) */
                 Bool islwdown    /**< LW radiation is downward (TRUE/FALSE) */
                )
{
  /*
   * Everything here depends on the cell and the day, not on the stand: the
   * solar geometry, the longwave correction and the slope of the saturation
   * vapour pressure curve.  Albedo is the only per-stand input, and it enters
   * only through swnet in petpar_stand().  Called once per stand this was two
   * sines, two cosines, an arccosine, a pow and an exp per stand per day.
   */
  Real delta,u,v,hh,s;
  delta=deg2rad(-23.4*cos(2*M_PI*(day+10.0)/NDAYYEAR));
  u=sin(deg2rad(lat))*sin(delta);
  v=cos(deg2rad(lat))*cos(delta);
  if(u>=v)
  {
    p->daylength=24;
  }
  else if(u<=-v)
    p->daylength=0;
  else
  {
    hh=acos(-u/v);
    p->daylength=24*hh*M_1_PI;
  }
  /* *par=dayseconds*swnet/2;*/
  p->par=dayseconds*swdown/2; /* MPAR based on SWdown instead of SWnet because albedo will be removed with PFT-dependent leaf albedo in water_stressed */

  if(islwdown)
    lw-=sigma_B*pow(degCtoK(temp),4);
  p->lw=lw;
  s=2.503e6*exp(17.269*temp/(237.3+temp))/
          ((237.3+temp)*(237.3+temp));
  p->eeqfac=s/(s+gamma_t(temp))/lambda(temp);
  p->swdown=swdown;
} /* of 'petpar_cell' */

void petpar_stand(Real *daylength, /**< daylength (h) */
                  Real *par,       /**< photosynthetic active radiation flux */
                  Real *eeq,       /**< equilibrium evapotranspiration */
                  const Petpar *p, /**< values shared by every stand of the cell */
                  Real beta        /**< Albedo */
                 )
{
  Real swnet;
  *daylength=p->daylength;
  *par=p->par;
  swnet=(1-beta)*p->swdown;  /* shortwave net flux, downward positive (W m-2) */
  /* calculation of PET respects only longwave flux during daylight  */
  *eeq=dayseconds*p->eeqfac*(swnet+p->lw*(*daylength)/24);
  if (*eeq<0) *eeq=0;
} /* of 'petpar_stand' */

void petpar(Real *daylength, /**< daylength (h) */
            Real *par,       /**< photosynthetic active radiation flux */
            Real *eeq,       /**< equilibrium evapotranspiration */
            Real lat,        /**< latitude (deg) */
            int day,         /**< day (1..365) */
            Real temp,       /**< temperature (deg C) */
            Real lw,         /**< longwave net/downward flux (W m-2) */
            Real swdown,     /**< shortwave downward flux (W m-2) */
            Bool islwdown,   /**< LW radiation is downward (TRUE/FALSE) */
            Real beta        /**< Albedo */
           )
{
  /* the whole calculation for a single stand, as before */
  Petpar p;
  petpar_cell(&p,lat,day,temp,lw,swdown,islwdown);
  petpar_stand(daylength,par,eeq,&p,beta);
} /* of 'petpar' */

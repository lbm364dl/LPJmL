/**************************************************************************************/
/**                                                                                \n**/
/**                 f  p  a  r  _  c  r  o  p  .  c                                \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
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

Real fpar_crop(const Pft *pft /**< pointer to PFT data */
              )               /** \return fractional PAR (0..1) */
{
  const Pftcrop *crop;
  crop=pft->data;
  if(pft->par->ismaize)
    return min(1,max(0,0.2558*max(0.01,crop->lai-crop->lai_nppdeficit)-0.0024))*(1-pft->snowcover);
  else
  {
    /* the same canopy on the same day for every stand that grows it */
    static Real lx=0,lv=0;
    static Bool have=FALSE;
    Real x=-pft->par->lightextcoeff*max(0,(crop->lai-crop->lai_nppdeficit));
    if(!have || x!=lx)
    {
      lx=x;
      lv=1-exp(x);
      have=TRUE;
    }
    return lv*(1-pft->snowcover);
  }
} /* of 'fpar_crop' */

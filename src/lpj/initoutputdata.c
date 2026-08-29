/**************************************************************************************/
/**                                                                                \n**/
/**                i  n  i  t  o  u  t  p  u  t  d  a  t  a  .  c                  \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Function initializes output data to zero                                   \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"

#define isinit(index) (isinit2(index,timestep,year,config))

static Bool isinit2(int index,int timestep,int year,const Config *config)
{
  if(timestep==ANNUAL && config->outnames[index].timestep>0)
    return ((year-config->outputyear) %  config->outnames[index].timestep==0);
  else
    return config->outnames[index].timestep==timestep;
} /* of 'isinit2' */

/* The scan below used to run the whole 397-entry output table for every cell on
   every day, to find the handful of entries that the timestep initialises: 435
   million tests in a three-year 1000-cell run, for five outputs.

   Which entries can be initialised at all is fixed once the configuration is
   read.  For DAILY and MONTHLY isinit2() reduces to outnames[index].timestep ==
   timestep, so those lists are exact; for ANNUAL it is the modulo test for
   entries with a positive timestep and always true for entries equal to ANNUAL,
   so the union of the two is the candidate list and isinit() still decides.
   Entries off the list cannot be true, so the same values are written. */

static int initidx[3][NOUT];
static int ninitidx[3];
static const Config *listconfig=NULL;

static int slotof(int timestep)
{
  switch(timestep)
  {
    case DAILY:   return 0;
    case MONTHLY: return 1;
    case ANNUAL:  return 2;
    default:      return -1;
  }
}

static void buildlists(const Config *config)
{
  static const int steps[3]={DAILY,MONTHLY,ANNUAL};
  int slot,index;
  for(slot=0;slot<3;slot++)
  {
    ninitidx[slot]=0;
    for(index=FPC;index<NOUT;index++)
      if(steps[slot]==ANNUAL
         ? (config->outnames[index].timestep>0 ||
            config->outnames[index].timestep==ANNUAL)
         : config->outnames[index].timestep==steps[slot])
        initidx[slot][ninitidx[slot]++]=index;
  }
  listconfig=config;
}

void initoutputdata(Output *output,      /**< output data */
                    int timestep,        /**< time step (ANNUAL, MONTHLY, DAILY) */
                    int year,            /**< simulation year */
                    const Config *config /**< LPJ configuration */
                   )
{
  int i,index,slot,k;
  slot=slotof(timestep);
  if(slot<0)
  {
    /* an unexpected timestep: fall back on the exhaustive scan */
    for(index=FPC;index<NOUT;index++)
      if(isinit(index))
        for(i=0;i<config->outputsize[index];i++)
          output->data[config->outputmap[index]+i]=0;
  }
  else
  {
    if(listconfig!=config)
      buildlists(config);
    /* set output data to zero */
    for(k=0;k<ninitidx[slot];k++)
    {
      index=initidx[slot][k];
      if(isinit(index))
        for(i=0;i<config->outputsize[index];i++)
          output->data[config->outputmap[index]+i]=0;
    }
  }
  if(isinit(PFT_GCGP))
    for(i=0;i<config->outputsize[PFT_GCGP];i++)
      output->data[config->outputmap[PFT_GCGP_COUNT]+i]=0;
  if(isinit(CFT_SWC))
    for(i=0;i<config->outputsize[CFT_SWC];i++)
      output->data[config->outputmap[NDAY_MONTH]+i]=0;
  /* set specific output data to one */
  if(isinit(DECAY_WOOD_AGR))
    output->data[config->outputmap[DECAY_WOOD_AGR]]=1;
  if(isinit(DECAY_WOOD_NV))
    output->data[config->outputmap[DECAY_WOOD_NV]]=1;
  if(isinit(DECAY_LEAF_AGR))
    output->data[config->outputmap[DECAY_LEAF_AGR]]=1;
  if(isinit(DECAY_LEAF_NV))
    output->data[config->outputmap[DECAY_LEAF_NV]]=1;
} /* of 'initoutputdata' */

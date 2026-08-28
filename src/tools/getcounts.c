/**************************************************************************************/
/**                                                                                \n**/
/**                  g  e  t  c  o  u  n  t  s  .  c                               \n**/
/**                                                                                \n**/
/**     Function computes counts and offsets used by MPI_Gatherv/                  \n**/
/**     MPI_Scatterv                                                               \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include <stdio.h>
#include "types.h"

void getcounts(int counts[],      /**< number of items for each task, filled in */
               int offsets[],     /**< item offsets for each task, filled in */
               const int cells[], /**< number of cells of each task */
               int n,             /**< number of items per cell */
               int tasks          /**< number of tasks */
              )
{
  /*
   * This used to recompute the uniform split from the total cell count, which
   * agreed with the model only while every task held the same number of cells.
   * It now takes the actual split, config->cellcounts, so that the gather
   * layout follows the decomposition wherever that comes from.
   */
  int i;
  for(i=0;i<tasks;i++)
    counts[i]=cells[i]*n;
  offsets[0]=0;
  for(i=1;i<tasks;i++)
    offsets[i]=offsets[i-1]+counts[i-1];
} /* of 'getcounts' */

/**************************************************************************************/
/**                                                                                \n**/
/**        p  n  e  t  _  s  o  r  t  c  o  n  n  e  c  t  .  c                    \n**/
/**                                                                                \n**/
/**     MPI-parallelization of networks                                            \n**/
/**                                                                                \n**/
/**     Function sorts each connection list by the global index of its source,     \n**/
/**     which makes the order the caller walks a node's incoming connections a     \n**/
/**     property of the network alone.                                             \n**/
/**                                                                                \n**/
/**     Without it the order comes out of pnet_reverse(): the reversed lists are    \n**/
/**     assembled from an MPI_Alltoallv, so a node's sources arrive grouped by      \n**/
/**     the task that owned them and the order changes whenever the                 \n**/
/**     decomposition does. Summing over them then sums in a different order,       \n**/
/**     and in LPJmL that is enough to move river discharge in the last bits and    \n**/
/**     flip threshold decisions -- irrigation, sowing, reservoir filling -- in     \n**/
/**     the cells that sit near one.                                                \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include "types.h"
#include "pnet.h"

static int compare(const void *a,const void *b)
{
  /* compare function to sort integers in ascending order, used by qsort() */
  return *((const int *)a)-*((const int *)b);
} /* of 'compare' */

int pnet_sortconnect(Pnet *pnet /**< Pointer to Pnet structure */
                    )           /** \return error code        */
{
  int i;
  if(pnet==NULL)
    return PNET_NULL_PTR_ERR;
  pnet_foreach(pnet,i)
    if(pnet->connect[i].n>1)
      qsort(pnet->connect[i].index,pnet->connect[i].n,sizeof(int),compare);
  return PNET_OK;
} /* of 'pnet_sortconnect' */

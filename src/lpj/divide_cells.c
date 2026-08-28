/**************************************************************************************/
/**                                                                                \n**/
/**              d  i  v  i  d  e  _  c  e  l  l  s  .  c                          \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Distribution of the cell grid over the parallel tasks.                     \n**/
/**                                                                                \n**/
/**     LPJmL hands every task the same *number* of cells. The grid runs from      \n**/
/**     north to south and the cost of a cell spans more than an order of          \n**/
/**     magnitude -- an ice cell carries one stand, a cropland cell carries one    \n**/
/**     per crop functional type and irrigation regime -- so an equal count is     \n**/
/**     not an equal share of the work. With river routing the tasks meet at a     \n**/
/**     barrier every simulated day, so the run proceeds at the pace of the        \n**/
/**     most heavily loaded task and the rest wait.                                \n**/
/**                                                                                \n**/
/**     Given a cell cost file the cells are split by cost instead, optionally     \n**/
/**     scaled by a per-task weight so a machine whose cores differ in speed can   \n**/
/**     be loaded in proportion to core speed. The blocks stay contiguous and in   \n**/
/**     rank order, which is what the routing network, the output gather and the   \n**/
/**     restart file all assume.                                                   \n**/
/**                                                                                \n**/
/**     Results do not depend on the decomposition. The inflow summation in        \n**/
/**     drain() walks the connection list of a cell in an order fixed by the       \n**/
/**     routing network, and pnet_setup() remaps only the values of those          \n**/
/**     indices, never their order. A stale or wrong cost file therefore costs     \n**/
/**     performance and nothing else.                                             \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"

#define CELLCOST_HEADER "LPJCOST"
#define CELLCOST_VERSION 1
#define MINCOST 1e-9 /* cost floor, so that no block can be empty of cost */

Real *cellcost=NULL; /* per-cell cost measured during the run, or NULL */

static void bounds_equal(int n,int ntask,int bound[])
{
  /* the historical split: the same number of cells for every task */
  int r;
  bound[0]=0;
  for(r=0;r<ntask;r++)
    bound[r+1]=bound[r]+n/ntask+((n % ntask>r) ? 1 : 0);
} /* of 'bounds_equal' */

static void bounds_cost(int n,int ntask,const Real cost[],const Real weight[],
                        int bound[])
{
  /* place the ntask-1 interior boundaries on the cumulative cost curve */
  int i,r;
  Real *prefix,*wsum,total,target;
  prefix=newvec(Real,n);
  wsum=newvec(Real,ntask+1);
  if(prefix==NULL || wsum==NULL)
  {
    /* running out of memory over a performance hint is not worth failing for */
    free(prefix);
    free(wsum);
    bounds_equal(n,ntask,bound);
    return;
  }
  total=0;
  for(i=0;i<n;i++)
  {
    /* without a cost file every cell counts the same, which still lets the
       task weights alone shift work towards the faster cores */
    total+=(cost==NULL) ? 1.0 : ((cost[i]>MINCOST) ? cost[i] : MINCOST);
    prefix[i]=total;
  }
  wsum[0]=0;
  for(r=0;r<ntask;r++)
    wsum[r+1]=wsum[r]+((weight!=NULL && weight[r]>0) ? weight[r] : 1.0);
  bound[0]=0;
  bound[ntask]=n;
  i=0;
  for(r=1;r<ntask;r++)
  {
    target=total*wsum[r]/wsum[ntask];
    while(i<n && prefix[i]<target)
      i++;
    bound[r]=i+1;
    /* every task must end up with at least one cell */
    if(bound[r]<bound[r-1]+1)
      bound[r]=bound[r-1]+1;
    if(bound[r]>n-(ntask-r))
      bound[r]=n-(ntask-r);
    i=bound[r];
  }
  free(prefix);
  free(wsum);
} /* of 'bounds_cost' */

void divide_counts(int n,           /**< number of grid cells to distribute */
                   int ntask,       /**< total number of tasks */
                   const Real cost[],  /**< cost of cell i, or NULL for equal counts */
                   const Real weight[],/**< relative speed of each task, or NULL */
                   int counts[]     /**< number of cells of each task, filled in */
                  )
{
  /*
   * The one place the decomposition is decided.  Everything that needs to know
   * how many cells a task holds -- the output gather, the routing network, the
   * reservoir network -- has to ask here rather than recompute it, or it will
   * disagree with the model the moment the split stops being uniform.
   */
  int r,*bound;
  bound=newvec(int,ntask+1);
  if(bound==NULL)
  {
    for(r=0;r<ntask;r++)
      counts[r]=n/ntask+((n % ntask>r) ? 1 : 0);
    return;
  }
  if(cost==NULL && weight==NULL)
    bounds_equal(n,ntask,bound);
  else
    bounds_cost(n,ntask,cost,weight,bound);
  for(r=0;r<ntask;r++)
    counts[r]=bound[r+1]-bound[r];
  free(bound);
} /* of 'divide_counts' */

void divide_cells(int *start,       /**< index of first grid cell, set to local first cell */
                  int *end,         /**< index of last grid cell, set to local last cell */
                  int rank,         /**< my rank id */
                  int ntask,        /**< total number of tasks */
                  const Real cost[],  /**< cost of cell *start+i, or NULL for equal counts */
                  const Real weight[] /**< relative speed of each task, or NULL for equal */
                 )
{
  int r,n,lo,*counts;
  n=*end-*start+1;
  counts=newvec(int,ntask);
  if(counts==NULL)
    return;               /* leave the range alone rather than corrupt it */
  divide_counts(n,ntask,cost,weight,counts);
  lo=*start;
  for(r=0;r<rank;r++)
    lo+=counts[r];
  *end=lo+counts[rank]-1;
  *start=lo;
  free(counts);
} /* of 'divide_cells' */

Real *readcellcost(const char *filename, /**< name of cell cost file */
                   int firstcell,        /**< index of first cell wanted */
                   int n,                /**< number of cells wanted */
                   Verbosity verbose     /**< error output enabled */
                  )                      /** \return cost array or NULL on error */
{
  FILE *file;
  char id[sizeof(CELLCOST_HEADER)-1];
  int version,first,ncell;
  Real *cost;
  file=fopen(filename,"rb");
  if(file==NULL)
  {
    if(verbose)
      printfopenerr(filename);
    return NULL;
  }
  cost=NULL;
  if(fread(id,sizeof(id),1,file)!=1 || memcmp(id,CELLCOST_HEADER,sizeof(id)) ||
     fread(&version,sizeof(int),1,file)!=1 ||
     fread(&first,sizeof(int),1,file)!=1 ||
     fread(&ncell,sizeof(int),1,file)!=1)
  {
    if(verbose)
      fprintf(stderr,"WARNING041: Cannot read header of cell cost file '%s', "
                     "cells are distributed by count.\n",filename);
  }
  else if(version!=CELLCOST_VERSION)
  {
    if(verbose)
      fprintf(stderr,"WARNING041: Version %d of cell cost file '%s' is not %d, "
                     "cells are distributed by count.\n",
              version,filename,CELLCOST_VERSION);
  }
  else if(first>firstcell || first+ncell<firstcell+n)
  {
    if(verbose)
      fprintf(stderr,"WARNING041: Cell cost file '%s' covers cells %d-%d, "
                     "not %d-%d, cells are distributed by count.\n",
              filename,first,first+ncell-1,firstcell,firstcell+n-1);
  }
  else
  {
    cost=newvec(Real,n);
    if(cost==NULL)
      printallocerr("cellcost");
    else if(fseek(file,(long)(firstcell-first)*sizeof(Real),SEEK_CUR) ||
            fread(cost,sizeof(Real),n,file)!=(size_t)n)
    {
      if(verbose)
        fprintf(stderr,"WARNING041: Cannot read %d values from cell cost file '%s', "
                       "cells are distributed by count.\n",n,filename);
      free(cost);
      cost=NULL;
    }
  }
  fclose(file);
  return cost;
} /* of 'readcellcost' */

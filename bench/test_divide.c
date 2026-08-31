/**************************************************************************************/
/**                                                                                \n**/
/**            t  e  s  t  _  d  i  v  i  d  e  .  c                               \n**/
/**                                                                                \n**/
/**     Unit test for divide_cells().                                              \n**/
/**                                                                                \n**/
/**     Two things have to hold. With no cost file and no weights the split must    \n**/
/**     be the one LPJmL has always made, byte for byte, so that existing runs      \n**/
/**     are untouched. With either of them the partition must still be a            \n**/
/**     partition: contiguous, complete, in rank order, and with no empty block,    \n**/
/**     because the routing network, the output gather and the restart file all     \n**/
/**     assume exactly that.                                                        \n**/
/**                                                                                \n**/
/**     Build and run:                                                             \n**/
/**       make -C .. && bench/build_test_divide.sh                                  \n**/
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
#include <math.h>
#include "lpj.h"

static int nfail=0;

static void divide_old(int *start,int *end,int rank,int ntask)
{
  /* the split as it was before divide_cells(), kept here as the reference */
  int i,lo,hi,n;
  n=*end-*start+1;
  lo=*start; hi=*start+n/ntask-1;
  if(n % ntask) hi++;
  for(i=1;i<=rank;i++)
  {
    lo=hi+1; hi=lo+n/ntask-1;
    if(n % ntask>i) hi++;
  }
  *start=lo; *end=hi;
}

/* Check that the ntask blocks tile [first,first+n) exactly once, in order. */
static void check_partition(const char *what,int first,int n,int ntask,
                            const Real cost[],const Real weight[])
{
  int rank,s,e,prev;
  prev=first-1;
  for(rank=0;rank<ntask;rank++)
  {
    s=first; e=first+n-1;
    divide_cells(&s,&e,rank,ntask,cost,weight);
    if(s!=prev+1)
    {
      printf("  FAIL %s: rank %d starts at %d, previous block ended at %d\n",
             what,rank,s,prev);
      nfail++;
      return;
    }
    if(e<s)
    {
      printf("  FAIL %s: rank %d has an empty block [%d,%d]\n",what,rank,s,e);
      nfail++;
      return;
    }
    prev=e;
  }
  if(prev!=first+n-1)
  {
    printf("  FAIL %s: last block ends at %d, expected %d\n",what,prev,first+n-1);
    nfail++;
  }
}

int main(void)
{
  int n,ntask,rank,checked=0,i,f,ni,it;
  int ns,ne,os,oe;
  Real *cost,*w,*load,mean,hi,lo;
  unsigned short seed[3]={7,11,13};
  int firsts[]={0,1,17,30000};
  int sizes[]={30,31,64,1000,58795};
  int tasks[]={1,2,3,7,8,24,30,64};

  puts("1. no cost file, no weights: must reproduce the historical split");
  for(f=0;f<4;f++)
    for(ni=0;ni<5;ni++)
      for(i=0;i<8;i++)
      {
        n=sizes[ni]; ntask=tasks[i];
        if(ntask>n) continue;
        for(rank=0;rank<ntask;rank++)
        {
          ns=os=firsts[f]; ne=oe=firsts[f]+n-1;
          divide_cells(&ns,&ne,rank,ntask,NULL,NULL);
          divide_old(&os,&oe,rank,ntask);
          checked++;
          if(ns!=os || ne!=oe)
          {
            printf("  FAIL n=%d ntask=%d rank=%d: new [%d,%d] old [%d,%d]\n",
                   n,ntask,rank,ns,ne,os,oe);
            nfail++;
          }
        }
      }
  printf("   %d partitions checked\n",checked);

  puts("2. degenerate cost distributions must still give a valid partition");
  n=1000;
  cost=newvec(Real,n);
  w=newvec(Real,64);
  load=newvec(Real,64);
  for(i=0;i<64;i++)
    w[i]=(i%3==0) ? 1.0 : 0.4;

  for(it=0;it<7;it++)
  {
    const char *what;
    switch(it)
    {
      case 0: what="all cells free";        for(i=0;i<n;i++) cost[i]=0;                       break;
      case 1: what="all cost in cell 0";    for(i=0;i<n;i++) cost[i]=0; cost[0]=1e9;          break;
      case 2: what="all cost in last cell"; for(i=0;i<n;i++) cost[i]=0; cost[n-1]=1e9;        break;
      case 3: what="negative cost";         for(i=0;i<n;i++) cost[i]=-1;                      break;
      case 4: what="one huge, rest tiny";   for(i=0;i<n;i++) cost[i]=1e-12; cost[n/2]=1e12;   break;
      case 5: what="uniform";               for(i=0;i<n;i++) cost[i]=1;                       break;
      default:what="random";                for(i=0;i<n;i++) cost[i]=erand48(seed)*100;       break;
    }
    for(i=0;i<8;i++)
    {
      ntask=tasks[i];
      if(ntask>n) continue;
      check_partition(what,0,n,ntask,cost,NULL);
      check_partition(what,30000,n,ntask,cost,w);
    }
    printf("   %-22s ok\n",what);
  }

  puts("3. weights alone, with no cost file, must shift work to the fast tasks");
  ntask=24;
  n=58795;
  for(i=0;i<ntask;i++)
    w[i]=(i<8) ? 1.0 : 0.45;
  mean=0;
  for(rank=0;rank<ntask;rank++)
  {
    ns=0; ne=n-1;
    divide_cells(&ns,&ne,rank,ntask,NULL,w);
    load[rank]=(ne-ns+1)/w[rank];   /* time = cells / speed */
    mean+=load[rank]/ntask;
  }
  lo=hi=load[0];
  for(rank=1;rank<ntask;rank++)
  {
    if(load[rank]<lo) lo=load[rank];
    if(load[rank]>hi) hi=load[rank];
  }
  printf("   slowest/mean = %.4f (want ~1)\n",hi/mean);
  if(hi/mean>1.01)
  {
    puts("  FAIL: weights did not balance the load");
    nfail++;
  }

  puts("4. a real cost distribution must end up balanced");
  free(cost);
  cost=newvec(Real,n);   /* the full grid, this time */
  for(i=0;i<n;i++)  /* cheap in the north, expensive in the middle, as the grid is */
    cost[i]=0.2+9.0*exp(-((Real)i/n-0.55)*((Real)i/n-0.55)/0.02);
  for(ni=0;ni<2;ni++)
  {
    ntask=(ni==0) ? 24 : 30;
    mean=0;
    for(rank=0;rank<ntask;rank++)
    {
      ns=0; ne=n-1;
      divide_cells(&ns,&ne,rank,ntask,cost,NULL);
      load[rank]=0;
      for(i=ns;i<=ne;i++) load[rank]+=cost[i];
      mean+=load[rank]/ntask;
    }
    lo=hi=load[0];
    for(rank=1;rank<ntask;rank++)
    {
      if(load[rank]<lo) lo=load[rank];
      if(load[rank]>hi) hi=load[rank];
    }
    /* what the equal-count split would have given, for comparison */
    {
      Real emean=0,ehi;
      Real *eload=newvec(Real,ntask);
      for(rank=0;rank<ntask;rank++)
      {
        os=0; oe=n-1;
        divide_old(&os,&oe,rank,ntask);
        eload[rank]=0;
        for(i=os;i<=oe;i++) eload[rank]+=cost[i];
        emean+=eload[rank]/ntask;
      }
      ehi=eload[0];
      for(rank=1;rank<ntask;rank++)
        if(eload[rank]>ehi) ehi=eload[rank];
      printf("   ntask=%2d: equal count %.3f -> by cost %.3f  (%.2fx)\n",
             ntask,ehi/emean,hi/mean,ehi/hi);
      free(eload);
    }
    if(hi/mean>1.01)
    {
      puts("  FAIL: cost split did not balance the load");
      nfail++;
    }
  }

  free(cost); free(w); free(load);
  puts(nfail ? "\nFAILED" : "\nOK");
  return nfail!=0;
}

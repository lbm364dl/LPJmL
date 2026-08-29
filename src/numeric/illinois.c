/**************************************************************************************/
/**                                                                                \n**/
/**                    i  l  l  i  n  o  i  s  .  c                                \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Root of a function by the Illinois method: regula falsi with the retained  \n**/
/**     endpoint's value halved, which keeps the bracket that plain regula falsi   \n**/
/**     keeps while removing the stalling that makes it converge as slowly as      \n**/
/**     bisection on one-sided functions.                                          \n**/
/**                                                                                \n**/
/**     Same signature as bisect(), and the same best-effort contract: if the      \n**/
/**     iteration limit is reached, the x with the smallest |f| seen is returned   \n**/
/**     rather than the last one tried.                                            \n**/
/**                                                                                \n**/
/**     water_stressed() solves for lambda 8.1 million times in a three-year run   \n**/
/**     and bisect() averages 13.6 evaluations doing it, a fifth of them running   \n**/
/**     to the iteration cap without reaching the tolerance at all.                \n**/
/**                                                                                \n**/
/** (C) Potsdam Institute for Climate Impact Research (PIK), see COPYRIGHT file    \n**/
/** authors, and contributors see AUTHORS file                                     \n**/
/** This file is part of LPJmL and licensed under GNU AGPL Version 3               \n**/
/** or later. See LICENSE file or go to http://www.gnu.org/licenses/               \n**/
/** Contact: https://github.com/PIK-LPJmL/LPJmL                                    \n**/
/**                                                                                \n**/
/**************************************************************************************/

#include "lpj.h"

Real illinois(Real (*fcn)(Real,void *), /**< function to find the zero of */
              Real xlow,                /**< lower end of the bracket */
              Real xhigh,               /**< upper end of the bracket */
              void *data,               /**< passed through to fcn */
              Real xacc,                /**< accuracy in x, 0 to ignore */
              Real yacc,                /**< accuracy in y */
              int maxit,                /**< maximum number of iterations */
              int *it                   /**< iterations performed */
             )                          /** \return position of zero of function */
{
  int i;
  Real flow,fhigh,x,fx,ymin,xmin;
  flow=(*fcn)(xlow,data);
  fhigh=(*fcn)(xhigh,data);
  ymin=fabs(flow);
  xmin=xlow;
  if(fabs(fhigh)<ymin)
  {
    ymin=fabs(fhigh);
    xmin=xhigh;
  }
  if(ymin<yacc)
  {
    *it=0;
    return xmin;
  }
  for(i=0;i<maxit;i++)
  {
    /* The secant step.  It is taken only when it lands strictly inside the
       bracket; otherwise -- equal or non-finite endpoint values, or a bracket
       that holds no sign change and so lets the secant point outside -- the
       midpoint is used, which is what bisect() would have done. */
    x=(xlow+xhigh)*0.5;
    if(fhigh!=flow)
    {
      Real xs;
      xs=xhigh-fhigh*(xhigh-xlow)/(fhigh-flow);
      if((xs>xlow && xs<xhigh) || (xs>xhigh && xs<xlow))
        x=xs;
    }
    fx=(*fcn)(x,data);
    if(fabs(fx)<ymin)
    {
      ymin=fabs(fx);
      xmin=x;
    }
    if(fabs(fx)<yacc || (xacc>0 && fabs(xhigh-xlow)<xacc))
    {
      *it=i+1;
      return x;
    }
    if(fx*fhigh<0)   /* the root is between x and the old xhigh */
    {
      xlow=xhigh;
      flow=fhigh;
    }
    else
      flow*=0.5;     /* the Illinois modification */
    xhigh=x;
    fhigh=fx;
  }
  *it=i;
  return xmin;
} /* of 'illinois' */

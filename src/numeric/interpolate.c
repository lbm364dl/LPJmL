/**************************************************************************************/
/**                                                                                \n**/
/**                  i  n  t  e  r  p  o  l  a  t  e  .  c                         \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Interpolation function in order to get daily values from monthly           \n**/
/**     data                                                                       \n**/
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
#include "date.h"

Real interpolate(const MReal mval, /**< monthly values to interpolate        */
                 int month,        /**< month (0..11)                        */
                 int dm            /**< day of month (0..ndaymonth[month]-1) */
                )                  /** \return interpolated value            */
{
  // fprintf(stderr, "monthhhh %f\n", month);
  int nextmonth;
  if(dm>=ndaymonth[month]/2)
  {
    dm-=ndaymonth[month]/2;
    nextmonth=(month<NMONTH-1) ? month+1 : 0;
  }
  else
  {
    nextmonth=month;
    if(month==0)
      month=NMONTH-1;
    else
      month--; 
    dm+=(ndaymonth[month]+1)/2; 
  }
  // fprintf(stderr, "nextmonth %f\n", nextmonth);
#ifdef DEBUG3
  printf("mval[%d]=%g,mval[%d]=%g,dm=%d,diffday=%g,%d\n",month,mval[month],
         nextmonth,mval[nextmonth],dm,diffday[month],ndaymonth[month]);
  printf("res=%g\n",mval[month]+dm*(mval[nextmonth]-mval[month])*diffday[month]);
#endif
  // fprintf(stderr, "diffday %f\n", diffday[month]);
  // fprintf(stderr, "mvals ");
  // for (int i = 0; i < NMONTH; ++i) {
  //   fprintf(stderr, "%f ", mval[i]);
  // }
  // fprintf(stderr, "\n");
  // fprintf(stderr, "mymonths %i %i %f %f\n", nextmonth, month, mval[nextmonth], mval[month]);
  // fprintf(stderr, "diff %f\n", mval[nextmonth]-mval[month]);
  // fprintf(stderr, "product %f\n", dm*(mval[nextmonth]-mval[month])*diffday[month]);
  // fprintf(stderr, "all %f\n", mval[month]+dm*(mval[nextmonth]-mval[month])*diffday[month]);

  return mval[month]+dm*(mval[nextmonth]-mval[month])*diffday[month];
} /* of 'interpolate' */

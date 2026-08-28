/**************************************************************************************/
/**                                                                                \n**/
/**               i  t  e  r  a  t  e  y  e  a  r  .  c                            \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Function performs iteration over the cell grid for one year with           \n**/
/**     river routing enabled.                                                     \n**/
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
extern void trace_day(const Cell [],int,int,const char *,const Config *);
extern void trace_annual(const Cell [],int,const char *,const Config *);
#endif

Bool iterateyear(Outputfile *output,  /**< Output file data */
                 Cell grid[],         /**< cell array */
                 Input *input,        /**< input data */
                 Real co2,            /**< atmospheric CO2 (ppmv) */
                 Real *ch4,           /**< CH4 (gC) */
                 Real *pch4,          /**< CH4 concentration (ppmv) */
                 int npft,            /**< number of natural PFTs */
                 int ncft,            /**< number of crop PFTs */
                 int year,            /**< simulation year (AD) */
                 const Config *config /**< LPJ configuration */
                )                     /** \return TRUE on error */
{
  Dailyclimate daily;
  Bool intercrop;
  int month,dayofmonth,day;
  int cell;
  intercrop=getintercrop(input->landuse);
#ifdef TRACE_DAILY
  trace_annual(grid,year,"yearin",config);
#endif
  if(setupannual_grid(output,grid,input,year,npft,ncft,intercrop,config))
    return TRUE;
#ifdef TRACE_DAILY
  trace_annual(grid,year,"setup",config);
#endif
  day=1;
  foreachmonth(month)
  {
    initmonthly_grid(grid,month,year,input->climate,config);
    foreachdayofmonth(dayofmonth,month)
    {
      if(cellcost==NULL)
        for(cell=0;cell<config->ngridcell;cell++)
        {
          update_daily_cell(grid+cell,cell,&daily,co2,*pch4,input,day,dayofmonth,month,year,
                            npft,ncft,intercrop,config);
        }
      else
        /* measuring for a cell cost file: time each cell so that a later run
           can split the grid by cost instead of by cell count */
        for(cell=0;cell<config->ngridcell;cell++)
        {
          double tcell;
          tcell=mrun();
          update_daily_cell(grid+cell,cell,&daily,co2,*pch4,input,day,dayofmonth,month,year,
                            npft,ncft,intercrop,config);
          cellcost[cell]+=mrun()-tcell;
        }
#ifdef TRACE_DAILY
      trace_day(grid,day,year,"cells",config);
#endif
      updatedaily_grid(output,grid,input->extflow,day,month,year,npft,ncft,config);
#ifdef TRACE_DAILY
      trace_day(grid,day,year,"routed",config);
#endif
      day++;
    } /* of 'foreachdayofmonth */
    update_monthly_grid(output,grid,input->climate,month,year,npft,ncft,config);
  } /* of 'foreachmonth */
#ifdef TRACE_DAILY
  trace_annual(grid,year,"preann",config);
#endif
  updateannual_grid(output,grid,input->landcover,co2,ch4,pch4,year,npft,ncft,intercrop,daily.isdailytemp,config);
#ifdef TRACE_DAILY
  trace_annual(grid,year,"postann",config);
#endif
  return FALSE;
} /* of 'iterateyear' */

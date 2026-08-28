/**************************************************************************************/
/**                                                                                \n**/
/**        w  r  i  t  e  c  e  l  l  c  o  s  t  .  c                             \n**/
/**                                                                                \n**/
/**     C implementation of LPJmL                                                  \n**/
/**                                                                                \n**/
/**     Function collects the measured per-cell cost from all tasks and writes     \n**/
/**     the cell cost file that a later run can be balanced with.                  \n**/
/**                                                                                \n**/
/**     Kept apart from divide_cells.c because it is the only part that calls      \n**/
/**     MPI: the utility programs link liblpj.a with plain gcc, and the linker      \n**/
/**     pulls in whole object files.                                               \n**/
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

Bool writecellcost(const char *filename, /**< name of cell cost file to create */
                   const Real cost[],    /**< cost of the local cells */
                   const Config *config  /**< LPJmL configuration */
                  )                      /** \return TRUE on error */
{
  FILE *file;
  Real *all;
  int version,first;
  Bool rc;
#ifdef USE_MPI
  int *counts,*offsets,i;
#endif
  all=NULL;
#ifdef USE_MPI
  counts=newvec(int,config->ntask);
  offsets=newvec(int,config->ntask);
  if(counts==NULL || offsets==NULL)
  {
    free(counts);
    free(offsets);
    return TRUE;
  }
  MPI_Allgather((void *)&config->ngridcell,1,MPI_INT,counts,1,MPI_INT,config->comm);
  offsets[0]=0;
  for(i=1;i<config->ntask;i++)
    offsets[i]=offsets[i-1]+counts[i-1];
  if(isroot(*config))
  {
    all=newvec(Real,offsets[config->ntask-1]+counts[config->ntask-1]);
    if(all==NULL)
      printallocerr("cellcost");
  }
  MPI_Gatherv((void *)cost,config->ngridcell,MPI_DOUBLE,all,counts,offsets,
              MPI_DOUBLE,0,config->comm);
  free(counts);
  free(offsets);
#else
  all=(Real *)cost;
#endif
  rc=FALSE;
  if(isroot(*config) && all!=NULL)
  {
    file=fopen(filename,"wb");
    if(file==NULL)
    {
      printfcreateerr(filename);
      rc=TRUE;
    }
    else
    {
      version=CELLCOST_VERSION;
      first=config->firstgrid;
      if(fwrite(CELLCOST_HEADER,sizeof(CELLCOST_HEADER)-1,1,file)!=1 ||
         fwrite(&version,sizeof(int),1,file)!=1 ||
         fwrite(&first,sizeof(int),1,file)!=1 ||
         fwrite(&config->nall,sizeof(int),1,file)!=1 ||
         fwrite(all,sizeof(Real),config->nall,file)!=(size_t)config->nall)
      {
        fprintf(stderr,"ERROR153: Cannot write cell cost file '%s': %s\n",
                filename,strerror(errno));
        rc=TRUE;
      }
      else
        printf("Cell cost written to '%s'.\n",filename);
      fclose(file);
    }
  }
#ifdef USE_MPI
  free(all);
  MPI_Bcast(&rc,1,MPI_INT,0,config->comm);
#endif
  return rc;
} /* of 'writecellcost' */

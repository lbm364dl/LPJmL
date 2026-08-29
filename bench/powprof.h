/* Count how often each pow()/exp()/log()/log10()/sqrt() call site is handed the
   same arguments it
   was handed last time.  A site with a high repeat rate is one a single-entry
   memo would collapse, which is exactly the transformation that is exact:
   the same arguments return the same bits.
 *
 * Build with:
 *   OPTFLAGS= -O3 -fno-math-errno -DPOW_PROFILE -include bench/powprof.h
 *
 * Each translation unit keeps its own table and prints its own sites at exit,
 * so nothing has to be linked in.
 */
#ifndef POWPROF_H
#define POWPROF_H
#ifdef POW_PROFILE

#include <math.h>
#include <stdio.h>

#define POWPROF_MAXSITE 96

typedef struct
{
  const char *file;
  int line;
  const char *fn;
  long calls,repeats,hits8;
  double lastx,lasty;
  double cx[8],cy[8];      /* an 8-entry cache, to see past a 1-entry miss */
  int nc,next;
} PowSite;

static PowSite powprof_sites[POWPROF_MAXSITE];
static int powprof_n;

static double powprof_note(double x,double y,const char *file,int line,
                           const char *fn,double result)
{
  int i;
  for(i=0;i<powprof_n;i++)
    if(powprof_sites[i].line==line && powprof_sites[i].file==file)
      break;
  if(i==powprof_n)
  {
    if(powprof_n>=POWPROF_MAXSITE)
      return result;
    powprof_sites[i].file=file;
    powprof_sites[i].line=line;
    powprof_sites[i].fn=fn;
    powprof_sites[i].lastx=x-1;   /* guarantee a miss on the first call */
    powprof_n++;
  }
  powprof_sites[i].calls++;
  if(powprof_sites[i].lastx==x && powprof_sites[i].lasty==y)
    powprof_sites[i].repeats++;
  powprof_sites[i].lastx=x;
  powprof_sites[i].lasty=y;
  {
    int k,found=0;
    for(k=0;k<powprof_sites[i].nc;k++)
      if(powprof_sites[i].cx[k]==x && powprof_sites[i].cy[k]==y)
      {
        found=1;
        break;
      }
    if(found)
      powprof_sites[i].hits8++;
    else
    {
      if(powprof_sites[i].nc<8)
        powprof_sites[i].nc++;
      powprof_sites[i].cx[powprof_sites[i].next]=x;
      powprof_sites[i].cy[powprof_sites[i].next]=y;
      powprof_sites[i].next=(powprof_sites[i].next+1)%8;
    }
  }
  return result;
}

__attribute__((destructor)) static void powprof_dump(void)
{
  int i;
  for(i=0;i<powprof_n;i++)
    if(powprof_sites[i].calls>0)
      fprintf(stderr,"POWPROF %-34s %s:%d calls=%ld repeat=%.1f%% cache8=%.1f%%\n",
              powprof_sites[i].fn,powprof_sites[i].file,powprof_sites[i].line,
              powprof_sites[i].calls,
              100.0*powprof_sites[i].repeats/powprof_sites[i].calls,
              100.0*powprof_sites[i].hits8/powprof_sites[i].calls);
}

#define pow(x,y) powprof_note((x),(y),__FILE__,__LINE__,__func__,(pow)((x),(y)))
#define exp(x)   powprof_note((x),0.0,__FILE__,__LINE__,__func__,(exp)((x)))
#define log(x)   powprof_note((x),1.0,__FILE__,__LINE__,__func__,(log)((x)))
#define log10(x) powprof_note((x),2.0,__FILE__,__LINE__,__func__,(log10)((x)))
#define sqrt(x)  powprof_note((x),3.0,__FILE__,__LINE__,__func__,(sqrt)((x)))

#endif
#endif

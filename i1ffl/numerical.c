#include <stdlib.h>
#include <math.h>
#include "RngStream.h" /*replace ran1 with parallel RNG*/
#include "netsim.h"


#define MAXIT 100          /* maximum number of iterations for Newtown-Raphson */
#define RT_SAFE_EPSILON 0.01 /* Minimal change in root values between iteration for Newtown-Raphson.
                              * This parameter has unit of minute in this program.
                              */

float gasdev(RngStream RS)
{
//   static int iset=0;
//   static float gset;
   float fac,r,v1,v2;

//   if (iset == 0) {
      do {
         v1=2.0*RngStream_RandU01(RS)-1.0;
         v2=2.0*RngStream_RandU01(RS)-1.0;
         r=v1*v1+v2*v2;
      } while (r >= 1.0);
      fac=sqrt(-2.0*log(r)/r);
//      gset=v1*fac;
//      iset=1;
      return v2*fac;
//   } else {
//      iset=0;
//      return gset;
//   }
}

float expdev(RngStream RS)
{
  float dum;
  do dum=RngStream_RandU01(RS);
  while (dum == 0.0);
  return -log(dum);
}

/* 
 * Newton-Raphson root-finding method with bisection steps, out of
 * Numerical Recipes function bracketed by x1 and x2. Returns root
 * within accuracy +/-xacc funcd is function of interest, returning
 * both function value and first deriv.x  
 */
float rtsafe(void (*funcd)(float, int, float, float*, float*, float*, float*, float*, int), 
             int n_copies, float RHS, float *p_i, float *as_i, float *c_i, float x1, float x2, int start)
{
    int j;
    float df,dx,dxold,f,fh,fl;
    float temp,xh,xl,rts;

    (*funcd)(x1, n_copies, RHS,p_i,as_i, c_i, &fl, &df, start);
    (*funcd)(x2, n_copies, RHS,p_i,as_i, c_i, &fh, &df, start); /* note df isn't used here */
    if (fabs(fl) < 1e-9) return x1;
    if (fabs(fh) < 1e-9) return x2;
    
    if ((fl > 0.0 && fh > 0.0) || (fl <0.0 && fh < 0.0))
    {
//        if (verbose) fprintf(fperrors,"warning in rtsafe: root should be bracketed\n");
//        if (fabs(fl) < fabs(fh)) return x1; else return x2;
    }
    
    if (fl < 0.0) 
    {
        xl=x1;
        xh=x2;
    } 
    else 
    {
        xh=x1;
        xl=x2;
    }
    
    rts=0.5*(x1+x2);
    dxold=fabs(x2-x1);
    dx=dxold;    
    (*funcd)(rts, n_copies, RHS, p_i,as_i, c_i, &f, &df,start);

//    done = 0;
    
    for (j=1;j<=MAXIT;j++)
    {
        if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (fabs(2.0*f) > fabs(dxold*df))) 
        {
//            done = 1;// modified: otherwise this bisection can mean 2 identical function calls for j=1
            dxold=dx;
            dx=0.5*(xh-xl);
            rts=xl+dx;      
            if (xl == rts) return rts;
        } 
        else 
        {
            dxold=dx;
            dx=f/df;
            temp=rts;
            rts -= dx;
            if (temp == rts) return rts;
        }
        if (fabs(dx) < RT_SAFE_EPSILON) return rts;
        
//        if (rts==0.0)
//        {
//            if (x1<x2) rts = fminf(2.0*x1,(xl+xh)/2.0);
//            else rts = fminf(2.0*x2,(xl+xh)/2.0);
//            fprintf(fperrors,"warning: dt=0 reset to %g\n",rts);
//        }        
//        if (j>1 || done==0) 
            (*funcd)(rts, n_copies, RHS, p_i,as_i, c_i, &f, &df,start);
        if (f < 0.0) 
            xl=rts;
        else 
            xh=rts;   
    }
//    fprintf(fperrors,"error in rtsafe: too many iterations\n");
    return 0.0;
}

void find_max(float *data, int start, int end, float *max, int *id_max)
{
    int i;    
    *max=data[start];
    *id_max=start;
    for(i=0;i<end-start;i++)
    {
        if(*max<data[i])
        {
            *max=data[i];
            *id_max=start+i;
        }            
    }      
}

void find_x(float *data, int start, int end, float x, float *id_x, int find_x_near_start)
{
    int i;  
    
    if(find_x_near_start)        
    {
        i=start;
        while(i<end-1)
        {
            if(data[i]>=x)
            {
                if(data[i+1]<=x)
                {
                    *id_x=(x-data[i+1])/(data[i]-data[i+1])+(float)(i);
                    return;
                }
                else
                    i++;
            }
        }  
        /*should find x by now. Otherwise,do an extrapolation*/
        if(data[end]>=data[end-1]) //don't know when 
            *id_x=TIME_INFINITY;
        else
            *id_x=log(data[end]/x)/log(data[end-1]/data[end])+(float)(end);      
    }
    else
    {
        i=end;        
        while(i>start)
        {
            if(data[i]>=x)
            {
                if(data[i-1]<=x)
                {
                    *id_x=(x-data[i-1])/(data[i]-data[i-1])+(float)(i-1);
                    return;
                }
                else
                    i--;
            }            
        }       
    }
}
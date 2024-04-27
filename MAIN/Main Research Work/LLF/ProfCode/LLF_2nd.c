#include<stdio.h>
#include<math.h>
#include<time.h>

int main(){
FILE *fa;
fa=fopen("LLF_2nd_test5_200.txt","w+");
double start = clock();
int i,n,nx;
static double x[100001],dx,dt,d[100001],E[100001],En[100001],a[100001],p[100001],pn[100001],u[100001],un[100001],U1[100001],U1n[100001],U2[100001],U2n[100001],U3n[100001],U3[100001],G1l[100001],G2l[100001],G3l[100001],G1r[100001],G2r[100001],G3r[100001],lam_1[100001],l1max,cfl,l,t,an[100001],dn[100001],delU1[100001],delU2[100001],delU3[100001],interval,leftpoint,e[100001],en[100001],lam_maxph[100001],lam_maxmh[100001],G1ph[100001],G2ph[100001],G3ph[100001],G1mh[100001],G2mh[100001],G3mh[100001],gamma=1.4,x0,tn,f1[100001],f2[100001],tn1,H[100001],M[100001],s1l[100001],s2l[100001],s3l[100001],s1r[100001],s2r[100001],s3r[100001],U1l[100001],U2l[100001],U3l[100001],U1r[100001],U2r[100001],U3r[100001],G1[100001],G2[100001],G3[100001],sigma1[100001],sigma2[100001],sigma3[100001],sig1[100001],sig2[100001],sig3[100001],d1l[100001],d2l[100001],d3l[100001],d1r[100001],d2r[100001],d3r[100001], e1=0.000000000000001;

cfl = 0.1;
nx = 200;
leftpoint = 0.0;
interval = 1.0;
dx = interval/nx;
x0 = 0.8;
tn = 0.012;

/* ........................GRID GENERATION.............................*/
x[1] = leftpoint + (dx/2.0);

for(int i=2 ; i<=nx ; i++){
  x[i] = x[i-1] + dx;
}

/*........................Initial Conditions.............................*/
#include "initial.h"
initial_toro5();
/* ........................Initial Conditions.............................*/
/*
for(i=1;i<=nx;i++)
{
if( x[i] < 0.3)
{

d[i] =  1.0;

p[i] =  1.0;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.75;

E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

G1[i] = d[i]*u[i];
G2[i] = d[i]*u[i]*u[i] + p[i];
G3[i] = (0.5*d[i]*u[i]*u[i]*u[i]) + ((gamma/(gamma - 1.0))*p[i]*u[i]);

}

else
{

d[i] =  0.125;

p[i] =  0.1;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.0;

E[i] = (p[i]/((gamma-1.0)*d[i])) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

G1[i] = d[i]*u[i];
G2[i] = d[i]*u[i]*u[i] + p[i];
G3[i] = (0.5*d[i]*u[i]*u[i]*u[i]) + ((gamma/(gamma - 1.0))*p[i]*u[i]);


//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;


}

}
*/
/*....................Time iterations.....................*/
t = 0.0;

for(n=1;n<=5000000;n++){
         /*  .......................    PROCEDURE FOR CALCULATING VALUE OF lam_MAX..................... */
    for(i=1;i<=nx;i++)
    {
      lam_1[i]= fabs(u[i]) + a[i];
    }
      l1max = fabs(u[1]) + a[1];
           for(i=1;i<=nx;i++)
            {
               if(l1max<lam_1[i])
                   l1max = lam_1[i];
               else
                   l1max = l1max;
            }
dt  = (cfl*dx)/(l1max);
l = dt/dx ;

/*...................B.C's........................*/

      //u[0] = -u[1];                 // Reflective B.C
      u[0] = u[1];                  // Transmissive B.C
      d[0] = d[1];
      p[0] = p[1];
      a[0] = a[1];
      //u[nx+1] = -u[nx];              // Reflective B.C
      u[nx+1] = u[nx];                // Transmissive B.C
      d[nx+1] = d[nx];
      p[nx+1] = p[nx];
      a[nx+1] = a[nx];
      U1[0] = d[0];
      U2[0] = d[0]*u[0];
      U3[0] = (p[0]/(gamma-1.0)) + (0.5*d[0]*u[0]*u[0]);
      U1[nx+1] = d[nx+1];
      U2[nx+1] = d[nx+1]*u[nx+1];
      U3[nx+1] = (p[nx+1]/(gamma-1.0)) + (0.5*d[nx+1]*u[nx+1]*u[nx+1]);
      G1[0] = d[0]*u[0];
      G2[0] = d[0]*u[0]*u[0] + p[0];
      G3[0] = (0.5*d[0]*u[0]*u[0]*u[0]) + ((gamma/(gamma - 1.0))*p[0]*u[0]);
      G1[nx+1] = d[nx+1]*u[nx+1];
      G2[nx+1] = d[nx+1]*u[nx+1]*u[nx+1] + p[nx+1];
      G3[nx+1] = (0.5*d[nx+1]*u[nx+1]*u[nx+1]*u[nx+1]) + ((gamma/(gamma - 1.0))*p[nx+1]*u[nx+1]);





/*......................Minmod limiter............................*/

for(i=1;i<=nx;i++)
    {

       sigma1[i]  =  (U1[i+1] - U1[i]);

       sigma2[i]  =  (U2[i+1] - U2[i]);
       
       sigma3[i]  =  (U3[i+1] - U3[i]);

    }

/*..............................slope1 at left of j+1/2  interface ...........................*/

for(i=3;i<=nx-2;i++)
    {
            if(  (fabs(sigma1[i-1]) < fabs(sigma1[i])) && (  (sigma1[i-1])*(sigma1[i]) > 0.0 ) )
               {

                  s1l[i] = sigma1[i-1];

                }

            else if(  (fabs(sigma1[i-1]) > fabs(sigma1[i])) && (  (sigma1[i-1])*(sigma1[i]) > 0.0 ) )

               {

                   s1l[i] = sigma1[i];

              }
       
           else if ( (sigma1[i-1])*(sigma1[i]) > 0.0  )

               {

                   s1l[i] = 0.0;
            
                }
   }

/*..............................s2l at j+1/2...........................*/


for(i=3;i<=nx-2;i++)
    {
            if(  (fabs(sigma2[i-1]) < fabs(sigma2[i])) && (  (sigma2[i-1])*(sigma2[i]) > 0.0 ) )
               {

                  s2l[i] = sigma2[i-1];

                }

             else if(  (fabs(sigma2[i-1]) > fabs(sigma2[i])) && (  (sigma2[i-1])*(sigma2[i]) > 0.0 ) )

               {

                   s2l[i] = sigma2[i];

              }
       
           else if ( (sigma2[i-1])*(sigma2[i]) > 0.0  )

               {

                   s2l[i] = 0.0;
            
                }
   }



/*..............................s3l at j+1/2...........................*/


for(i=3;i<=nx-2;i++)
    {
            if(  (fabs(sigma3[i-1]) < fabs(sigma3[i])) && (  (sigma3[i-1])*(sigma3[i]) > 0.0 ) )
               {

                  s3l[i] = sigma3[i-1];

                }

           else if(  (fabs(sigma3[i-1]) > fabs(sigma3[i])) && (  (sigma3[i-1])*(sigma3[i]) > 0.0 ) )

               {

                   s3l[i] = sigma3[i];

              }
       
           else if ( (sigma3[i-1])*(sigma3[i]) > 0.0  )

               {

                   s3l[i] = 0.0;
            
                }
   }

/*...........2nd order reconstruction for conserved variables  at left of j+1/2...................*/

for(i=3;i<=nx-2;i++)
    {

         U1l[i]  =  U1[i]  +  0.5*s1l[i];

         U2l[i]  =  U2[i]  +  0.5*s2l[i];

         U3l[i]  =  U3[i]  +  0.5*s3l[i]; 
     }    

      


/*..............................s1r at right of  j+1/2...........................*/

for(i=3;i<=nx-2;i++)
    {
            if(  (fabs(sigma1[i]) < fabs(sigma1[i+1])) && (  (sigma1[i])*(sigma1[i+1]) > 0.0 ) )
               {

                  s1r[i] = sigma1[i];

                }

            else if(  (fabs(sigma1[i]) > fabs(sigma1[i+1])) && (  (sigma1[i])*(sigma1[i+1]) > 0.0 ) )

               {

                   s1r[i] = sigma1[i+1];

              }
       
           else if ( (sigma1[i])*(sigma1[i+1]) > 0.0  )

               {

                   s1r[i] = 0.0;
            
                }
   }

/*..............................s2r at j+1/2...........................*/


for(i=3;i<=nx-2;i++)
    {
            if(  (fabs(sigma2[i]) < fabs(sigma2[i+1])) && (  (sigma2[i])*(sigma2[i+1]) > 0.0 ) )
               {

                  s2r[i] = sigma2[i];

                }

           else if(  (fabs(sigma2[i]) > fabs(sigma2[i+1])) && (  (sigma2[i])*(sigma2[i+1]) > 0.0 ) )

               {

                   s2r[i] = sigma2[i+1];

              }
       
          else if ( (sigma2[i])*(sigma2[i+1]) > 0.0  )

               {

                   s2r[i] = 0.0;
            
                }
   }


/*..............................s3r at j+1/2...........................*/

for(i=3;i<=nx-2;i++)
    {
            if(  (fabs(sigma3[i]) < fabs(sigma3[i+1])) && (  (sigma3[i])*(sigma3[i+1]) > 0.0 ) )
               {

                  s3r[i] = sigma3[i];

                }

           else if(  (fabs(sigma3[i]) > fabs(sigma3[i+1])) && (  (sigma3[i])*(sigma3[i+1]) > 0.0 ) )

               {

                   s3r[i] = sigma3[i+1];

              }
       
           else if ( (sigma3[i])*(sigma3[i+1]) > 0.0  )

               {

                   s3r[i] = 0.0;
            
                }
   }


/*...........2nd order reconstruction for conserved variables  at right of j+1/2...................*/

for(i=3;i<=nx-2;i++)
    {

         U1r[i]  =  U1[i+1]  -  0.5*s1r[i];

         U2r[i]  =  U2[i+1]  -  0.5*s2r[i];					

         U3r[i]  =  U3[i+1]  -  0.5*s3r[i]; 
     }    


     /* .................. COMPUTING FLUXES at left of j+1/2 PART .......................*/
for(i=3;i<=nx-2;i++)
 {
    G1l[i] = U2l[i];
    G2l[i] = ( (U2l[i]*U2l[i])/(U1l[i] + e1) ) + (gamma-1.0)*( U3l[i]  -  ((0.5*U2l[i]*U2l[i])/(U1l[i] + e1)) ) ;
    G3l[i] = (U3l[i]*U2l[i])/(U1l[i] + e1) +  (gamma-1.0)*( U3l[i]  -  ((0.5*U2l[i]*U2l[i])/(U1l[i] + e1)) )*(U2l[i]/(U1l[i] + e1));
 }


     /* .................. COMPUTING FLUXES at right of j+1/2 PART .......................*/
for(i=3;i<=nx-2;i++)
 {
    G1r[i] = U2r[i];
    G2r[i] = ( (U2r[i]*U2r[i])/(U1r[i] + e1) ) + (gamma-1.0)*( U3r[i]  -  ((0.5*U2r[i]*U2r[i])/(U1r[i] + e1)) ) ;
    G3r[i] = (U3r[i]*U2r[i])/(U1r[i] + e1) +  (gamma-1.0)*(  U3r[i]  -  ((0.5*U2r[i]*U2r[i])/(U1r[i] + e1))  )*(U2r[i]/(U1r[i] + e1));
 }

    /*............................Coefficient of diffusion LLF J+1/2................................*/

for(i=3;i<=nx-2;i++)
{
  f1[i] = fabs(u[i]) + a[i];
  f2[i] = fabs(u[i+1]) + a[i+1];
}
for(i=3;i<=nx-2;i++)
 {
  if(f1[i]>f2[i])
          lam_maxph[i] = f1[i];
     else 
         lam_maxph[i] = f2[i];
 }
/* .................................*/
for(i=3;i<=nx-2;i++)
 {
   delU1[i] = (U1r[i] - U1l[i]);
   delU2[i] = (U2r[i] - U2l[i]);
   delU3[i] = (U3r[i] - U3l[i]);

 /*  delU1[i] = (U1[i+1] - U1[i]);
   delU2[i] = (U2[i+1] - U2[i]);
   delU3[i] = (U3[i+1] - U3[i]);*/
 }
  
/* ............................ COMPUTING TOTAL FLUXES AT  J+1/2 INTERFACE...................................*/
for(i=3;i<=nx-2;i++)
  {
    G1ph[i] = 0.5*(G1l[i] + G1r[i]) - 0.5*(lam_maxph[i])*delU1[i];
    G2ph[i] = 0.5*(G2l[i] + G2r[i]) - 0.5*(lam_maxph[i])*delU2[i];
    G3ph[i] = 0.5*(G3l[i] + G3r[i]) - 0.5*(lam_maxph[i])*delU3[i];
  }

/*............................lam_ MIN TOTAL J-1/2................................*/
for(i=3;i<=nx-2;i++)
  {
    f1[i]  =  fabs(u[i-1]) + a[i-1];
    f2[i]  =  fabs(u[i]) + a[i];
  }
for(i=3;i<=nx-2;i++)
{
   if(f1[i]>=f2[i])
        lam_maxmh[i] = f1[i];
       else 
           lam_maxmh[i] = f2[i];
}


/*......................Minmod limiter............................*/

for(i=1;i<=nx;i++)
    {

       sigma1[i]  =  (U1[i] - U1[i-1]);

       sigma2[i]  =  (U2[i] - U2[i-1]);
       
       sigma3[i]  =  (U3[i] - U3[i-1]);


    }


/*..............................s1l at j-1/2...........................*/

for(i=3;i<=nx-2;i++)
    {
            if(  (fabs(sigma1[i-1]) < fabs(sigma1[i])) && (  (sigma1[i-1])*(sigma1[i]) > 0.0 ) )
               {

                  s1l[i] = sigma1[i-1];

                }

            else if(  (fabs(sigma1[i-1]) > fabs(sigma1[i])) && (  (sigma1[i-1])*(sigma1[i]) > 0.0 ) )

               {

                   s1l[i] = sigma1[i];

              }
       
           else if( (sigma1[i-1])*(sigma1[i]) > 0.0 ) 

               {

                   s1l[i] = 0.0;
            
                }
   }

/*..............................s2l at j-1/2...........................*/


for(i=3;i<=nx-2;i++)
    {
            if(  (fabs(sigma2[i-1]) < fabs(sigma2[i])) && (  (sigma2[i-1])*(sigma2[i]) > 0.0 ) )
               {

                  s2l[i] = sigma2[i-1];

                }

           else if(  (fabs(sigma2[i-1]) > fabs(sigma2[i])) && (  (sigma2[i-1])*(sigma2[i]) > 0.0 ) )

               {

                   s2l[i] = sigma2[i];

              }
       
           else if( (sigma2[i-1])*(sigma2[i]) > 0.0  )

               {

                   s2l[i] = 0.0;
            
                }
   }



/*..............................s3l at j-1/2...........................*/


for(i=3;i<=nx-2;i++)
    {
            if(  (fabs(sigma3[i-1]) < fabs(sigma3[i])) && (  (sigma3[i-1])*(sigma3[i]) > 0.0 ) )
               {

                  s3l[i] = sigma3[i-1];

                }

           else if(  (fabs(sigma3[i-1]) > fabs(sigma3[i])) && (  (sigma3[i-1])*(sigma3[i]) > 0.0 ) )

               {

                   s3l[i] = sigma3[i];

              }
       
           else if( (sigma3[i-1])*(sigma3[i]) > 0.0 ) 

               {

                   s3l[i] = 0.0;
            
                }
   }

/*...........2nd order reconstruction for conserved variables  at left of j-1/2...................*/

for(i=3;i<=nx-2;i++)
    {

         U1l[i]  =  U1[i-1]  +  0.5*s1l[i];

         U2l[i]  =  U2[i-1]  +  0.5*s2l[i];

         U3l[i]  =  U3[i-1]  +  0.5*s3l[i]; 
     }    

      



/*..............................s1r at j-1/2...........................*/

for(i=3;i<=nx-2;i++)
    {
            if(  (fabs(sigma1[i]) < fabs(sigma1[i+1])) && (  (sigma1[i])*(sigma1[i+1]) > 0.0 ) )
               {

                  s1r[i] = sigma1[i];

                }

          else if(  (fabs(sigma1[i]) > fabs(sigma1[i+1])) && (  (sigma1[i])*(sigma1[i+1]) > 0.0 ) )

               {

                   s1r[i] = sigma1[i+1];

              }
       
          else if( (sigma1[i])*(sigma1[i+1]) > 0.0 ) 

               {

                   s1r[i] = 0.0;
            
                }
   }

/*..............................s2r at j-1/2...........................*/


for(i=3;i<=nx-2;i++)
    {
            if(  (fabs(sigma2[i]) < fabs(sigma2[i+1])) && (  (sigma2[i])*(sigma2[i+1]) > 0.0 ) )
               {

                  s2r[i] = sigma2[i];

                }

           else if(  (fabs(sigma2[i]) > fabs(sigma2[i+1])) && (  (sigma2[i])*(sigma2[i+1]) > 0.0 ) )

               {

                   s2r[i] = sigma2[i+1];

              }
       
           else if( (sigma2[i])*(sigma2[i+1]) > 0.0 ) 

               {

                   s2r[i] = 0.0;
            
                }
   }


/*..............................s3r at j-1/2...........................*/

for(i=3;i<=nx-2;i++)
    {
            if(  (fabs(sigma3[i]) < fabs(sigma3[i+1])) && (  (sigma3[i])*(sigma3[i+1]) > 0.0 ) )
               {

                  s3r[i] = sigma3[i];

                }

             else if(  (fabs(sigma3[i]) > fabs(sigma3[i+1])) && (  (sigma3[i])*(sigma3[i+1]) > 0.0 ) )

               {

                   s3r[i] = sigma3[i+1];

              }
       
           else if( (sigma3[i])*(sigma3[i+1]) > 0.0 ) 

               {

                   s3r[i] = 0.0;
            
                }
   }


/*...........2nd order reconstruction for conserved variables  at right of j-1/2...................*/

for(i=3;i<=nx-2;i++)
    {

         U1r[i]  =  U1[i]  -  0.5*s1r[i];

         U2r[i]  =  U2[i]  -  0.5*s2r[i];

         U3r[i]  =  U3[i]  -  0.5*s3r[i]; 
     }    





     /* .................. COMPUTING FLUXES at left of j-1/2 PART .......................*/
for(i=3;i<=nx-2;i++)
 {
    G1l[i] = U2l[i];
    G2l[i] = ( (U2l[i]*U2l[i])/(U1l[i] + e1) ) + (gamma-1.0)*( U3l[i]  -  ((0.5*U2l[i]*U2l[i])/(U1l[i] + e1)) ) ;
    G3l[i] = (U3l[i]*U2l[i])/(U1l[i] + e1) +  (gamma-1.0)*( U3l[i]  -  ((0.5*U2l[i]*U2l[i])/(U1l[i] + e1)) )*(U2l[i]/(U1l[i] + e1));
 }


     /* .................. COMPUTING FLUXES at right of j+1/2 PART .......................*/
for(i=3;i<=nx-2;i++)
 {
    G1r[i] = U2r[i];
    G2r[i] = ( (U2r[i]*U2r[i])/(U1r[i] + e1) ) + (gamma-1.0)*( U3r[i]  -  ((0.5*U2r[i]*U2r[i])/(U1r[i] + e1)) ) ;
    G3r[i] = (U3r[i]*U2r[i])/(U1r[i] + e1) +  (gamma-1.0)*( U3r[i]  -  ((0.5*U2r[i]*U2r[i])/(U1r[i] + e1)) )*(U2r[i]/(U1r[i] + e1));
 }

/* ...............................*/
for(i=3;i<=nx-2;i++)
 {
   delU1[i] = (U1r[i] - U1l[i]);
   delU2[i] = (U2r[i] - U2l[i]);
   delU3[i] = (U3r[i] - U3l[i]);

  /* delU1[i] = (U1[i] - U1[i-1]);
   delU2[i] = (U2[i] - U2[i-1]);
   delU3[i] = (U3[i] - U3[i-1]);*/
 }

/* ............................ COMPUTING FLUXES  AT  J-1/2 INTERFACE...................................*/
for(i=3;i<=nx-2;i++)
  {
    G1mh[i] = 0.5*(G1l[i] + G1r[i]) - 0.5*(lam_maxmh[i])*delU1[i];
    G2mh[i] = 0.5*(G2l[i] + G2r[i]) - 0.5*(lam_maxmh[i])*delU2[i];
    G3mh[i] = 0.5*(G3l[i] + G3r[i]) - 0.5*(lam_maxmh[i])*delU3[i];
  }
/* ........................UPDATE FORMULA (Euler forward differences) ...............................*/
for(i=3;i<=nx-2;i++)
{
    U1n[i] = U1[i] - l*(G1ph[i] - G1mh[i]);
    U2n[i] = U2[i] - l*(G2ph[i] - G2mh[i]);
    U3n[i] = U3[i] - l*(G3ph[i] - G3mh[i]);
    dn[i] = U1n[i];
    un[i] = U2n[i]/U1n[i];
    pn[i] = (gamma-1.0)*(U3n[i] - (0.5*U2n[i]*U2n[i]/U1n[i]));
    En[i] = (pn[i]/((gamma-1.0)*dn[i])) + (0.5*un[i]*un[i]);
    an[i] = sqrt(gamma*pn[i]/dn[i]);
    en[i] = pn[i]/((gamma-1.0)*dn[i]);
}
/* ---------------- INITIALIZING -------------------*/
for(i=3;i<=nx-2;i++)
  {
    U1[i] = U1n[i];
    U2[i] = U2n[i];
    U3[i] = U3n[i];
    d[i] = dn[i];
    u[i] = un[i];
    p[i] = pn[i];
    E[i] = En[i];
    a[i] = an[i];
    e[i] = en[i];
 }

t=t+dt;

printf("t=%f | dt=%f\n", t, dt);

if(t>=tn) 
break;
}

for(i=1;i<=nx;i++)
  {
      //  fprintf(fa,"%lf %lf %lf %lf %lf\n",x[i],d[i],u[i],p[i],e[i]);
       fprintf(fa,"%f %f %f %f %f\n",x[i],d[i],u[i],p[i],e[i]);
  }
double end = clock();
// printf("%lf\n", (end - start)/CLOCKS_PER_SEC);


return 0;
}

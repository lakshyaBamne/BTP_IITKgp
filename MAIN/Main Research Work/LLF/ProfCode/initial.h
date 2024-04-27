void initial_wes3(){

      for(int i=1;i<=nx;i++){

            if( x[i] < 0.4){

                  d[i] =  3.857;

                  p[i] =  10.333;

                  a[i] =  sqrt((gamma)*p[i]/d[i]);

                  u[i] =  0.92;

                  E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

                  U1[i] = d[i];

                  U2[i] = d[i]*u[i];

                  U3[i] = d[i]*E[i];

                  H[i] = (U3[i] + p[i])/U1[i] ;

                  M[i]  =  (u[i]/a[i]);

            }

            else{

                  d[i] =  1.0;

                  p[i] =  1.0;

                  a[i] =  sqrt((gamma)*p[i]/d[i]);

                  u[i] =  3.55;

                  E[i] = (p[i]/((gamma-1.0)*d[i])) + 0.5*u[i]*u[i];

                  U1[i] = d[i];

                  U2[i] = d[i]*u[i];

                  U3[i] = d[i]*E[i];

                  p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

                  a[i]=sqrt((gamma)*p[i]/d[i]);

                  H[i] = (U3[i] + p[i])/d[i] ;

                  M[i]  =  (u[i]/a[i]);

            }

      }

}
                                   /*........................................................*/
void initial_wes1()
{

for(i=1;i<=nx;i++)
{
if( x[i] < 0.5)
{

d[i] =  1.0;

p[i] =  1.0;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.0;

//E[i] = (27.0913/3.86);

E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

H[i] = (U3[i] + p[i])/U1[i] ;

M[i]  =  (u[i]/a[i]);

}

else
{

d[i] =  0.125;

p[i] =  0.1;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.0;

E[i] = (p[i]/((gamma-1.0)*d[i])) + 0.5*u[i]*u[i];

//E[i]  =  8.4168;

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

a[i]=sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

M[i]  =  (u[i]/a[i]);

}

}

}


void initial_wes2()
{

for(i=1;i<=nx;i++)
{
if( x[i] < 0.5)
{

d[i] =  0.445;

p[i] =  3.528;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.698;

//E[i] = (27.0913/3.86);

E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

H[i] = (U3[i] + p[i])/U1[i] ;

M[i]  =  (u[i]/a[i]);

}

else
{

d[i] =  0.5;

p[i] =  0.571;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.0;

E[i] = (p[i]/((gamma-1.0)*d[i])) + 0.5*u[i]*u[i];

//E[i]  =  8.4168;

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

a[i]=sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

M[i]  =  (u[i]/a[i]);

}

}

}

/*...................................................................*/
void initial_toro1()
{
/* ........................Initial Conditions.............................*/

for(i=1;i<=nx;i++)
{
if( x[i] < 0.3)
{

d[i] =  1.0;

p[i] =  1.0;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.75;

//E[i] = (27.0913/3.86);

E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

H[i] = (U3[i] + p[i])/U1[i] ;

M[i]  =  (u[i]/a[i]);

}

else
{

d[i] =  0.125;

p[i] =  0.1;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.0;

E[i] = (p[i]/((gamma-1.0)*d[i])) + 0.5*u[i]*u[i];

//E[i]  =  8.4168;

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

a[i]=sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

M[i]  =  (u[i]/a[i]);

}

}

}

/*.........................................................*/

void initial_toro2()
{

for(i=1;i<=nx;i++)
{
if( x[i] < 0.5)
{

d[i] =  1.0;

p[i] =  0.4;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  -2.0;

//E[i] = (27.0913/3.86);

E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

H[i] = (U3[i] + p[i])/U1[i] ;

M[i]  =  (u[i]/a[i]);

}

else
{

d[i] =  1.0;

p[i] =  0.4;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  2.0;

E[i] = (p[i]/((gamma-1.0)*d[i])) + 0.5*u[i]*u[i];

//E[i]  =  8.4168;

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

a[i]=sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

M[i]  =  (u[i]/a[i]);

}

}

}

/*.....................................................*/

/* ......................................STEADY SHOCK..........................*/
void initial_SS()
{


for(i=1;i<=nx;i++)
{

if( x[i] < 0.5)
{

d[i] = 1.0;

u[i] = 1.0;

p[i] = (1.0/(1.4*4));

E[i] = (p[i]/(d[i]*0.4))  +  (0.5*u[i]*u[i]) ;

a[i]  =  sqrt(1.4*p[i]/d[i]);

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

H[i] = ((gamma/(gamma-1.0))*(p[i]/d[i])) + (u[i]*u[i]/2.0);

}

else
{
d[i] =  (56.0/21.0);

u[i] =   sqrt((1.8*9*0.25*21)/(10.8*56));

p[i] = (4.5*0.25)/(1.4);

E[i] = p[i]/(d[i]*0.4)  +  (0.5*u[i]*u[i]) ;

a[i]  =  sqrt(1.4*p[i]/d[i]);

U1[i]  = d[i];

U2[i]  = d[i]*u[i];

U3[i] = d[i]*E[i];

H[i] = ((gamma/(gamma-1.0))*(p[i]/d[i])) + (u[i]*u[i]/2.0);

}

}

}
      /*...................................................*/

void initial_MS()
{


for(i=1;i<=nx;i++)
{
if( x[i] < 0.1)
{

d[i] =  3.86;

u[i] =  -(3.1266/3.86);

//p[i]  =  0.1;

E[i] =   (27.0913/3.86);  

//E[i] = p[i]/((gamma-1.0)*d[i]) + (0.5*u[i]*u[i]);

p[i] =  (gamma - 1.0)*(d[i])*(E[i] - (0.5*u[i]*u[i]));

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

a[i] = sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

e[i] = p[i]/((gamma -1.0)*d[i]);
}

else
{
d[i] =  1.0;

u[i] =  -3.44;

//p[i]  =  0.1;

E[i] = 8.4168;

p[i] = (gamma - 1.0)*(d[i])*(E[i] - (0.5*u[i]*u[i]));

//E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

a[i] = sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

e[i] = p[i]/((gamma -1.0)*d[i]);

}

}

}
                             /*...........................................................*/
void initial_toro3()
{

for(i=1;i<=nx;i++)
{
if( x[i] < 0.5)
{

d[i] =  1.0;

p[i] =  1000.0;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.0;

//E[i] = (27.0913/3.86);

E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

H[i] = (U3[i] + p[i])/U1[i] ;

M[i]  =  (u[i]/a[i]);

}

else
{

d[i] =  1.0;

p[i] =  0.01;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.0;

E[i] = (p[i]/((gamma-1.0)*d[i])) + 0.5*u[i]*u[i];

//E[i]  =  8.4168;

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

a[i]=sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

M[i]  =  (u[i]/a[i]);

}

}

}

/*...........................................................*/
void initial_toro5()
{

for(i=1;i<=nx;i++)
{
if( x[i] < 0.8)
{

d[i] =  1.0;

p[i] =  1000.0;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  -19.59745;

//E[i] = (27.0913/3.86);

E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

H[i] = (U3[i] + p[i])/U1[i] ;

M[i]  =  (u[i]/a[i]);

}

else
{

d[i] =  1.0;

p[i] =  0.01;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  -19.59745;

E[i] = (p[i]/((gamma-1.0)*d[i])) + 0.5*u[i]*u[i];

//E[i]  =  8.4168;

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

a[i]=sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

M[i]  =  (u[i]/a[i]);

}

}

}

/*...........................................................*/
void initial_toro6()
{

for(i=1;i<=nx;i++)
{
if( x[i] < 0.5)
{

d[i] =  1.4;

p[i] =  1.0;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.0;

//E[i] = (27.0913/3.86);

E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

H[i] = (U3[i] + p[i])/U1[i] ;

M[i]  =  (u[i]/a[i]);

}

else
{

d[i] =  1.0;

p[i] =  1.0;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.0;

E[i] = (p[i]/((gamma-1.0)*d[i])) + 0.5*u[i]*u[i];

//E[i]  =  8.4168;

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

a[i]=sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

M[i]  =  (u[i]/a[i]);

}

}

}


/*...........................................................*/
void initial_toro4()
{

for(i=1;i<=nx;i++)
{
if( x[i] < 0.4)
{

d[i] =  5.99924;

p[i] =  460.894;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  19.5975;

//E[i] = (27.0913/3.86);

E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

H[i] = (U3[i] + p[i])/U1[i] ;

M[i]  =  (u[i]/a[i]);

}

else
{

d[i] =  5.99242;

p[i] =  46.0950;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  -6.19633;

E[i] = (p[i]/((gamma-1.0)*d[i])) + 0.5*u[i]*u[i];

//E[i]  =  8.4168;

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

a[i]=sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

M[i]  =  (u[i]/a[i]);

}

}

}

/* ......................................Woodward and Collela Test with Two Discontinuities..........................*/
void initial_blast()

{

for(i=1;i<=nx;i++)
{

if( x[i] < 0.1)
{

d[i] = 1.0;

u[i] = 0.0;

p[i] = 1000.0;

E[i] = (p[i]/(d[i]*0.4))  +  (0.5*u[i]*u[i]) ;


e[i] = (E[i]) - (0.5*u[i]*u[i]);

a[i]  =  sqrt(1.4*p[i]/d[i]);

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

}

if(x[i] > 0.1 && x[i] < 0.9)

{

d[i] =  1.0;

u[i] =  0.0;

p[i] = 0.01;

E[i] = p[i]/(d[i]*0.4)  +  (0.5*u[i]*u[i]) ;

e[i] = (E[i]) - (0.5*u[i]*u[i]);

a[i]  =  sqrt(1.4*p[i]/d[i]);

U1[i]  = d[i];

U2[i]  = d[i]*u[i];

U3[i] = d[i]*E[i];

}


if(x[i] > 0.9)
{

d[i] =  1.0;

u[i] =  0.0;

p[i] = 100.0;

E[i] = p[i]/(d[i]*0.4)  +  (0.5*u[i]*u[i]) ;

e[i] = (E[i]) - (0.5*u[i]*u[i]);

a[i]  =  sqrt(1.4*p[i]/d[i]);

U1[i]  = d[i];

U2[i]  = d[i]*u[i];

U3[i] = d[i]*E[i];

}

}

}

/*...........................................................*/
void initial_toro7()
{

for(i=1;i<=nx;i++)
{
if( x[i] < 0.45)
{

d[i] =  1.4;

p[i] =  1.0;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.1;

//E[i] = (27.0913/3.86);

E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

H[i] = (U3[i] + p[i])/U1[i] ;

M[i]  =  (u[i]/a[i]);

}

else
{

d[i] =  1.0;

p[i] =  1.0;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.1;

E[i] = (p[i]/((gamma-1.0)*d[i])) + 0.5*u[i]*u[i];

//E[i]  =  8.4168;

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

a[i]=sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

M[i]  =  (u[i]/a[i]);

}

}

}


/*...........................................................*/






/*...........................................................*/
void initial_toro8()
{

for(i=1;i<=nx;i++)
{
if( x[i] < 0.5)
{

d[i] =  1.0;

p[i] =  0.1;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  2.0;

//E[i] = (27.0913/3.86);

E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

H[i] = (U3[i] + p[i])/U1[i] ;

M[i]  =  (u[i]/a[i]);

}

else
{

d[i] =  1.0;

p[i] =  0.1;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  -2.0;

E[i] = (p[i]/((gamma-1.0)*d[i])) + 0.5*u[i]*u[i];

//E[i]  =  8.4168;

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

a[i]=sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

M[i]  =  (u[i]/a[i]);

}

}

}


/*...........................................................*/
void initial_over_heating()
{

for(i=1;i<=nx;i++)
{
if( x[i] < 0.5)
{

d[i] =  1.0;

p[i] =  1.0;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  -0.8*a[i];

//E[i] = (27.0913/3.86);

E[i] = p[i]/((gamma-1.0)*d[i]) + 0.5*u[i]*u[i];

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

H[i] = (U3[i] + p[i])/U1[i] ;

M[i]  =  (u[i]/a[i]);

}

else
{

d[i] =  1.0;

p[i] =  1.0;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.8*a[i];

E[i] = (p[i]/((gamma-1.0)*d[i])) + 0.5*u[i]*u[i];

//E[i]  =  8.4168;

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

a[i]=sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

M[i]  =  (u[i]/a[i]);

}

}

}


/*....................................................................*/
void initial_shock_entropy()

{

for(i=1;i<=nx;i++)
{
if( x[i] < -0.8)
{

d[i] =  3.857143;

p[i] =  10.3333;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  2.629369;

//E[i] = (27.0913/3.86);

E[i] = p[i]/((gamma-1.0)*d[i]) + (0.5*u[i]*u[i]);

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

H[i] = (U3[i] + p[i])/d[i] ;

M[i]  =  (u[i]/a[i]);

}

else
{

d[i] =  1.0 + 0.2*sin(5*3.14*x[i]);

p[i] =  1.0;

a[i] =  sqrt((gamma)*p[i]/d[i]);

u[i] =  0.0;

E[i] = (p[i]/((gamma-1.0)*d[i])) + (0.5*u[i]*u[i]);

//E[i]  =  8.4168;

U1[i] = d[i];

U2[i] = d[i]*u[i];

U3[i] = d[i]*E[i];

//p[i] =  (gamma-1.0)*(U3[i]  -  ((0.5*U2[i]*U2[i])/(U1[i])))  ;

//a[i]=sqrt((gamma)*p[i]/d[i]);

H[i] = (U3[i] + p[i])/d[i] ;

M[i]  =  (u[i]/a[i]);


}

}

}

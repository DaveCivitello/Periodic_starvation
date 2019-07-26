/* file IndividualModel_YVE.c */
#include <R.h>
#include <math.h>
static double parms[25];

#define iM parms[0]
#define k parms[1]
#define M parms[2]
#define EM parms[3]
#define Fh parms[4]
#define muD parms[5]
#define DR parms[6]
#define fe parms[7]
#define yRP parms[8]
#define ph parms[9]
#define yPE parms[10]
#define iPM parms[11]
#define eh parms[12]
#define mP parms[13]
#define alpha parms[14]
#define yEF parms[15]
#define LM parms[16]
#define kd parms[17]
#define z parms[18]
#define kk parms[19]
#define hb parms[20]
#define theta parms[21]
#define mR parms[22]
#define yVE parms[23]
#define startage parms[24]


/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=25;
odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
double *yout, int *ip)
{
if (ip[0] <1) error("nout should be at least 1");

double Chi = M/(1 + EM);
double fH = y[0]/(y[0]+Fh);
double fP = y[2]/(y[2] + eh);
double L = y[1];
double VOL = L*L*L;
double SA = L*L;
double fbf = 0;

if( t[0] < startage){fH = fe;}else{fbf = 1;};
double Dens = y[5]/(Chi*VOL);
double kstar = fmin(k + y[5]*alpha, 1);
double aM = iM*yEF;
double mV = aM*k/(LM*Chi);
double mD = muD*mV;
double rp = Dens*Dens/(ph*ph + Dens*Dens);
double Jec = y[2]/(1 + yVE*kstar*EM*y[2])*(aM*SA + yVE*EM*(mV+mR*EM*y[7])*Chi*VOL);

double Delta = fmax(0,(mV+mR*EM*y[7])*Chi*VOL - kstar*Jec);
double Delta2 = (mV+mR*EM*y[7])*Chi*VOL - Jec;

ydot[0] = -iM*SA*fH*fbf;
ydot[1] = fmax(yVE/(3*Chi)*(kstar*aM*y[2] - (mV+mR*EM*y[7])*Chi*L)/(1 + yVE*kstar*EM*y[2]),0);

ydot[2] = aM/(Chi*EM*L)*(fH - y[2]) - iPM*Dens*fP/EM;
ydot[5] = yPE*iPM*fP*(1 - rp)*y[5] - mP*y[5];
ydot[6] = yRP*yPE*iPM*fP*rp*y[5];
if(y[3] < DR){
  ydot[3] = (1 - kstar)*Jec - fmax(mD*y[3],0) - Delta;
  ydot[4] = 0;}else{
  ydot[3] = fmin(0, (1 - kstar)*Jec - mD*DR - Delta);
  ydot[4] = fmax((1 - kstar)*Jec - mD*DR  - Delta, 0);}
if(Delta2 > 0){
  ydot[2] = 1/(EM*Chi*VOL)*(aM*SA*fH - ((mV+mR*EM*y[7])*Chi*VOL  + iPM*fP*y[5]));
  ydot[3] = -mD*fmax(0,fmin(y[3],DR));
  ydot[4] = 0;}
ydot[7] = theta/(Chi*VOL)*ydot[6] + kd*(1-y[2]) - kd*y[7] - 3*y[7]*ydot[1]/L;
ydot[8] = kk*fmax(y[7] - z, 0) + hb;
if(y[2] <= 0){
  ydot[0] = 0;
  ydot[1] = 0;
  ydot[2] = 0;
  ydot[3] = 0;
  ydot[4] = 0;
  ydot[5] = 0;
  ydot[6] = 0;
  ydot[7] = 0;
  ydot[8] = hb;}

  yout[0] = exp(-y[8]);

}

/* END file IndividualModel_YVE.c */ 

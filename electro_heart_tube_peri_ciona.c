
//----------------------------------------------------------------------------------------------------------------//
// Include statements, don't worry about changing these.  //
//----------------------------------------------------------------------------------------------------------------//

#include <stdio.h>
#include <math.h>  
#include <stdlib.h>

//#include "myfft.h"  
//#include <fftw.h>  
#include <fftw3.h>
  

//----------------------------------------------------------------------------------------------------------------//
// Declaring variables and allocating memory, don't need to change any of this! //
//----------------------------------------------------------------------------------------------------------------//



double *u, *v, *w1,  *w2, *f1, *f2, *a, *b, *c, *d;
double *theta_1, *theta_2;
double bptime, bdptime, ap;  
double *unew, *vnew, *psi, *prest, *btu, *btv, *ut, *vt, *dt1p, *dt2p, *dt1p2, *dt2p2, damp;  
double *d13, *d23, *dt13, *dt23, *dt1p3, *dt2p3, *f1t3, *f2t3, *f1t4, *f2t4, *d14, *d24, *dt14, *dt24, *dt1p4, *dt2p4; 
double *dtr, *d1, *d2, *d12, *d22, *f1r, *f2r, *f1t, *f2t;  
double *dfinal, *deltax, *deltay, *d_bp;  
double  *bt1, *bt2, *f1b, *f2b, *b1, *b2, *force1r, *force2r;  
double *force1l, *force2l, *tforce1r, *tforce2r, *force1m, *force2m, *marker1, *marker2;  
double *markeru, *markerv, *du, *dv, *du2, *dv2, *du3, *dv3, *du4, *dv4, *bu, *bv, *dt1, *dt2, *dt12, *dt22, *ut2, *vt2;  
double *vorticity, *mesh, *dmesh, *force1r, *force2r;
double inertialf, viscousf, tlift, ttranslation;  
int counter, plan, itop;  
double tr1, da, tm, dm, pos, ave_vel;
double atau1, atau2;
double thetan1, thetan2;
double FHpulse, tscale, *kelect;
int Nact, Np, WI3, WI4;
double *restds, *velect, *welect, *f1act, *f2act, *f12act, *f22act;
double FHpulse, tscale, pulse_time, I_mag, lengthFH, dxFH, dtFH, dpFH, aFH, bFH, gammaFH, DFH;
double *welect_old, *velect_old, *II;
  
//Functions that allocate memory  
int *alloc_1d_int(int n1);  
void free_1d_int(int *i);  
double *alloc_1d_double(int n1);  
void free_1d_double(double *d);  
int **alloc_2d_int(int n1, int n2);  
void free_2d_int(int **ii);  
double **alloc_2d_double(int n1, int n2);  
void free_2d_double(double **dd);  
  
fftw_complex *mymatrix, *uft, *vft, *uftout, *vftout, *w1ft, *w2ft, *w1ftout, *w2ftout, *w1ftin, *w2ftin, *pft; 
fftw_plan plu1, plv1, plw1, plw2, pl, pl1, pl2, pl3, pl4;  
  
int n; 
int markersize1, markersize2, markersize3; 
  
double mo;  
double mn;  
int isn;  
int i;  
int j;  
int n1, n2, n3;  
int no;  
int int1;  
int int2;  
int i_1, ip1, j_1, jp1;  
int M;  
int N;  
int ii;  
int jj;  
double dt_pdx2;  
double dtnu;  
double nu;  
int x1;  
int x2;  
int ns;  
double hi;  
double ptime;  
double width;  
double dptime;  
double rspeed;  
  
double thetar, ES;  
double thetal;  
double omega;  
double ds;  
double length;  
double EI, stiff;  
double transdistance;  
double timeturn;  
double timetranslate;  
double thetaup;  
double thetadown;  
double nu;  
double flu_err;
  
double scale;  
double p;	//density of fluid  
double mu; //viscosity  
double pi;  
double Do1;  
double Do12;  
double Do2;  
double Do21;  
double Dpm1;  
double Dpm12;  
double Dpm2;  
double Dpm21;  
double dt;  
double dx;  
int k;  
int h;  
double ttime;  
double time;  
double dt_dx;  
  
double sign;  
double ks, kb;  
double d_dx;  
double dx_dt;  
double dt_pdx;  
double dtmu;  
double dt_p;  
double d_2dx;  
double d_dx2;  
double num;  
double lift;  
double translation;  
double err;  
double EIds4; 
double Re; 
double *d1s, *d2s, *f1s, *f2s, *f1s2, *f2s2, *f1t2, *f2t2, *f1r2, *f2r2;  
double time_final;
double changept;
double porosity;
int corner;

//----------------------------------------------------------------------------------------------------------------//
// Declaring number of points in fluid grid and along the boundary, shouldn't need to be changed//
//----------------------------------------------------------------------------------------------------------------//

  
int M = 512;  //number of fluid grid pts
int N = 512;  //number of fluid grid pts
int WI; //number of boundary pts on one side 
int WI2; //number of boundary pts on the other side
int Nact; //number of boundary points in active region.
int Nside; //number of boundary points on each side
    
int FHN_fp; //beginning of excitable domain
int FHN_ip; //end of excitable domain
int p_start; //beginning of pacemaker
int p_end;   //end of pacemaker

//Main IBM functions  
void mvwing(double array1[], double array2[], double a_wing[], double angle, int WI, int WI2, double length, double sgn, double ds);  
void mv2wing(double array2[], int WI, double rspeed, double mo);  

void framp(double f1[], double f2[], double u[], double v[], double k, double ttime, double ramptime, double ramp, double rampmax, int M, int N);
void pulse(double dt1[], double dt12[], double freq, double amp, double diameter, double ttime);
void fcalc(double f_1[], double array1[], double array_t1[], double vel1[], double vel_t1[],double stiff, int L);   
void fcalc1(double f_1[], double array1[], double array_t1[], double EIds4, int WI);  
void fcalc2(double f_1[], double f_2[], double array1[], double array2[], double stiff, double ds, int L); 
void fcalc3(double f_1[], double f_2[], double f_12[], double f_22[], double array1[], double array2[], double array12[], double array22[], double stiff[], double dsv[], int L);  
void fcalc4(double f_1[], double array1[], double array_t1[], double stiff, int L); 
void fcalc_tree(double f1t[], double f2t[], double f1t2[], double f2t2[], double d1[], double d2[], double d12[], double d22[], double ks, double ds, int WI);
void fspread(double *f1, double *f2, double *array1, double *array2,  double *a_f1, double *a_f2, double ds, double dx, int L);  
void intwing(double *u, double *v, double array1[], double array2[], double *uarray, double *varray, double dx, double dt, int L);  
void intpwing(double *u, double *v, double arrayf1[], double arrayf2[],  double array1[], double array2[], double *uarray, double *varray, double dx, double dt, int L, double porosity, double angle); 
void flu(double *u, double *v, double *f1, double *f2, double *a, double *b, double *c, double dx, double dt);  
void makew(double *u, double *v, double *f1, double *f2, double dx, double dt, double *w1, double *w2, int M, int N);  
void reporting(double zz);  
double showforce(double *force, int M, int N, double ttime);  
  
//Miscallaneous  
void vortex(double *u, double *v, double *vorticity, int M, int N, double dx);  
double trans, transn, chord, tau, d_tt, d_tr, time1, time2, time3, timer, d_taut, d_taur, ramp, ramptime, rampmax;
double tau1, tau2, tau3, taur, middlex, middley, middley2, theta, theta2, tau_final, omegan, omegao;
double diameter, R1, R2, Ltube, thetac, center1x, center1y, center2x, center2y, freq, percont, amp, angfreq, Wo;

//zeroing functions  
void zero2(double *array, int L, int O);  
void zero1(double array[], int L);  
  
//record to file functions  
double recordf(double dptime, double *v, double f1r[], double f2r[], double f1s[], double f2s[], double f1t[], double f2t[], double force1r[], double force2r[], double f1[], double f2[], double ttime);  
double record(double *u, double *v, double *vorticity, double dptime, double ttime); 
double recordu(double *u, double dptime, double ttime); 
double recordv(double *v, double dptime, double ttime);  
double recordo(double dptime, double d1[], double d2[], double d12[], double d22[], double marker1[], double marker2[], double ttime); 
double recordp2(double dptime, double d13[], double d23[], double d14[], double d24[], double ttime); 
double recordavevel(double *u, double *v, double dptime, double ttime);  
double recordp(double dptime, double *prest, double ttime);
void record1real(double *array, int L);  
double record2(double *array, int M, int N);  
void record1(double *d2, double *d1, double L, double ttime);  
double recordpart(double dptime, double marker1[], double marker2[], double ttime);
  
int row(int ir, int jr);  
double zz;  
int vortest;
  
//checks
int check_fluid, check_forcespread, check_velocity, norm_run, fftcheck;
double a11, a12, a21, a22, b_1, b_2, bk_11, bk_12, bk_21, bk_22;
double w_1, w_2, k_11, k_12, k_21, k_22;

FILE *out;  
FILE *out1;  
FILE *out2;  
FILE *out3;  
FILE *out4;
FILE *out5;
FILE *out6;
FILE *out11;
  

//----------------------------------------------------------------------------------------------------------------//
// beginning of the main file, this is where all changable parameters will be! //
//----------------------------------------------------------------------------------------------------------------//


int main()  
{  
//Opens files to write to
  out5=fopen("particles_valveless_test","w"); //x,y coordinates of particles
  out6 = fopen("ave_vel_test","w"); //average velocity along tube vs. time
  out=fopen("markers_valveless_test","w");  //x,y coordinates of boundaries
  out1=fopen("vorticity_valveless_test","w");  //vorticity
  out2=fopen("forces_valveless_test", "w");  //forces on walls (ignore for now)
  out3=fopen("velu_valveless_test","w");  //velocity in u direction
  out4=fopen("velv_valveless_test","w"); //velocity in v direction
  out11=fopen("pmarkers_valveless_test","w"); //pericardial markers_valveless_test
  
  //Declare variables  
  plan = 1;
  check_fluid=0;
  check_forcespread=0;
  check_velocity=0;
  norm_run=1;
  fftcheck=1;

  //Initializations 
  pi=4*atan(1);  
  ttime = 0.000;       //real time  
  lift=0;				//zero initial lift and translation  
  translation =0;   
  vortest = 1; 
   
  //Kinematic parameters -- don't remember why these are what they are, but they work out in the old paper.
 

//----------------------------------------------------------------------------------------------------------------//
// Time can be changed if we want longer simulations. Wo is the womersley number (a dimensionless parameter which describes the fluid dynamics of the system) 
// and will need to be changed for certain simulations. Time_final doesn't need to be changed. 
// To change Wo: just fill in a different number! Be sure to document all changes made in a separate
// README file. For example:  Wo = 10, changed from 0.9 
//----------------------------------------------------------------------------------------------------------------//

  freq = 0.554; //heart beat frequency of typical Ciona
  time = 8*(1/freq);  //seconds
  //Wo = 0.1; //Womersley number. Default is 0.082. Old:0.9
  time_final = time;

//----------------------------------------------------------------------------------------------------------------//
// dt will need to be changed (made smaller) if you notice that the simulation isn't running for it's entirety //
//----------------------------------------------------------------------------------------------------------------//

  //Timing variables
  dptime = time/1800;   //print forces every dptime. Ok to keep this number big  
  bdptime = time/160;  //print velocities every bdptime. Don't make this too big
  time_final = time;
  dt = 5*4*4*20*time/(5*5*2*2*4*1.2*2*2.5*3.5*200000);     //Make dt smaller if things blow up adding a *4 on the denominator 
  changept = 1;

//----------------------------------------------------------------------------------------------------------------//
// Diameter can be changed to alter the diameter of the tube (not the pericardium). Freq controls frequency 
// of pumping and can also be changed to test fluid velocity v. frequency.  BE SURE TO DOCUMENT EVERYTHING // 
//----------------------------------------------------------------------------------------------------------------//

  //Spatial variables for heart tube
  length = 0.015;     //length of the box
  width = 0.015; 
  diameter = 0.0006957; //diameter of tube old: 0.000035
  R2 = 0.0024;  //outer radius
  R1 = R2-diameter; //inner radius
  Ltube = 0.00791; //length straight part

  //freq = 1/0.9; //sec^-1 old: 2.3
  angfreq = 2*pi*freq;
  percont = 0.749; //percent of tube contracted
  amp = 0.5*(diameter*percont); //amplitude contraction

  dx = length/(N);        //spatial step
  ds = length/(2*(N));   //boundary step  
 
  Nact = floor(Ltube/ds); //number of points straight part of tube
  Nside = floor(pi*R2/ds); //number of points each side

  WI = Nact+Nact+Nside+Nside;
  WI2 = WI;

//----------------------------------------------------------------------------------------------------------------//
// ap is distance the pericardium is expanded above and below the pumping tube. ap will need to be changed if testing
// how the diameter of the pericarium is affecting the fluid dynamics. The total diameter of the pericardium is: 
// 2*ap + diameter. The other parameters here shouldn't need to be changed. 
//----------------------------------------------------------------------------------------------------------------//

  //pericardial stuff
  ap = 0.7*diameter;
  Np = floor(ap/ds);
  WI3 = Nact+Np+Np;
  WI4 = Nact+Np+Np;
  
  //--------------

  //Mechanical variables
  //Re = 10; //Reynolds number for simulation
  p= 1025.0;  //kg/m^3 
  //nu = angfreq*(0.5*diameter/Wo)*(0.5*diameter/Wo); //this gets modified for changing Wo
  mu=0.001002;   //dynamic viscosity in N s/m^2//p*nu; 
  nu = mu/p; //kinematic viscosity
  
  EI = 0.00000000000233795; // (1000*10*0.0027*5/(ds*100))*ds*ds*ds*ds; //this is picked as a good guess
  kb = 410848865.94/100;  //boundary stiffness  
  ks = 410848865.94;     //stretching stiffness  
  ES = kb;

//----------------------------------------------------------------------------------------------------------------//
// These are the fitzHugh parameters and shouldn't need to be changed. If you want to change the dynamics of the pumping
// then they can be changed. BE SURE TO DOCUMENT EVERYTHING IN A README FILE! //
//----------------------------------------------------------------------------------------------------------------//

  //----------------------------
  //FitzHugh Nagumo parameters
  FHpulse = 128; // this is used to scale time so it pulses at right frequency
  tscale = 0; //
  pulse_time = 0; //keeps track of when to pulse, set to zero initially
  I_mag = 0.75; //Magnitude of current injected
  lengthFH = 1000; //length of domain from nagumo1.m
  dxFH = lengthFH/Nact; //spatial step size for FH solver
  dtFH = dt*FHpulse/(1/freq); //time step size for FH solver
  aFH = 0.01; //threshold potential
  bFH = 0.01; //epsilon in FH model, sets time scale
  gammaFH = 0.5; //gamma used in FH model, inverse of strength of blocking
  DFH = 100; //diffusion coefficient from nagumo1.m
  dpFH = FHpulse/200;
  //--------------------------
  //kelect = ks/100;
  damp = 0.0; //50000/(16*8*8);
  d_dx=1/dx;  
  dt_pdx2=dt/(p*dx*dx);  
  d_dx = 1.0/dx;  
  dt_p=dt/p; 
  EIds4 = EI/(ds*ds*ds*ds);  

  //Printing variables
  ptime= 0.000;     //keeps track of when to print  
  bptime = 0.000;

  //-----------------------------------------------------------
  //------------------------------------------------------------
  //set values for the check
//  if(check_fluid==1){

    //Mechanical variables
//    mu = 3.0;  
//    p= 2.0; 
//    nu = mu/p;

    //New spatial variables
//    length = 0.4;  
//    width = 0.4;  
//    dx = length/N;        //absolute length of grid  
//    ds = length/(2*N);   //boundary step  

    //new timing parameter values
//    time = 0.2; 	//total time
//    dptime = time/10;   //print every dptime  
//    bdptime = time/60;  
//    time_final = time;
//    dt = dx*dx/nu;      //time step  

//    d_dx=1/dx;  
//    dt_pdx2=dt/(p*dx*dx);  
//    d_dx = 1.0/dx;  
//    dt_p=dt/p; 
//  }
  //-------------------------------------------------------------  
  //-------------------------------------------------------------



  //allocate memory  
   
  n2 = M;  
  n1 = N;  
  prest = (double *) malloc(n1*n2*sizeof(double));  
  unew = (double *) malloc(n1*n2*sizeof(double));  
  vnew = (double *) malloc(n1*n2*sizeof(double));  
  psi = (double *) malloc((n1/2)*(n2/2)*sizeof(double)); 
  mesh = alloc_1d_double(256);  
  u = (double *) malloc(n1*n2*sizeof(double));  
  v = (double *) malloc(n1*n2*sizeof(double));  
  d = (double *) malloc(n1*n2*sizeof(double));  
  a = (double *) malloc(n1*n2*sizeof(double));  
  b = (double *) malloc(n1*n2*sizeof(double));  
  c = (double *) malloc(n1*n2*sizeof(double));  
  f1 = (double *) malloc(n1*n2*sizeof(double));  
  f2 = (double *) malloc(n1*n2*sizeof(double));  
  vorticity = (double *) malloc(n1*n2*sizeof(double)); 
  force1r = alloc_1d_double(WI);  
  force2r = alloc_1d_double(WI);  
  tforce1r = alloc_1d_double(WI);  
  tforce2r = alloc_1d_double(WI);  
  dtr = alloc_1d_double(WI);  
  dt1 = alloc_1d_double(WI);  
  dt2 = alloc_1d_double(WI);  
  dt12 = alloc_1d_double(WI2);  
  dt22 = alloc_1d_double(WI2); 
  d1 = alloc_1d_double(WI);  
  d2 = alloc_1d_double(WI); 
  d12 = alloc_1d_double(WI2);  
  d22 = alloc_1d_double(WI2); 
  dt13 = alloc_1d_double(WI3);  
  dt23 = alloc_1d_double(WI3);  
  dt14 = alloc_1d_double(WI4);  
  dt24 = alloc_1d_double(WI4); 
  d13 = alloc_1d_double(WI3);  
  d23 = alloc_1d_double(WI3); 
  d14 = alloc_1d_double(WI4);  
  d24 = alloc_1d_double(WI4); 
  
  btu = alloc_1d_double(4*(N-30)); 
  btv = alloc_1d_double(4*(N-30)); 
  ut = alloc_1d_double(WI); 
  vt = alloc_1d_double(WI); 
  ut2 = alloc_1d_double(WI2); 
  vt2 = alloc_1d_double(WI2); 
  dt1p = alloc_1d_double(WI);
  dt2p = alloc_1d_double(WI);  
  dt1p2 = alloc_1d_double(WI2);
  dt2p2 = alloc_1d_double(WI2);  
  dt1p3 = alloc_1d_double(WI3);
  dt2p3 = alloc_1d_double(WI3);  
  dt1p4 = alloc_1d_double(WI4);
  dt2p4 = alloc_1d_double(WI4);  
  f1r = alloc_1d_double(WI);  
  f2r = alloc_1d_double(WI);  
  f1t = alloc_1d_double(WI);  
  f2t = alloc_1d_double(WI);  
  deltax = alloc_1d_double(WI);  
  deltay = alloc_1d_double(WI);  
  d_bp = alloc_1d_double(WI);  
 

  dfinal = alloc_1d_double(WI);  
  
  bt1 = alloc_1d_double(4*(N-30));  
  bt2 = alloc_1d_double(4*(N-30));  
  b1 = alloc_1d_double(4*(N-30));  
  b2 = alloc_1d_double(4*(N-30));  
  f1b = alloc_1d_double(4*(N-30));  
  f2b = alloc_1d_double(4*(N-30));  
  du = alloc_1d_double(WI);  
  dv = alloc_1d_double(WI);  
  du2 = alloc_1d_double(WI2);  
  dv2 = alloc_1d_double(WI2);  
  du3 = alloc_1d_double(WI3);  
  dv3 = alloc_1d_double(WI3);   
  du4 = alloc_1d_double(WI4);  
  dv4 = alloc_1d_double(WI4); 
  bv = alloc_1d_double(4*(N-30));  
  bu = alloc_1d_double(4*(N-30));  
  markeru = alloc_1d_double(30*12);  
  markerv = alloc_1d_double(30*12);  
  dmesh = alloc_1d_double(M/2);  
  d1s = alloc_1d_double(WI);  
  d2s = alloc_1d_double(WI);  
  f1s = alloc_1d_double(WI);  
  f2s = alloc_1d_double(WI);  
  f1s2 = alloc_1d_double(WI2);  
  f2s2 = alloc_1d_double(WI2);  
  f1t2 = alloc_1d_double(WI2);  
  f2t2 = alloc_1d_double(WI2); 
  f1r2 = alloc_1d_double(WI2);  
  f2r2 = alloc_1d_double(WI2); 
  f1t3 = alloc_1d_double(WI3);  
  f2t3 = alloc_1d_double(WI3);
  f1t4 = alloc_1d_double(WI4);  
  f2t4 = alloc_1d_double(WI4);  

  restds = alloc_1d_double(Nact);
  velect = alloc_1d_double(Nact);
  welect = alloc_1d_double(Nact);
  velect_old = alloc_1d_double(Nact);
  welect_old = alloc_1d_double(Nact);
  kelect = alloc_1d_double(Nact);
  II = alloc_1d_double(Nact);
  f1act = alloc_1d_double(Nact);
  f2act = alloc_1d_double(Nact);
  f12act = alloc_1d_double(Nact);
  f22act = alloc_1d_double(Nact);

  w1 = (double *) malloc(n1*n2*sizeof(double));  
  w2 = (double *) malloc(n1*n2*sizeof(double));  
  
  marker1 = alloc_1d_double(30*12);  
  marker2 = alloc_1d_double(30*12);  
 
 
  //in = fftw_malloc(sizeof(fftw_complex) * N);
  //out = fftw_malloc(sizeof(fftw_complex) * N);
  uft = fftw_malloc(n2*n1*sizeof(fftw_complex));  
  vft = fftw_malloc(n2*n1*sizeof(fftw_complex));
  uftout = fftw_malloc(n2*n1*sizeof(fftw_complex));  
  vftout = fftw_malloc(n2*n1*sizeof(fftw_complex));  
  w1ft = fftw_malloc(n2*n1*sizeof(fftw_complex));  
  w2ft = fftw_malloc(n2*n1*sizeof(fftw_complex));  
  w1ftout = fftw_malloc(n2*n1*sizeof(fftw_complex));  
  w2ftout = fftw_malloc(n2*n1*sizeof(fftw_complex));  

  mymatrix = fftw_malloc(3*5*sizeof(fftw_complex));  
  pft = fftw_malloc(n2*n1*sizeof(fftw_complex));

  //makes plans
  pl1 = fftw_plan_dft_2d(M, N, w1ft, w1ftout, FFTW_FORWARD, FFTW_MEASURE);
  pl2 = fftw_plan_dft_2d(M, N, w2ft, w2ftout, FFTW_FORWARD, FFTW_MEASURE);
  pl3 = fftw_plan_dft_2d(M, N, uft, uftout, FFTW_BACKWARD, FFTW_MEASURE);
  pl4 = fftw_plan_dft_2d(M, N, vft, vftout, FFTW_BACKWARD, FFTW_MEASURE);
  
  //--------------------------------------------------------------  
  //Set initial conditions  
  zero2(u,M,N);            //for fluid velocity  
  zero2(v,M,N);  
  zero2(vorticity,M,N);  
  zero1(force1r, WI);  
  zero1(force2r, WI);  

  zero2(f1,M,N);            //zero forces applied to fluid  
  zero2(f2,M,N);  
    
     
  //----------------------------------------------  
  //set inital values for boundary points and boundary target points  
  corner = floor(0.5*M);

  //bottom straight piece
  for(i=0;i<Nact;i++){  
    d1[i] = corner*ds;  
    d2[i] = (i+corner)*ds;  
    d12[i] = (corner)*ds+diameter;  
    d22[i] = (i+corner)*ds; 
  }  
  
  //pericardium
  for(i=0;i<Np;i++){
    d13[i] = corner*ds-i*ds;  
    d23[i] = (corner)*ds;  
    d14[i] = (corner)*ds+diameter+i*ds;  
    d24[i] = (corner)*ds; 
  }
  
  for(i=Np;i<Nact+Np;i++){
    d13[i] = corner*ds-Np*ds;  
    d23[i] = (corner)*ds+(i-Np)*ds;  
    d14[i] = (corner)*ds+diameter+Np*ds;  
    d24[i] = (corner)*ds+(i-Np)*ds; 
  }
  
  for(i=Np+Nact;i<Nact+Np+Np;i++){
    d13[i] = corner*ds-Np*ds+(i-Np-Nact)*ds;  
    d23[i] = (corner)*ds+Nact*ds;  
    d14[i] = (corner)*ds+diameter+Np*ds-(i-Np-Nact)*ds;  
    d24[i] = (corner)*ds+Nact*ds; 
  }
  
  center1x = corner*ds+R2;
  center1y = (corner+Nact)*ds;
  
  //right end
  for(i=Nact;i<Nact+Nside;i++){
    thetac = (i-Nact)*pi/Nside-pi/2;
    d1[i] = center1x+R2*sin(thetac);  
    d2[i] = center1y+R2*cos(thetac); 
    d12[i] = center1x+R1*sin(thetac);  
    d22[i] = center1y+R1*cos(thetac);  
  }  

  //straight top piece
  for(i=Nact+Nside;i<Nact+Nact+Nside;i++){  
    d1[i] = corner*ds+2*R2;  
    d2[i] = (corner+Nact-(i-Nact-Nside))*ds;  
    d12[i] = (corner)*ds+2*R1+diameter;  
    d22[i] = (corner+Nact-(i-Nact-Nside))*ds;
  }  

  center2x = corner*ds+R2;
  center2y = corner*ds;

  //left end
  for(i=Nact+Nact+Nside;i<Nact+Nact+Nside+Nside;i++){
    thetac = pi/2+(i-Nact-Nact-Nside)*pi/Nside;
    d1[i] = center2x+R2*sin(thetac);  
    d2[i] = center2y+R2*cos(thetac); 
    d12[i] = center2x+R1*sin(thetac);  
    d22[i] = center2y+R1*cos(thetac);  
  } 

  //set target points to initial boundary position
  for(i=0;i<WI;i++){  
    dt1[i] = d1[i];  
    dt2[i] = d2[i];
    dt1p[i] = dt1[i];
    dt2p[i] = dt2[i];  
  }  
  
  for(i=0;i<WI2;i++){  
    dt12[i] = d12[i];  
    dt22[i] = d22[i];
    dt1p2[i] = dt12[i];
    dt2p2[i] = dt22[i];  
  }  
  
  for(i=0;i<WI3;i++){  
    dt13[i] = d13[i];  
    dt23[i] = d23[i];
    dt1p3[i] = dt13[i];
    dt2p3[i] = dt23[i];  
  }
  
  for(i=0;i<WI4;i++){  
    dt14[i] = d14[i];  
    dt24[i] = d24[i];
    dt1p4[i] = dt14[i];
    dt2p4[i] = dt24[i];  
  }
  
  for(i=0;i<4*(N-30);i++){
    btu[i]=0;
    btv[i]=0;
  }

for(i=0;i<Nact;i++){
kelect[i] = ks/100;
velect[i] = 0;
welect[i] = 0;
velect_old[i]=0;
welect_old[i]=0;
 restds[i]=(1-percont)*diameter;
II[i]=0;
}

  //-----------------------------------------------------------  
  //Set initial values for particles 
  for(i=0;i<12;i++)  
    for(j=0;j<30;j++){  
      marker1[i*30+j] = (corner+(diameter/(6*length))*Nact*i)*ds;  
      marker2[i*30+j] = (corner+j*Nact*(Ltube/(length*15)))*ds;  
    }  
  
  //------------------------------------------------------------------  
  //Set values for matrix used in fluid solver.  
   for(i=0;i<M;i++)  
    for(j=0;j<N;j++){  
      if(((i==0) && (j==0))||((i==0) && (j==N/2))||((i==M/2) && (j==0))||((i==M/2) && (j==N/2))){  
	a[row(i,j)]=0;  
	b[row(i,j)]=0;  
	c[row(i,j)]=0;  



      }  
      else{  
	a[row(i,j)]= sin(2*pi*i/M)*sin(2*pi*i/M)/((sin(2*pi*i/M)*sin(2*pi*i/M)+sin(2*pi*j/N)*sin(2*pi*j/N)));  
	b[row(i,j)]= sin(2*pi*i/M)*sin(2*pi*j/N)/((sin(2*pi*i/M)*sin(2*pi*i/M)+sin(2*pi*j/N)*sin(2*pi*j/N)));  
	c[row(i,j)]= sin(2*pi*j/N)*sin(2*pi*j/N)/((sin(2*pi*i/M)*sin(2*pi*i/M)+sin(2*pi*j/N)*sin(2*pi*j/N)));  
      }  
    }  
    
  for(i=0;i<M;i++)  
    for(j=0;j<N;j++){  
      d[row(i,j)]=1/(1+(4*nu*dt/(dx*dx))*(sin(pi*i/M)*sin(pi*i/M)+sin(pi*j/N)*sin(pi*j/N)));  
    } 
  mo = -1;  
  counter = 0;
 
  //-------------------------------------------------------
  //-------------------------------------------------------

  //------------------------------------------------------
  //------------------------------------------------------
  

  if(norm_run==1){
    while(ttime<time)  
      {  
  
	//print to output file every dptime.
	//recordp(bdptime, prest, ttime);  
        recordpart(dptime, marker1, marker2, ttime);  
	recordp2(dptime, d13, d23, d14, d24, ttime); 
	recordo(dptime, d1, d2, d12, d22, marker1, marker2, ttime);  
	recordu(u, bdptime, ttime);  
	recordv(v, bdptime, ttime); 
	record(u, v, vorticity, bdptime, ttime);   
        recordavevel(u, v, dptime, ttime);  
	recordf(dptime, v, f1r, f2r, f1s, f2s, f1t, f2t, force1r, force2r, f1, f2, ttime); 
       

	//zero matrices  
	zero1(f1s, WI);  
	zero1(f2s, WI);  
	zero1(f1s2, WI2);  
	zero1(f1t2, WI2);  
	zero1(f2t2, WI2);
	zero1(f1r2, WI2);  
	zero1(f2s2, WI2);
	zero1(f2r2, WI2);
	zero1(f1r, WI);  
	zero1(f2r, WI);  
	zero1(f1t, WI);  
	zero1(f2t, WI);  
	zero1(f1b, (4*(N-30)));  
	zero1(f2b, (4*(N-30)));
	zero1(f1act, Nact);
	zero1(f2act, Nact);  
	zero1(f12act, Nact);
	zero1(f22act, Nact);
	zero2(f1, M, N);  
	zero2(f2, M, N);
	
/*	while ('\n' != getchar ());
	printf ("Press enter to continue...");
	getchar ()*/;

	zero1(f1t3, WI3);  
	zero1(f2t3, WI3);
	zero1(f1t4, WI4);  
	zero1(f2t4, WI4);
	//These markers are here to move around in the fluid. Never use them.

	//scaled time for FH solver
	tscale = ttime*(FHpulse/(1/freq));

        //this function drives the section of the tube 
        //pulse(dt1, dt12, freq, amp, diameter, ttime);
	//FH solver
	FHN_ip = floor(0.075*Nact); //beginning of excitable domain
	FHN_fp = floor(0.85*Nact); //end of excitable domain
	p_start = floor(0.175*Nact); //beginning of pacemaker
	p_end = floor(0.225*Nact);   //end of pacemaker

	for(i=0;i<FHN_ip-1;i++){
	  velect[i] = 0;
	  welect[i] = 0;
	}
	velect[FHN_ip-1] = 0;
	welect[FHN_ip-1] = welect_old[FHN_ip-1]+dtFH*(bFH*(velect_old[FHN_ip-1]-gammaFH*welect_old[FHN_ip-1]));

    if(tscale>(pulse_time)){
        for(i=p_start;i<p_end;i++){
            II[i] = I_mag;
        }}
    if(tscale>(pulse_time+dpFH)){
        pulse_time = pulse_time+FHpulse;
        }
	if(tscale<(pulse_time)){
        for(i=0;i<Nact;i++){
            II[i] = 0;
        }
    }


	for(i=FHN_ip;i<FHN_fp;i++){
        velect[i] = velect_old[i]+dtFH*(((DFH/(dxFH*dxFH))*(velect_old[i+1]-2*velect_old[i]+velect_old[i-1]))-velect_old[i]*(velect_old[i]-1)*(velect_old[i]-aFH)-welect_old[i]+II[i]);
        welect[i] = welect_old[i]+dtFH*(bFH*(velect_old[i]-gammaFH*welect_old[i]));
    }
    
	velect[FHN_fp] = 0;
	welect[FHN_fp] = welect_old[FHN_fp]+dtFH*(bFH*(velect_old[FHN_fp]-gammaFH*welect_old[FHN_fp]));
	
	for(i=FHN_fp;i<Nact;i++){
	velect[i]=0;
	welect[i]=0;
	}

	for(i=0;i<Nact;i++){
	velect_old[i]=velect[i];
	welect_old[i]=welect[i];
	}
	
	for(i=0;i<Nact;i++){
		//restds[i]=percont*diameter*(1-velect[i]);
	  kelect[i]=(ks/100)*(velect[i])*(velect[i])*(velect[i])*(velect[i]);
		if(velect[i]<0){
			kelect[i]=0;
		}
	}
	//---------------------------------------------------  
	
	//For now, the only force on the wings is going to be the target force.	  
	//Calculate forces.  

        //Uncomment these if you want bending stiffness
        fcalc1(f1r, d1, dt1, EIds4, WI);           //elastic force on wing  
	fcalc1(f2r, d2, dt2, EIds4, WI);  
	fcalc1(f1r2, d12, dt12, EIds4, WI);           //elastic force on wing  
	fcalc1(f2r2, d22, dt22, EIds4, WI); 
	fcalc(f1s, d1, dt1, du, ut, kb, WI);		  //target force on wing  
	fcalc(f2s, d2, dt2, dv, vt, kb, WI);  
	fcalc(f1s2, d12, dt12, du2, ut2, kb, WI2);		  //target force on wing  
	fcalc(f2s2, d22, dt22, dv2, vt2, kb, WI2); 
	fcalc2(f1t, f2t, d1, d2, ks, ds, WI);     //stretching force on wing  
	fcalc2(f1t2, f2t2, d12, d22, 4*ks, ds, WI);  
	fcalc3(f1act, f2act, f12act, f22act, d1, d2, d12, d22, kelect, restds, Nact);  
	fcalc4(f1t3, d13, dt13, kb, WI3);		  //target force on wing  
	fcalc4(f2t3, d23, dt23, kb, WI3); 
	fcalc4(f1t4, d14, dt14, kb, WI4);		  //target force on wing  
	fcalc4(f2t4, d24, dt24, kb, WI4);
	
	//zero old velocities
	zero1(du, WI);  
	zero1(dv, WI);    
    zero1(du2, WI2);  
    zero1(dv2, WI2); 
	zero1(du3, WI3);  
    zero1(dv3, WI3); 
	zero1(du4, WI4);  
    zero1(dv4, WI4); 
	zero1(bu, (4*(N-30)));  
	zero1(bv, (4*(N-30))); 
  
	//This commented stuff down here is used when you want to remove the springs and let the wing be fexible in certain parts.
	//----------------------------------------------------  
	//remove target springs from bottom strip, not removing all pieces. The 2 prevents breakage.  
	for(i=2; i<Nact-1;i++){
	  f1s[i]=0;
	  //f2s[i]=0;
	  	f1s2[i] = 0;
	 	//f2s2[i] = 0;
	 }

	//framp(f1, f2, u, v, ks, ttime, ramptime, ramp, rampmax, M, N);
	//Spread all of the different forces to the fluid grid.
	fspread(f1, f2, d1, d2, f1r, f2r, ds, dx, WI);  
	fspread(f1, f2, d12, d22, f1r2, f2r2, ds, dx, WI);  
	fspread(f1, f2, d1, d2, f1t, f2t, ds, dx, WI);  
	fspread(f1, f2, d12, d22, f1t2, f2t2, ds, dx, WI2);
	fspread(f1, f2, d1, d2, f1act, f2act, ds, dx, Nact);  
	fspread(f1, f2, d12, d22, f12act, f22act, ds, dx, Nact);
	fspread(f1, f2, d1, d2, f1s, f2s, ds, dx, WI);  
	fspread(f1, f2, d12, d22, f1s2, f2s2, ds, dx, WI2);  
 	fspread(f1, f2, d13, d23, f1t3, f2t3, ds, dx, WI3); 
	fspread(f1, f2, d14, d24, f1t4, f2t4, ds, dx, WI4); 
		
	//Solve fluid  
	flu(u, v, f1, f2, a, b, c, dx, dt);  

  			  
	//Interpolate and move boundary.  
	//This commented part gets used when you don't want the wing to be porous.
	//intwing(u, v, d1, d2, du, dv, dx, dt, WI);         //wing  
        //intwing(u, v, d12, d22, du2, dv2, dx, dt, WI2); 
	thetan1=pi/2+theta;
	thetan2=pi/2+theta2;
	intpwing(u, v, f1s, f2s, d1, d2, du, dv, dx, dt, WI, porosity, thetan1); 
	intpwing(u, v, f1s2, f2s2, d12, d22, du2, dv2, dx, dt, WI2, porosity, thetan2);
	intwing(u, v, marker1, marker2, markeru, markerv, dx, dt,(12*30));          //boundary  
 	intwing(u, v, d13, d23, du3, dv3, dx, dt, WI3);  
	intwing(u, v, d14, d24, du4, dv4, dx, dt, WI4);  	
	//----------------------------------------------------  
 
	//Update time, do other stuff to get ready for the next time step.  
	ttime = ttime+dt; 
 
  	for(i=0;i<WI;i++){
	  dt1p[i] = dt1[i];
	  dt2p[i] = dt2[i];
	}
  	for(i=0;i<WI2;i++){
		dt1p2[i] = dt12[i];
		dt2p2[i] = dt22[i];
	}

	//zero again 
	markersize1 = 30*12;  
	zero1(markeru, markersize1);  
	zero1(markerv, markersize1);  
      } } 

    
  fclose(out);  
  fclose(out1);  
  fclose(out2);  
  fclose(out3);  
}  
  
void reporting(double zz){  
  fprintf(out3,"\n");  
  fprintf(out3, "there was an error in %g", zz);  
  fclose(out);  
  exit(-1);  
}  
  
int row(int ir, int jr){  
  int ij;  
  ij = ir*N+jr;  
  return ij;  
}  
  
  
void mvwing(double array1[], double array2[], double a_wing[], double angle, int WI, int WI2, double length, double sgn, double ds){  
//spin wing  

tau = ttime*trans/chord;

//acceleration at the beginning of the stroke
if((tau>0) && (tau<tau1)){
transn = trans*.5*(1+cos(pi+pi*((tau)/d_taut)));
}

//constant translation upstroke
if((tau>tau1) && (tau<tau2)){
transn = trans;
}

//deceleration at end of upstroke
if((tau >tau2) && (tau<tau3)){
transn = trans*.5*(1+cos(pi*(tau - tau2)/d_taut));
}

//constant upstroke
if(tau<(atau1)){
omegan = 0;
}

if((tau > atau1) && (tau < (atau2))){
omegan = omega*0.5*(1-cos(2*pi*(tau-atau1)/d_taur));
}

middley = middley + transn*dt;
middley2 = middley2 - transn*dt;
theta = theta + omegan*dt;
theta2 = theta2-omegan*dt;

if(tau<atau1){
  for (i=0; i<WI; i++){                  
    dt1[i] = middlex + i*ds*sin(theta);  
    dt2[i] = middley + i*ds*cos(theta);  
  }  

  for (i=0; i<WI2; i++){                  
    dt12[i] = middlex + i*ds*sin(theta2);  
    dt22[i] = middley2 + i*ds*cos(theta2);  
  } 
}


if((tau>atau1) && (tau<atau2)){
   if(changept==1){ 
      middlex = middlex + (WI-1)*ds*sin(theta);
      middley = middley + (WI-1)*ds*cos(theta);
      middley2 = middley2 - (WI-1)*ds*cos(theta); 
      changept=0;
 }
  for (i=0; i<WI; i++){                  
    dt1[WI-1-i] = middlex - i*ds*sin(theta);  
    dt2[WI-1-i] = middley - i*ds*cos(theta);  
  }  

  for (i=0; i<WI2; i++){                  
    dt12[WI2-1-i] = middlex - i*ds*sin(theta2);  
    dt22[WI2-1-i] = middley2 - i*ds*cos(theta2);  
  } 
}

  
  for(i=0;i<WI;i++){
  ut[i] = (dt1[i]-dt1p[i])/dt;
  vt[i] = (dt2[i]-dt2p[i])/dt;
  }

  for(i=0;i<WI2;i++){
  ut2[i] = (dt12[i]-dt1p2[i])/dt;
  vt2[i] = (dt22[i]-dt2p2[i])/dt;
  }
 } 
  
  
  
void mv2wing(double array2[], int WI, double rspeed, double mo){  
  //translate  
  for (i=0;i<WI;i++){  
    array2[i] = array2[i]+mo*rspeed*dt;  
    vt[i] = mo*rspeed;
    ut[i] = 0;
  }  
}  

void pulse(double dt1[], double dt12[], double freq, double amp, double diameter, double ttime){
  for(i=150; i<180;i++){
    dt1[i] = 300*ds+amp*0.5*(1+sin(2*pi*freq*ttime-pi/2)); 
    dt12[i] = (300)*ds+diameter-amp*0.5*(1+sin(2*pi*freq*ttime-pi/2)); 
  }
}

void framp(double f1[], double f2[], double u[], double v[], double ks, double ttime, double ramptime, double ramp, double rampmax, int M, int N){

	if (ttime<ramptime){
		ramp=rampmax*ttime/ramptime;
	}
		if (ttime>ramptime){
		ramp=rampmax;
	}
	for(i=15;i<M-15;i++){
		for(j=60;j<65;j++){
			pos = i-15;
			f2[row(i,j)]=(ks/10)*(ramp*(1-(1-pos/300)*(1-pos/300)) - v[row(i,j)]);
			f1[row(i,j)]=(ks/10)*(0 -u[row(i,j)]);
			//f1[row(i,j)]=k*(ramp*(1-(1-pos/300)*(1-pos/300)) - u[row(i,j)]);
			//f2[row(i,j)]=k*(0 -v[row(i,j)]);
		}
	}
}

void fcalc1(double f_1[], double array1[], double array_t1[], double EIds4, int WI){  
  //caculate elastic forces  
  
  for(i=0;i<WI;i++){  
    f_1[i]=0;}  
  
  for(i=1; i<WI-1;i++){  
    f_1[i+1] = f_1[i+1]-EIds4*(array1[i+1]+array1[i-1]-2*array1[i]);  
    f_1[i-1] = f_1[i-1]-EIds4*(array1[i+1]+array1[i-1]-2*array1[i]);  
    f_1[i] = f_1[i]+2*EIds4*(array1[i+1]+array1[i-1]-2*array1[i]);	  
  }  
}  
  
void fcalc2(double f_1[], double f_2[], double array1[], double array2[], double stiff, double ds, int L){  
  //calculate forces from stretching  

  for(i=0; i<L-1;i++){  
    deltax[i] = array1[i+1]-array1[i];  
    deltay[i] = array2[i+1]-array2[i];  
    d_bp[i] = sqrt(deltax[i]*deltax[i]+deltay[i]*deltay[i]);  
  }  
  f_1[0] = -(deltax[0]/d_bp[0])*stiff*(ds-d_bp[0]);  
  f_2[0] = -(deltay[0]/d_bp[0])*stiff*(ds-d_bp[0]);  
  
  for(i=1; i<L-1; i++){  
    f_1[i] = (deltax[i-1]/d_bp[i-1])*stiff*(ds-d_bp[i-1]) - (deltax[i]/d_bp[i])*stiff*(ds-d_bp[i]);  
    f_2[i] = (deltay[i-1]/d_bp[i-1])*stiff*(ds-d_bp[i-1]) - (deltay[i]/d_bp[i])*stiff*(ds-d_bp[i]);  
  }  
  f_1[L-1] = (deltax[L-2]/d_bp[L-2])*stiff*(ds-d_bp[L-2]);  
  f_2[L-1] = (deltay[L-2]/d_bp[L-2])*stiff*(ds-d_bp[L-2]);  
}  

void fcalc3(double f_1[], double f_2[], double f_12[], double f_22[], double array1[], double array2[], double array12[], double array22[], double stiff[], double dsv[], int L)
{
  for(i=0; i<L;i++){  
    deltax[i] = array1[i]-array12[i];  
    deltay[i] = array2[i]-array22[i];  
    d_bp[i] = sqrt(deltax[i]*deltax[i]+deltay[i]*deltay[i]);  
  } 
  for(i=0; i<L; i++){  
    f_1[i] = (deltax[i]/d_bp[i])*stiff[i]*(dsv[i]-d_bp[i]);  
    f_2[i] = (deltay[i]/d_bp[i])*stiff[i]*(dsv[i]-d_bp[i]);  
    f_12[i] = -(deltax[i]/d_bp[i])*stiff[i]*(dsv[i]-d_bp[i]);  
    f_22[i] = -(deltay[i]/d_bp[i])*stiff[i]*(dsv[i]-d_bp[i]); 
  }
}  

void fcalc_tree(double f1t[], double f2t[], double f1t2[], double f2t2[], double d1[], double d2[], double d12[], double d22[], double ks, double ds, int WI){
	deltax[1] = d12[2] - d1[WI/2 + 2];
	deltay[1] = (d22[2] - d2[WI/2 + 2]);
	d_bp[1] = sqrt(deltax[1]*deltax[1]+deltay[1]*deltay[1]); 
	
	f1t2[2] =f1t2[2]+(deltax[1]/d_bp[1])*ks*(sqrt(2)*2*ds-d_bp[1]);
	f2t2[2] = f2t2[2]+(deltay[1]/d_bp[1])*ks*(sqrt(2)*2*ds-d_bp[1]);
	f1t[WI/2 + 2] = f1t[WI/2 + 2]-(deltax[1]/d_bp[1])*ks*(sqrt(2)*2*ds-d_bp[1]);
	f2t[WI/2 + 2] = f2t[WI/2 + 2]-(deltay[1]/d_bp[1])*ks*(sqrt(2)*2*ds-d_bp[1]);
	
	deltax[2] = d12[2] - d1[WI/2 - 2];
	deltay[2] = (d22[2] - d2[WI/2 - 2]);
	d_bp[2] = sqrt(deltax[2]*deltax[2]+deltay[2]*deltay[2]); 
	
	f1t2[2] =f1t2[2]+(deltax[2]/d_bp[2])*ks*(sqrt(2)*2*ds-d_bp[2]);
	f2t2[2] = f2t2[2]+(deltay[2]/d_bp[2])*ks*(sqrt(2)*2*ds-d_bp[2]);
	f1t[WI/2 - 2] = f1t[WI/2 - 2]-(deltax[2]/d_bp[2])*ks*(sqrt(2)*2*ds-d_bp[2]);
	f2t[WI/2 - 2] = f2t[WI/2 - 2]-(deltay[2]/d_bp[2])*ks*(sqrt(2)*2*ds-d_bp[2]);

}

  
void fcalc(double f_1[], double array1[], double array_t1[], double vel1[], double vel_t1[], double stiff, int L){  
  //calculate target forces  
  for(i=0; i<L;i++){  
    f_1[i] = stiff*(array_t1[i] - array1[i]) + damp*(vel_t1[i]-vel1[i]);  
  }  
}  

void fcalc4(double f_1[], double array1[], double array_t1[], double stiff, int L){
   for(i=0; i<L;i++){  
    f_1[i] = stiff*(array_t1[i] - array1[i]);  
  }  
}  

void fspread(double *f1, double *f2, double *array1, double *array2,  double *a_f1, double *a_f2, double ds, double dx, int L){  
//spread forces from wing   
  for(i=0;i<L;i++){  
    x1 = (int)ceil(array1[i]/dx);  
    x2 = (int)ceil(array2[i]/dx);  
    
    for(k=0; k<4; k++){  
      for(h=0; h<4; h++){  
	ii=x1+1-k;  
	jj=x2+1-h;  
	  
	f1[row(ii,jj)] = f1[row(ii,jj)]+a_f1[i]*(1/(4*dx))*(1+cos((pi/(2*dx))*(dx*(x1+1-k)-array1[i])))*(1/(4*dx))*(1+cos((pi/(2*dx))*(dx*(x2+1-h)-array2[i])))*ds;  
	f2[row(ii,jj)] = f2[row(ii,jj)]+a_f2[i]*(1/(4*dx))*(1+cos((pi/(2*dx))*(dx*(x1+1-k)-array1[i])))*(1/(4*dx))*(1+cos((pi/(2*dx))*(dx*(x2+1-h)-array2[i])))*ds;  
      }}}  
}  
  

  
void record1(double *d1, double *d2,  double L, double ttime){  
  //record 1-d array  
  fprintf(out,"\n");  
  for (i=0; i<L; i++){  
    fprintf(out, "%g\t%g\t%g", d1[i], d2[i], ttime);  
    fprintf(out,"\n");  
  }  
}  
  
void intwing(double *u, double *v, double array1[], double array2[], double *uarray, double *varray, double dx, double dt, int L){  
//move wing at local fluid velocity  
  for(i=0; i<L; i++){  
    x1 = (int)ceil(array1[i]/dx);  
    x2 = (int)ceil(array2[i]/dx);  
      
    for(k=0; k<4; k++){  
      for(h=0; h<4; h++){  
	uarray[i] = uarray[i]+u[row(x1+1-k, x2+1-h)]*(1/(4*dx))*(1+cos((pi/(2*dx))*(dx*(x1+1-k)-array1[i])))*(1/(4*dx))*(1+cos((pi/(2*dx))*(dx*(x2+1-h)-array2[i])))*dx*dx;  
	varray[i] = varray[i]+v[row(x1+1-k, x2+1-h)]*(1/(4*dx))*(1+cos((pi/(2*dx))*(dx*(x1+1-k)-array1[i])))*(1/(4*dx))*(1+cos((pi/(2*dx))*(dx*(x2+1-h)-array2[i])))*dx*dx;  
      }}  
      
    array1[i] = array1[i] + uarray[i]*dt;  
    array2[i] = array2[i] + varray[i]*dt;  
  }  
  
}
  
void intpwing(double *u, double *v, double arrayf1[], double arrayf2[],  double array1[], double array2[], double *uarray, double *varray, double dx, double dt, int L, double porosity, double angle){ 
  for(i=0; i<L; i++){  
    x1 = (int)ceil(array1[i]/dx);  
    x2 = (int)ceil(array2[i]/dx);  
      
    for(k=0; k<4; k++){  
      for(h=0; h<4; h++){  
	uarray[i] = uarray[i]+u[row(x1+1-k, x2+1-h)]*(1/(4*dx))*(1+cos((pi/(2*dx))*(dx*(x1+1-k)-array1[i])))*(1/(4*dx))*(1+cos((pi/(2*dx))*(dx*(x2+1-h)-array2[i])))*dx*dx;  
	varray[i] = varray[i]+v[row(x1+1-k, x2+1-h)]*(1/(4*dx))*(1+cos((pi/(2*dx))*(dx*(x1+1-k)-array1[i])))*(1/(4*dx))*(1+cos((pi/(2*dx))*(dx*(x2+1-h)-array2[i])))*dx*dx;  
      }}  
      
    array1[i] = array1[i] + (uarray[i]+porosity*(arrayf1[i]*cos(angle)+arrayf2[i]*sin(angle))*cos(angle))*dt;  
    array2[i] = array2[i] + (varray[i]+porosity*(arrayf1[i]*cos(angle)+arrayf2[i]*sin(angle))*sin(angle))*dt;  
  }  
  
}

double recordpart(double dptime, double marker1[], double marker2[], double ttime){  
  //record motion of wing and markers  
  if(ttime>bptime)  
    {  
      fprintf(out5,"\n");  
      for (i=0; i<360; i++){  
	fprintf(out5, "%g\t%g\t%g", marker1[i], marker2[i], ttime);  
	fprintf(out5,"\n");  
      }  
  
      fflush(out5);
      //ptime = ptime + dptime;  
    }  
    
  return bptime;  
} 

double recordo(double dptime, double d1[], double d2[], double d12[], double d22[], double marker1[], double marker2[], double ttime){  
  //record motion of wing and markers  
  if(ttime>bptime)  
    {  
      fprintf(out,"\n");  
      for (i=0; i<WI; i++){  
	fprintf(out, "%g\t%g\t%g\t%g\t%g", d1[i], d2[i], d12[i], d22[i], ttime);  
	fprintf(out,"\n");  
      }  
  
      fflush(out);
      //ptime = ptime + dptime;  
    }  
    
  return bptime;  
}  

double recordp2(double dptime, double d13[], double d23[], double d14[], double d24[], double ttime){
  //record motion of wing and markers  
  if(ttime>bptime)  
    {  
      fprintf(out11,"\n");  
      for (i=0; i<WI3; i++){  
	fprintf(out11, "%g\t%g\t%g\t%g\t%g", d13[i], d23[i], d14[i], d24[i], ttime);  
	fprintf(out11,"\n");  
      }  
  
      fflush(out11);
      //ptime = ptime + dptime;  
    }  
    
  return bptime;  
}   
  
double showforce(double *force, int M, int N, double ttime){  
  if(ttime>ptime)  
    {  
      fprintf(out,"\n");  
      for (i=0; i<M; i++){  
	for (j=0; j<N; j++){  
	  fprintf(out, "%g\t%g\t%g\t%g", force[row(i,j)], i*ds, j*ds, ttime);  
	  fprintf(out,"\n");  
	}}  
      fflush(out);  
    }  
  return ptime;  
}  
  
double recordavevel(double *u, double *v, double dptime, double ttime){  
  if(ttime>ptime){
    itop = 150+floor((diameter/length)*600);
    for(i=150; i<itop; i++){
    ave_vel=v[row(i,450)];
    fprintf(out6,"\n");  
    fprintf(out6, "%g\t%g", ave_vel, ttime);
    }
    //ave_vel=ave_vel/(itop-150);  
    fflush(out6);
    ave_vel=0;
}}


double recordf(double dptime, double *v, double f1r[], double f2r[], double f1s[], double f2s[], double f1t[], double f2t[], double force1r[], double force2r[], double f1[], double f2[], double ttime){  
  //record lift and translation  
  //for(i=0;i<WI;i++){  
  //  tforce1r[i] = f1[i]; //f1s[i] + f1t[i]+ tforce1r[i];  
  //  tforce2r[i] = f2[i]; //f2s[i] + f2t[i]+ tforce2r[i];  
  //}  
  
  //counter = counter+1;  
  
  if(ttime>ptime){  

    for(i=9*WI/10;i<WI;i++){  
      force1r[i] = f1r[i]+f1s[i]+f1t[i];  
      force2r[i] = f2s[i]+f2r[i]+f2t[i];
      tforce1r[i] = f1s[i]+f1s2[i]; //tforce1r[i]/counter;  
      tforce2r[i] = f2s[i]+f2s2[i]; //tforce2r[i]/counter;
    }  
 
    counter = 0; 
    lift = 0;
    translation = 0; 
    tlift = 0;
    ttranslation = 0;

    for(i=0;i<WI;i++){  
      lift = lift - force1r[i];  
      translation = translation - force2r[i];  
      tlift = tlift - tforce1r[i];
      ttranslation = ttranslation -tforce2r[i];
    }  

    fprintf(out2,"\n");  
    fprintf(out2, "%g\t%g\t%g\t%g\t%g", lift, translation, tlift, ttranslation, ttime);  
    ptime= ptime+dptime;  

    //zero the forces
    for(i=0;i<WI;i++){  
      force2r[i] = 0;  
      force1r[i] = 0;  
      tforce2r[i] = 0;
      tforce1r[i] = 0;
    }  
    lift = 0;  
    translation = 0;  
    tlift = 0;
    ttranslation = 0;

    fflush(out2);  
    ptime= ptime+dptime;  
  }  
  return ptime;  
} 

double recordp(double dptime, double *prest, double ttime)
{
  if((ttime>bptime))  
    {  

      for(i=0;i<M/2;i++){
	for(j=0;j<(N/2);j++){
	  fprintf(out4, "%g\t", prest[row(2*i, 2*j)]);
	}
	fprintf(out4, "\n");
      }
      fflush(out4);  
    }    
  return bptime;
} 
  

//double record(double *u, double *v, double *vorticity, double dptime, double ttime){  
//record vorticity  
//   if((ttime>bptime))  
//    {  
//for(j=1;j<(N/2);j++){
//psi[row(0,j)] = psi[row(0,j-1)] + 2*dx*u[row(0, (2*j))];
//}

//for(i=1;i<M/2;i++){
//for(j=0;j<(N/2);j++){
//psi[row(i,j)] = psi[row(i-1,j)] - 2*dx*v[row((2*i), (2*j))];
//}}

//     for(i=0;i<M/2;i++){  
//	  for(j=0;j<N/2;j++){  
//	    fprintf(out1, "%g\t", psi[row(i,j)]);  
//         }
//	    fprintf(out1, "\n");  
//	    }  
//      fflush(out1);  
//      bptime = bdptime + bptime;
//    }  
//  return bptime;  
//}   

double record(double *u, double *v, double *vorticity, double dptime, double ttime){  
  //record vorticity  
//Note that this does not record vorticity for the whole domain. I only plotted some of it because back then there wasn't much hard drive space.
   if((ttime>bptime))  
    {  

      vortex(u, v, vorticity, M, N, dx);  
      fprintf(out1,"\n");  
        
      for(i=0;i<M;i++){  
	  for(j=0;j<N;j++){  
	    fprintf(out1, "%g\t", vorticity[row(i,j)]);  
         }
	    fprintf(out1, "\n");  
	    }  
     fflush(out1);  
      bptime = bdptime + bptime;
    }  
  return bptime;  
}    
  
double recordu(double *u, double dptime, double ttime){  
  //recordu-velocity

   if((ttime>bptime))  
    {  
  
      fprintf(out3,"\n");  
        
      for(i=0;i<M;i++){  
	  for(j=0;j<N;j++){  
	    fprintf(out3, "%g\t", u[row(i,j)]);  
         }
	    fprintf(out3, "\n");  
	    }  
     fflush(out3);  
    }  
  return bptime;  
}  

  
double recordv(double *v, double dptime, double ttime){  
  //record v-velocity

   if((ttime>bptime))  
    {  
  
      fprintf(out4,"\n");  
        
      for(i=0;i<M;i++){  
	  for(j=0;j<N;j++){  
	    fprintf(out4, "%g\t", v[row(i,j)]);  
         }
	    fprintf(out4, "\n");  
	    }  
     fflush(out4);  
    }  
  return bptime;  
}  
void vortex(double *u, double *v, double *vorticity, int M, int N, double dx){  
  for(i=1;i<(M-1);i++){  
    for(j=1;j<(N-1);j++){  
      vorticity[row(i,j)] = (v[row(i+1,j)]-v[row(i-1,j)])/(2*dx) - (u[row(i,j+1)]-u[row(i,j-1)])/(2*dx);  
    }}  
}  

  
double record2(double *array, int M, int N){  
  //record 2-d array  
  fprintf(out,"\n");  
  for (i=0; i<M; i++){   
    for (j=0;j<N;j++){  
      fprintf(out, "%g\t",array[row(i,j)]);  
    }  
    fprintf(out,"\n");  
  }  
  fflush(out);  
}  
  
  
void record1real(double *array, int L){  
  //record 1-d real array  
  fprintf(out,"\n");  
  for (i=0; i<L; i++){  
    fprintf(out, "%g\n",array[i]);  
  }  
  fflush(out);  
}  
  
void zero2(double *array, int L, int O){  
  //zero2d array  
  for(i=0; i<L; i++)  
    for(j=0;j<O;j++){  
      array[row(i,j)]=0.0;  
    }  
}  
  
void zero1(double array[], int L){  
  //zero 1-d array  
  for(i=0; i<L; i++){  
    array[i]=0.0;  
  }  
}  
  
void flu(double *u, double *v, double *f1, double *f2, double *a, double *b, double *c, double dx, double dt)  
//fluid solver  
{  
  
  pi=4*atan(1);  
  mn=M*N;  
  scale = 1.0/(mn);  
  
  //-----------------------------------------------------------------  
   
  zero2(w1, M, N);  
  zero2(w2, M, N);  
    
  makew(u, v, f1, f2, dx, dt, w1, w2, M, N);  
  for (i=0; i<M; i++)  
    for (j = 0; j<N; j++){  
      w1ft[i*N+j][0] = w1[row(i,j)];  
      w2ft[i*N+j][0] = w2[row(i,j)];  
      w1ft[i*N+j][1] = 0;  
      w2ft[i*N+j][1] = 0;  
    }  
  //if(fftcheck==1){  
  //pl = fftw2d_create_plan(M,N,FFTW_FORWARD, FFTW_MEASURE|FFTW_IN_PLACE);}
  
  //take forward transforms  
  fftw_execute(pl1);  
  fftw_execute(pl2);  
  //fftwnd_destroy_plan(pl);  
  
  //----------------------------------------------------------  
  
  
  for(i=0;i<M;i++)  
    for(j=0;j<N;j++){  
      uft[i*N+j][0] = (w1ftout[i*N+j][0]-(a[row(i,j)]*w1ftout[i*N+j][0]+b[row(i,j)]*w2ftout[i*N+j][0]))*d[row(i,j)];  
      uft[i*N+j][1] = (w1ftout[i*N+j][1]-(a[row(i,j)]*w1ftout[i*N+j][1]+b[row(i,j)]*w2ftout[i*N+j][1]))*d[row(i,j)];  
      vft[i*N+j][0] = (w2ftout[i*N+j][0]-(b[row(i,j)]*w1ftout[i*N+j][0]+c[row(i,j)]*w2ftout[i*N+j][0]))*d[row(i,j)];  
      vft[i*N+j][1] = (w2ftout[i*N+j][1]-(b[row(i,j)]*w1ftout[i*N+j][1]+c[row(i,j)]*w2ftout[i*N+j][1]))*d[row(i,j)];  
    }  
  
  
  //fftwnd_destroy_plan(pl); 
  //if(fftcheck==1){ 
  //  pl2 = fftw2d_create_plan(M,N,FFTW_BACKWARD, FFTW_MEASURE|FFTW_IN_PLACE);
  //  fftcheck=0;}
  
  //take backward transforms  
  fftw_execute(pl3);  
  fftw_execute(pl4); 
  //fftwnd_destroy_plan(pl);  
  
  for(i=0;i<M;i++)  
    for(j=0;j<N;j++){  
      u[row(i,j)] = scale*uftout[i*N+j][0];  
      v[row(i,j)] = scale*vftout[i*N+j][0];  
	}
   

}  
  
void makew(double *u, double *v, double *f1, double *f2, double dx, double dt, double *w1, double *w2, int M, int N)  
{  
  //Functions calculate the value of w1.   
  //pi=4*atan(1);  

  
  //----------------------------------------------------------------  
  
  //Fill in boudary values, using periodic domain.  
  //Fill in top    

  for(i=0;i<M;i++){
    for(j=0;j<N;j++){
jp1 = j+1;
ip1=i+1;
i_1 = i-1;
j_1 = j-1;
      if(i==(M-1)){
	ip1 = 0;}
      if(i==0){
	i_1 = M-1;}
      if(j==(N-1)){
	jp1 = 0;}
      if(j==0){
	j_1=N-1;}

         w1[row(i,j)] = u[row(i,j)] - 0.25*dt*d_dx*(u[row(i,j)]*(u[row(ip1,j)]-u[row(i_1,j)])+v[row(i,j)]*(u[row(i,jp1)]-u[row(i,j_1)]) +u[row(ip1,j)]*u[row(ip1,j)]-u[row(i_1,j)]*u[row(i_1,j)]+v[row(i,jp1)]*u[row(i,jp1)]-v[row(i,j_1)]*u[row(i,j_1)])+(dt/p)*f1[row(i,j)];
         w2[row(i,j)] = v[row(i,j)] - 0.25*dt*d_dx*(u[row(i,j)]*(v[row(ip1,j)]-v[row(i_1,j)])+v[row(i,j)]*(v[row(i,jp1)]-v[row(i,j_1)]) +v[row(ip1,j)]*u[row(ip1,j)]-v[row(i_1,j)]*u[row(i_1,j)]+v[row(i,jp1)]*v[row(i,jp1)]-v[row(i,j_1)]*v[row(i,j_1)])+(dt/p)*f2[row(i,j)];
    }}
}  
  
  
int *alloc_1d_int(int n1)  
{  
  int *i;  
      
  i = (int *) malloc(n1 * sizeof(int));  
  if (i == NULL) {  
    printf("Allocation Failure!\n");  
    exit(1);  
  }  
  return i;  
}  
  
  
void free_1d_int(int *i)  
{  
  free(i);  
}  
  
  
double *alloc_1d_double(int n1)  
{  
  double *d;  
      
  d = (double *) malloc(n1 * sizeof(double));  
  if (d == NULL) {  
    printf("Allocation Failure!\n");  
    exit(1);  
  }  
  return d;  
}  
  
  
void free_1d_double(double *d)  
{  
  free(d);  
}  
  
  
int **alloc_2d_int(int n1, int n2)  
{  
  int **ii, *i;  
  int j;  
      
  ii = (int **) malloc(n1 * sizeof(int *));  
  if (ii == NULL) {  
    printf("Allocation Failure!\n");  
    exit(1);  
  }  
  i = (int *) malloc(n1 * n2 * sizeof(int));  
  if (i == NULL) {  
    printf("Allocation Failure!\n");  
    exit(1);  
  }  
  ii[0] = i;  
  for (j = 1; j < n1; j++) {  
    ii[j] = &i[n2 * j];  
  }  
  return ii;  
}  
  
  
void free_2d_int(int **ii)  
{  
  free(ii[0]);  
  free(ii);  
}  
  
  
double **alloc_2d_double(int n1, int n2)  
{  
  double **dd, *d;  
  int j;  
      
  dd = (double **) malloc(n1 * sizeof(double *));  
  if (dd == NULL) {  
    printf("Allocation Failure!\n");  
    exit(1);  
  }  
  d = (double *) malloc(n1 * n2 * sizeof(double));  
  if (d == NULL) {  
    printf("Allocation Failure!\n");  
    exit(1);  
  }  
  dd[0] = d;  
  for (j = 1; j < n1; j++) {  
    dd[j] = &d[n2 * j];  
  }  
  return dd;  
}  
  
  
void free_2d_double(double **dd)  
{  
  free(dd[0]);  
  free(dd);  
}  
        
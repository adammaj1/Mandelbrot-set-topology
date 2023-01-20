/*



  Adam Majewski
 
  
  console program in c programing language 
===============================================================

 
  Structure of a program or how to analyze the program 
  
  
  
  Creating graphic:
  * memory array
  * save it to the disk as a pgm file
  * convert pgm file to png usnigng Image Magic convert
  
 Algorithms
 * interior detection ( speed up computations)
 
 
 
 https://mathr.co.uk/blog/2017-05-17_periodicity_scan.html
  // find the period of a nucleus within a large box uses Robert P. Munafo's Jordan curve method
      int p = m_d_box_period_do(c0, 4.0 * cabs(dc0), maxiters);
      if (p > 0)
        // refine the nucleus location (uses Newton's method)
        if (m_converged == m_d_nucleus(&c0, c0, p, 16))
        {
          // verify the period with a small box
          // if the period is wrong, the size estimates will be way off
          as[atoms].period = m_d_box_period_do(c0, 0.001 * cabs(dc0), 2 * p);
          if (as[atoms].period > 0)
          {
            as[atoms].nucleus = c0;
            // size of component using algorithm from ibiblio.org M-set e-notes
            as[atoms].size = cabs(m_d_size(c0, as[atoms].period));
            // size of atom domain using algorithm from an earlier blog post of mine
            as[atoms].domain_size = m_d_domain_size(c0, as[atoms].period);
            // shape of component (either cardioid or disc) after Dolotin and Morozov (2008 eq. 5.8)
            as[atoms].shape = m_d_shape_discriminant(m_d_shape_estimate(c0, as[atoms].period));
            atoms++; 
  
  

   
  ==========================================

  
  ---------------------------------
  indent m.c 
  default is gnu style 
  -------------------



  c console progam 
  
	export  OMP_DISPLAY_ENV="TRUE"	
  	gcc d.c -lm -Wall -march=native -fopenmp
  	time ./a.out > b.txt


  gcc m.c -lm -Wall -Wextra -march=native -fopenmp


  time ./a.out

  time ./a.out >a.txt
  
  ./g.sh
  
  ============================
  
  gcc m.c -lm -Wall -march=native -fopenmp -pg
  gprof ./a.out > p.txt
  
   

  

*/

#include <stdio.h>
#include <stdlib.h>		// malloc
#include <string.h>		// strcat
#include <math.h>		// M_PI; needs -lm also
#include <complex.h> 		// complex numbers : https://stackoverflow.com/questions/6418807/how-to-work-with-complex-numbers-in-c
#include <stdbool.h>
#include <omp.h>		// OpenMP
#include <time.h>



#define kMax 21 // number of examples, see line 211 plane_examples
#define OrbitLength 50001 // redefine in local_setup





// https://sourceforge.net/p/predef/wiki/Standards/

#if defined(__STDC__)
#define PREDEF_STANDARD_C_1989
#if defined(__STDC_VERSION__)
#if (__STDC_VERSION__ >= 199409L)
#define PREDEF_STANDARD_C_1994
#endif
#if (__STDC_VERSION__ >= 199901L)
#define PREDEF_STANDARD_C_1999
#endif
#endif
#endif





/* --------------------------------- global variables and consts ------------------------------------------------------------ */


// each typedef should have different range !!!



/* Representation FunctionType 
	https://mrob.com/pub/muency/representationfunction.html
	function defining relation between data and the image
*/
typedef enum  {
			LSM =100, 
			LCM = 101,
			DEM = 102, 
			Unknown = 103, 
			BD = 104, 
			MBD = 105, 
			BET = 106,
			Period = 107,
			LastIteration = 108,
			AtomDomains = 109,
			Test = 110,
			PeriodLCM = 111,
			MultiplierMap = 112,  // the best method
			MultiplierMapLCM = 113,
			MultiplierMapAll = 114,
			SAC, 
			DLD, 
			ND, 
			NP, 
			POT, 
			Blend
		
		} RepresentationFunctionTypeT; 


// virtual 2D array and integer ( screen) coordinate
// Indexes of array starts from 0 not 1 
//unsigned int ix, iy; // var
static int ixMin = 0;	// Indexes of array starts from 0 not 1
static int ixMax;	//
static int iWidth;	// horizontal dimension of array

static int iyMin = 0;	// Indexes of array starts from 0 not 1
static int iyMax;	//

static int iHeight = 1000;	//  
// The size of array has to be a positive constant integer 
static int iSize;	// = iWidth*iHeight; 


// ----------memmory 1D arrays ==================
// unsigned char = for 1 byte ( 8 bit) colors 
unsigned char *data;
unsigned char *edge;

 

// unsigned int i; // var = index of 1D array
//static unsigned int iMin = 0; // Indexes of array starts from 0 not 1
static int iMax;	// = i2Dsize-1  = 
// The size of array has to be a positive constant integer 
// unsigned int i1Dsize ; // = i2Dsize  = (iMax -iMin + 1) =  ;  1D array with the same size as 2D array






// on the initial plane , before transformation

double  DisplayAspectRatio = 1.5 ; // https://en.wikipedia.org/wiki/Aspect_ratio_(image)



const complex double critical_point = 	0.0; //




// parameter plane 
double xMin ;	//-0.05;
double xMax ;	//0.75;
double yMin ;	//-0.1;
double yMax ;	//0.7;


double PixelWidth;	// =(CxMax-CxMin)/ixMax;
double PixelHeight;	// =(CyMax-CyMin)/iyMax;

double plane_radius;
complex double plane_center;
double zoom;



// important points
complex double pseudocardioid_nucleus;
complex double pseudocardioid_cusp; // 
complex double pseudocardioid_root_half; // common point between pseudocardioid and period*2 componnet 

complex double second_bulb_nucleus;




/* 

examples containing islands = minibrots
note that plane_center does not equal to hyperbolic component  center
plane radius is adjusted to show only 2 main components of the island
https://mrob.com/pub/muency/mainsequence.html
https://mrob.com/pub/muency/largestislands.html
https://mrob.com/pub/mu-data/largest-islands.txt

plane : plane_center_x,  		plane_center_y,  		plane_radius, 		period
*/
double plane_examples[kMax][4] = {
	{-0.4,				+0.0,				0.8,			1}, 
	{+0.29254,			-0.01497, 			0.00015,		16}, 
	{-1.763,  			+0.0,				0.016,			3}, 
	{-0.15842, 			+1.03335, 			0.008,			4},  
	{+0.358431,			+0.643507,			0.005,			5},  
	{+0.442990,			+0.374,				0.003,			6},  
	{+0.432259,			+0.2275,			0.002,			7}, 
	{+0.404879,			+0.146216,			0.0015,			8}, 
	{+0.378631,			+0.098841,			0.001,			9}, 
	{+0.356854, 			+0.069659,			0.0007,			10},//
	{+0.339454,			+0.050823,			0.0005,			11}, // 10
	{+0.325631,			+0.038164,			0.0005,			12},
	{+0.260270,		 	+0.00167,			0.00002, 		32 },
	{+0.2524945,			+0.0001973,			0.0000025,		64},
	{+0.25061329,			+0.00002399,			0.0000003,		128},
	{+0.250151979,			+0.000002959,			0.000000036,		256},
	{+0.250037823,			+0.0000003673,			0.000000004,		512}, 
	{+0.2500094342,			+0.000000045694786520646,	0.0000000005,		1024}, //17  1024 	2.500094340031833728e-01 + 4.569478652064613658e-08 i 	1.144e-09
	{+0.25000235583,		+5.701985912706845832e-09,	0.00000000007,		2048}, // 
	{+0.2500005886144,		+0.0000000007122,		0.0000000000083,	4096}, // 19  8087162
	{+2.500001471109009610e-01,	+8.897814201389663379e-11,	1.1e-12,		8192}  // 20 , long time 
	 
};


int iterMax_Period ; // see local_setup 
const int iterMax_LastIteration = 100000;
const int iterMax_LSM = 1001;
const int iterMax_DEM = 1001;
const int iterMax_BET = 1000;

// EscapeRadius for bailout test 
double ER = 2.0;
double ER2;	
double ER_DEM = 100.0;
double ER2_DEM; 

//double precision= 1.0E-10; // relate with zoom, PixelWidth in setup
int iPixelSpacingBits ;  // precision in binary digits or bits 
double eps=0.001; // precision
double eps2; 






// pixel counters
int iUnknownTopology = -1; // -1 means that procedure was not used 
int iInterior = 0;
int iExterior = 0;
int iUnknownPeriod = -1;



/* colors = shades of gray from 0 to 255 */
unsigned char iColorOfExterior 		= 255;
unsigned char iColorOfInterior 		= 150;
unsigned char iColorOfPseudoCardioid 	= 150;
unsigned char iColorOfSecondBulb 	= 110;
unsigned char iColorOfBoundary 		= 0;
unsigned char iColorOfUnknown 		= 30;


// to measure time taken by a function 
time_t start_LastIteration;
time_t end_LastIteration;

time_t start_Period;
time_t end_Period;

time_t start_MultiplierMap;
time_t end_MultiplierMap;






// This function swaps values pointed by xp and yp
// https://www.geeksforgeeks.org/c-program-swap-two-numbers/
void swap(double *xp, double *yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}




// Function to find the angle between two lines. Each line is defined by 2 points
double findAngleBetweenLines(const complex double z1a, const complex double z1b, const complex double z2a, const complex double z2b)
{
	// slope of line 1 and 2 
	double m1 = (cimag(z1a) - cimag(z1b))/(creal(z1a) - creal(z1b));
	double m2 = (cimag(z2a) - cimag(z2b))/(creal(z2a) - creal(z2b));
	
	double angle = fabs((m2 - m1)/ (1 + m1 * m2));

	// Calculate tan inverse of the angle
	double ret = atan(angle);

	// Convert the angle from radian to degree
	double val = (ret * 180) / M_PI;

	// return the result
	return val;
}




// the hyperbolic center being the pixel with the smallest multiplier

double complex give_nucleus( const double complex c0, const int period, const int mMax)
{
	double complex c = c0;
	for (int m = 0; m < mMax; ++m)
	{
		double complex z = 0;
		double complex dc = 0;
		for (int i = 0; i < period; ++i)
		{
			dc = 2 * z * dc + 1;
			z = z * z + c;
		}
		c = c - z / dc;
	}
	return c;
}




//  https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/qpolynomials
complex double fc( const double complex z , const complex double c ){

	return z*z +c;
}

// -----------------------------------from screen to world coordinate ; linear mapping ---------------------------------------
// uses global cons
double Give_x (const int ix)
{
	return (xMin + ix * PixelWidth);
}

// uses global cons
double Give_y (const int iy) {
	return (yMax - iy * PixelHeight);  // reverse y axis
}				



// screen to world coordinate mapping
complex double Give_c (const int ix, const int iy)
{
	double x = Give_x (ix);
	double y = Give_y (iy);

	return x + y * I;

}


/* ------------------------------------------ functions -------------------------------------------------------------*/

double clamp(double x, double lo, double hi) {
	return fmin(fmax(x, lo), hi);
}

//------------------complex numbers -----------------------------------------------------




static inline int sgn(double z) {
	if (z > 0) { return  1; }
	if (z < 0) { return -1; }
	return 0;
}

static inline bool odd(int a) {
	return a & 1;
}


static inline bool cisfinite(double _Complex z) {
	return isfinite(creal(z)) && isfinite(cimag(z));
}





 double cabs2(complex double z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}


int SameValue(complex double Z1, complex double Z2)
{
    if (cabs2(Z1- Z2) < eps2 ) 
       {return 1; /* true */ }
       else return 0; /* false */
    }
    

int escapes(complex double z){
	if (cabs2(z) > ER2) 
		return 1; // escapes 
	return 0; // not escapes

}

/* -----------  array functions = drawing -------------- */

/* gives position of 2D point (ix,iy) in 1D array  ; uses also global variable iWidth */
unsigned int Give_i (const int ix, const int iy)
{
  return ix + iy * iWidth;
}


// ***********************************************************************************************
// ********************** Give c from multiplier         ***************************************
// ***************************************************************************************************


/*
https://mathr.co.uk/blog/2013-04-01_interior_coordinates_in_the_mandelbrot_set.html
There is an algorithm to find points on the boundary of the Mandelbrot set, 
given a particular hyperbolic component and the desired internal angle. 
It involves Newton's method in two complex variables to solve
F^p(z,c)=z
∂F^p(z,c)/∂z =b
where 
F0(z,c)=z and 
Fq+1(z,c)=Fq(F(z,c)2+c), 
p is the period of the target component, 
and m=e2πiθ with the θ the desired internal angle. 
The resulting c is the coordinates of the point on the boundary. 
It can also be modified to find points in the interior, simply set m=re2πiθ with |r|≤1.
*/



// mandelbrot-numerics/c/include/mandelbrot-numerics.h


enum m_newton { m_failed, m_stepped, m_converged };
typedef enum m_newton m_newton;

// mandelbrot-numerics/c/bin/m-util.h

// static const double twopi = 6.283185307179586;
// epsilon^2
static const double epsilon2 = 1.9721522630525295e-31;





// mandelbrot-numerics/c/lib/m_d_interior.c 
// double precision: m_d_*()  


m_newton m_d_interior_step(double complex *z_out, double complex *c_out, double complex z_guess, double complex c_guess, double complex multiplier, int period) {
  double complex c = c_guess;
  double complex z = z_guess;
  double complex dz = 1;
  double complex dc = 0;
  double complex dzdz = 0;
  double complex dcdz = 0;
  for (int p = 0; p < period; ++p) {
    dcdz = 2 * (z * dcdz + dc * dz);
    dzdz = 2 * (z * dzdz + dz * dz);
    dc = 2 * z * dc + 1;
    dz = 2 * z * dz;
    z = z * z + c;
  }
  double complex det = (dz - 1) * dcdz - dc * dzdz;
  double complex z_new = z_guess - (dcdz * (z - z_guess) - dc * (dz - multiplier)) / det;
  double complex c_new = c_guess - ((dz - 1) * (dz - multiplier) - dzdz * (z - z_guess)) / det;
  if (cisfinite(z_new) && cisfinite(c_new)) {
    *z_out = z_new;
    *c_out = c_new;
    if (cabs2(z_new - z_guess) <= epsilon2 && cabs2(c_new - c_guess) <= epsilon2) {
      return m_converged;
    } else {
      return m_stepped;
    }
  } else {
    *z_out = z_guess;
    *c_out = c_guess;
    return m_failed;
  }
}


/* 
  usage: ./m_d_interior z-guess c-guess interior period maxsteps
  
  input: 
  	z-guess 	=  	the initial guess z_0 for Newton method, use z = 0
  	c-guess 	= 	the initial guess c_0 for Newton method, use nucleus of given hyperbolic component
  	interior	=  multiplier 
	period	= 	Period of given hyperbolic component
	maxsteps	= 
  
  
*/
m_newton m_d_interior(double complex *z_out, double complex *c_out, double complex z_guess, double complex c_guess, double complex multiplier, int period, int maxsteps) {

	m_newton result = m_failed;
  	double complex z = z_guess;
  	double complex c = c_guess;
  
  	for (int i = 0; i < maxsteps; ++i) {
    		if (m_stepped != (result = m_d_interior_step(&z, &c, z, c, multiplier, period))) 
    			{ break;  }
  		}
  	// 	
  	*z_out = z;
  	*c_out = c;
  	return result;
}




complex double aproximate_c( const int p, const complex double center, const complex double multiplier){


	complex double c = 0.0;
	complex double z = 0;
	int maxsteps = 100;
	
	m_newton result;
	
	result = m_d_interior(&z,  &c, 0.0, center, multiplier, p, maxsteps);
    	if (result != m_converged) 
    		{return -1000;}
	
	
	
	return c;
}







complex double give_c_from_multiplier(const int p, const complex double center, const double angle, const double r )
{
	/*
	input:
	Internal Radius = r in [0,1] ; double
  	Internal Angle In Turns = t  or theta in range [0,1) ; double  
  	p = period ; int 
  	
  	output = c = complex point of 2D parameter plane  
  	*/
  	

	complex double m = 0.0; // multiplier
	complex double c = 0.0; // result 
	double t = angle;
	
	t = t*2*M_PI; // from turns to radians
	m = r* cexp(I*t); // point of unit circle = multiplier
  		
	// map circle to component
	switch (p){
	
		case 1: c = (2.0*m - m*m)/4.0; break;
		case 2: c = (m -4.0)/ 4.0; break;
		default : c = aproximate_c( p, center, m); // for higher periods there is no exact method; use numerical aproximation	 
	
  
	}
	return c; 
}










// ***********************************************************************************************
// ********************** edge detection usung Sobel filter ***************************************
// ***************************************************************************************************

// from Source to Destination
int ComputeBoundaries(const unsigned char S[], unsigned char D[])
{
 
  int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  int i; /* index of 1D array  */
  /* sobel filter */
  unsigned char G, Gh, Gv; 
  // boundaries are in D  array ( global var )
 
  // clear D array
  memset(D, iColorOfExterior, iSize*sizeof(*D)); // 
 
  fprintf(stderr, "\tfind boundaries in S array using  Sobel filter\n");   
#pragma omp parallel for schedule(dynamic) private(i,iY,iX,Gv,Gh,G) shared(iyMax,ixMax)
  for(iY=1;iY<iyMax-1;++iY){ 
    for(iX=1;iX<ixMax-1;++iX){ 
      Gv= S[Give_i(iX-1,iY+1)] + 2*S[Give_i(iX,iY+1)] + S[Give_i(iX-1,iY+1)] - S[Give_i(iX-1,iY-1)] - 2*S[Give_i(iX-1,iY)] - S[Give_i(iX+1,iY-1)];
      Gh= S[Give_i(iX+1,iY+1)] + 2*S[Give_i(iX+1,iY)] + S[Give_i(iX-1,iY-1)] - S[Give_i(iX+1,iY-1)] - 2*S[Give_i(iX-1,iY)] - S[Give_i(iX-1,iY-1)];
      G = sqrt(Gh*Gh + Gv*Gv);
      i= Give_i(iX,iY); /* compute index of 1D array from indices of 2D array */
      if (G==0) {D[i]=255;} /* background */
      else {D[i]=0;}  /* boundary */
    }
  }
 
   
 
  return 0;
}



// copy from Source to Destination
int CopyBoundaries(const unsigned char S[],  unsigned char D[])
{
 
  int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  int i; /* index of 1D array  */
 
 
  fprintf(stderr, "\tcopy boundaries from S array to D array \n");
  for(iY=1;iY<iyMax-1;++iY)
    for(iX=1;iX<ixMax-1;++iX)
      {
		i= Give_i(iX,iY); 
		if(i>-1 && i< iMax)
			{if (S[i]==0) {D[i]=0;} }
			else { fprintf(stderr, "point is out viewport \n");}
      }
 
 
 
  return 0;
}


 
 

// ***************************************************************************************************************************
// ************************** Last Iteration = Fast Iteration =  Interior detection  *****************************************
// ****************************************************************************************************************************


// gives last iterate = escape time
// https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/Mandelbrot_set/mandelbrot
//
 int give_last_iteration(double complex C )
  {
   int i=0;
   int iMax = iterMax_LastIteration;
   
   // note that we start with c instead of 0, to avoid multiplying the derivative by 0
   double complex Z = C; // initial value for iteration Z0 
   complex double D = 1.0; // derivative with respect to z 
   
   
   
   for(i=0;i<iMax;i++)
    { if(cabs2(Z) > ER2) 
    	{ return i;} // exterior
      if(cabs2(D) < eps2) 
      	{return -i; } // interior
      	
      D = 2.0*D*Z; // derivative
      Z = Z*Z+C; // iteration of complex quadratic polynomial:  z = f(z) 
      
    }
    
   // fprintf(stderr,"unknown topology c = %.16f %+.16f \n", creal(c), cimag(c) ); // print for analysis
   iUnknownTopology +=1;  // update pixel counters
   return 0; // unknown
 }
 
 

/* show only 2 xomponent of the island : main pseudocardioid and it's main bulb */
unsigned char ComputeColorOf_last_iteration(complex double c){

	unsigned char color;
	int last_iteration = give_last_iteration(c);
	
	
   		
   	if (last_iteration > 0)
   		{ color = iColorOfExterior; }
   		else {
   			if (last_iteration < 0) 
   				{ color = iColorOfInterior; }
   				else {color = iColorOfUnknown;} // i ==0
		}
	
	return color;


}


 

 





// ***************************************************************************************************************************
// ************************** PER = Period ****************************************************************** ************
// ****************************************************************************************************************************





//**************** Box period *************************************




static double cross(double _Complex a, double _Complex b) {
  return cimag(a) * creal(b) - creal(a) * cimag(b);
}

static bool crosses_positive_real_axis(double _Complex a, double _Complex b) {
  if (sgn(cimag(a)) != sgn(cimag(b))) {
    double _Complex d = b - a;
    int s = sgn(cimag(d));
    int t = sgn(cross(d, a));
    return s == t;
  }
  return false;
}

static bool surrounds_origin(double _Complex a, double _Complex b, double _Complex c, double _Complex d) {
  return odd
    ( crosses_positive_real_axis(a, b)
    + crosses_positive_real_axis(b, c)
    + crosses_positive_real_axis(c, d)
    + crosses_positive_real_axis(d, a)
    );
}

typedef struct  {
  double _Complex c[4];
  double _Complex z[4];
  int p;
} m_d_box_period ;


m_d_box_period *m_d_box_period_new(double _Complex center, double radius) {
  m_d_box_period *box = (m_d_box_period *) malloc(sizeof(*box));
  if (! box) {
    return 0;
  }
  box->z[0] = box->c[0] = center + ((-radius) + I * (-radius));
  box->z[1] = box->c[1] = center + (( radius) + I * (-radius));
  box->z[2] = box->c[2] = center + (( radius) + I * ( radius));
  box->z[3] = box->c[3] = center + ((-radius) + I * ( radius));
  box->p = 1;
  return box;
}

void m_d_box_period_delete(m_d_box_period *box) {
  if (box) {
    free(box);
  }
}

bool m_d_box_period_step(m_d_box_period *box) {
  if (! box) {
    return false;
  }
  bool ok = true;
  for (int i = 0; i < 4; ++i) {
    box->z[i] = box->z[i] * box->z[i] + box->c[i];
    ok = ok && cisfinite(box->z[i]);
  }
  box->p = box->p + 1;
  return ok;
}

bool m_d_box_period_have_period(const m_d_box_period *box) {
  if (! box) {
    return true;
  }
  return surrounds_origin(box->z[0], box->z[1], box->z[2], box->z[3]);
}

int m_d_box_period_get_period(const m_d_box_period *box) {
  if (! box) {
    return 0;
  }
  return box->p;
}

int m_d_box_period_do(double _Complex center, double radius, int maxperiod) {
  m_d_box_period *box = m_d_box_period_new(center, radius);
  if (! box) {
    return 0;
  }
  int period = 0;
  for (int i = 0; i < maxperiod; ++i) {
    if (m_d_box_period_have_period(box)) {
      period = m_d_box_period_get_period(box);
      break;
    }
    if (! m_d_box_period_step(box)) {
      break;
    }
  }
  m_d_box_period_delete(box);
  return period;
}

//**************** period by iteration *************************************

int GivePeriodByIteration (const double complex c ){

	// int period = 0;
	int iMax = iterMax_Period;
	int i;
	complex double orbit[OrbitLength]; // length(orbit) = iMax + 1 ( see local_setup)
	complex double z = 0.0; // critical point
	
	
	
	// iteration without saving points 
	for(i=0; i<iMax; ++i)
    	{ 
    		z  = fc(z, c); 
      		 if (escapes(z) ) {return 0; } // escaping = exterior of M set  so break the procedure
   	} // for(i
	
	
	// iteration = computing the orbit = fiil the array
	 orbit[0] = z; 
  	for(i=1; i<iMax+1; ++i)
    	{ 
    		z  = fc(z, c); 
      		 if ( escapes(z)) {return 0; } // escaping = exterior of M set  so break the procedure
   		orbit[i] = z;
   		//fprintf(stderr," i = %d z = %f+%f \n", i, creal(orbit[i]), cimag(orbit[i]));
      
    	} // for(i=0
	
	// look for similar points = attractor 
	// go from the last point of the orbit 
	//z = orbit[0];
	for(i=iMax-1; i>0; --i){
		if ( SameValue( z, orbit[i]))
			{//printf(" z = %f+%f diff = %e\n", creal(orbit[i]), cimag(orbit[i]), cabs(z - orbit[i])); 
			return iMax - i;} // period
			//else printf(" z = %f+%f diff = %e\n", creal(orbit[i]), cimag(orbit[i]), cabs(z - orbit[i]));
	
	
	
	}
	
	
	
	
	return -1; // period not found 

}



int GivePeriod(complex double c){




	if (cabs2(c)>4.0) {return 0;} // exterior
	if (cabs2(1.0 - csqrt(1.0-4.0*c))<=1.0 ) {return 1;} // main cardioid
	if (cabs2(4.0*c + 4)<=1.0){return 2;} // period 2 component
	
	int period =  GivePeriodByIteration(c); 
	
	if ( period < 0) // last chance
		{ 
			iUnknownPeriod +=1;
			//fprintf(stderr,"unknown period = %d < 0 c = %.16f %+.16f \n", period, creal(c), cimag(c) ); // print for analysis
			//period = m_d_box_period_do(c, 0.5, iterMax_LastIteration);
			//fprintf(stderr,"box period = %d  \n", period);
			
		}
	
	// period > 0 means is periodic 
	// period = 0 means is not periodic = exterior = escaping to infinity
	// period < 0 means period not found, maybe increase global variable iterMax_Period ( see local_setup)

	return period; 

}





/* show only 2 component of the island : main pseudocardioid and it's main bulb */
unsigned char ComputeColorOf_Period(const int period_of_pseudocardioid,  complex double c){

	unsigned char color;
	int period = GivePeriod(c);
	
	
   		
   	if (period == period_of_pseudocardioid)
   		{ color = iColorOfPseudoCardioid ; }
   		else {
   			if (period == 2*period_of_pseudocardioid) 
   				{color= iColorOfSecondBulb; }
   				else { color = iColorOfExterior; }
		}
	
	return color;


}

// ***************************************************************************************************************************
// ************************** MultiplierMap ****************************************************************** ************
// ****************************************************************************************************************************

/* newton function : N(z) = z - (fp(z)-z)/f'(z)) */

complex double N( const complex double c, const complex double zn , const int pMax )
{
//, const double ER2


  
complex double z = zn;
complex double d = 1.0; /* d = first derivative with respect to z */
int p; 

for (p=0; p < pMax; p++){

   //printf ("p = %d ;  z = %f ; %f ;  d = %f ; %f \n", p, creal(z), cimag(z), creal(d), cimag(d)); 
   d = 2*z*d; /* first derivative with respect to z */
   z = z*z +c ; /* complex quadratic polynomial */
   //if (cabs(z) >er2) break;
    
}

 
 
 //printf (" next \n\n"); 
 //if ( cabs2(d) > 2) 
     z = zn - (z - zn)/(d - 1) ;
    
 

return z;
}

/* 
compute periodic point 
of complex quadratic polynomial
using Newton iteration = numerical method 
, double er2

*/

complex double GivePeriodic(complex double c, complex double z0, int period, double eps2){

complex double z = z0;
complex double zPrev = z0; // prevoiuus value of z
int n ; // iteration
const int nMax = 64;

for (n=0; n<nMax; n++) {
     
    z = N( c, z, period); // , er2
    if (cabs2(z - zPrev)< eps2) break;
    
    zPrev = z; }

return z;
}

complex double AproximateMultiplierMap(complex double c, int period, double eps2, double er2){
     
     complex double z;  // variable z 
     complex double zp ; // periodic point 
     complex double zcr = 0.0; // critical point
     complex double d = 1;
     //complex double w;
     
     int p;
     
     zp =  GivePeriodic( c, zcr, period,  eps2); // , er2      Find periodic point z0 such that f^p(z0,c)=z0 using Newton's method in one complex variable
     //zp = find(-1, 0, period, c); 
     //printf (" zp = %f ; %f p = %d \n", creal(zp), cimag(zp), period); 
     
     // Find w by evaluating first derivative with respect to z of f^p at z0 
     if ( cabs2(zp)<er2) {
     
     //printf (" zp = %f ; %f p = %d \n", creal(zp), cimag(zp), period); 
     z = zp;
     for (p=0; p < period; p++){
        d = 2*z*d; /* first derivative with respect to z */
        z = z*z +c ; /* complex quadratic polynomial */
     
          }
          
          
     
     
     }
     else d= 10000;

return d;
}

complex double GiveMultiplierMap( const int period, const complex double c){

	complex double w;
	switch(period){
		case 1: w = 1.0 - csqrt(1.0-4.0*c); 	break; // explicit
		case 2: w = 4.0*c + 4; 			break; //explicit
		default:w = AproximateMultiplierMap(c, period, eps2, ER2); 	break; // 

	}

	return w;

}


unsigned char ComputeColorOf_MultiplierMap(const int period_of_pseudocardioid,  complex double c){

	unsigned char color = 0;
	complex double multiplier = GiveMultiplierMap(period_of_pseudocardioid, c);
	double iRadius = cabs2(multiplier);
	if (iRadius <= 1.0)
		{ // plot pseudocardioid
			color = iColorOfPseudoCardioid;
						
			if (iRadius < 0.005) {
				pseudocardioid_nucleus = give_nucleus(c, period_of_pseudocardioid, 60);
				
			}
			
		} 
		else { // plot second main bulb
			multiplier = GiveMultiplierMap(2*period_of_pseudocardioid, c);
			iRadius = cabs2(multiplier);
			if (iRadius <= 1.0) 
				{ 
					color = iColorOfSecondBulb ; 
					if (iRadius < 0.005) {
						second_bulb_nucleus = give_nucleus(c, 2*period_of_pseudocardioid, 60);
											
						}
					
				}
				else {color = iColorOfExterior;}
				
				
				
			
		} 
		
		
	
	
	return color;
	
}

/* ==================================================================================================
 ============================= Draw functions ===============================================================
 =====================================================================================================
*/ 
unsigned char ComputeColor(const int period_of_pseudocardioid, const RepresentationFunctionTypeT RepresentationFunctionType, const complex double c ){

	unsigned char iColor= 0;
	
	switch(RepresentationFunctionType){
	
		case LastIteration : 	{iColor = ComputeColorOf_last_iteration(c); 				break;}
		case Period: 		{iColor = ComputeColorOf_Period(period_of_pseudocardioid,  c); 		break;}
		case MultiplierMap: 	{iColor = ComputeColorOf_MultiplierMap(period_of_pseudocardioid,  c); 	break;}
		
		default: {}
	}
	return iColor;



}



unsigned char  GiveColor(const int period, const RepresentationFunctionTypeT RepresentationFunctionType, const int ix, const int iy){

	complex double c = Give_c(ix,iy);
	unsigned char iColor = ComputeColor(period, RepresentationFunctionType, c);
	
	return iColor; 

}


// plots  raster point (ix,iy) = computes it's color and save it to the array A
int DrawPoint (const int period, const RepresentationFunctionTypeT RepresentationFunctionType,  const int ix, const int iy, unsigned char A[])
{
	
	
	unsigned char iColor = GiveColor(period,  RepresentationFunctionType,  ix, iy);
	int i = Give_i (ix, iy);	/* compute index of 1D array from indices of 2D array */
    	A[i] = iColor ; // 
  
  	return 0;
}




// fill array 
// uses global var :  ...
// scanning complex plane 
int DrawImage (const int k, const RepresentationFunctionTypeT RepresentationFunctionType,   unsigned char A[])
{
  	int ix, iy;		// pixel coordinate 
  	
  	const int period = (int) plane_examples[k][3];
  	// start computing unknown pixels , if iUnknown <0 it means that it was not used !
  	//if (RepresentationFunctionType==LastIteration ) {iUnknownTopology = 0;}
  	//if (RepresentationFunctionType==Period ) 	{iUnknownPeriod = 0;} 

  	fprintf(stderr, "\tRepresentationFunctionType = %d period = %d  \n", RepresentationFunctionType, period);
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic,1) private(ix,iy) shared(A, ixMax , iyMax)
	// #pragma omp parallel for schedule(dynamic, 1)
  	for (iy = iyMin; iy <= iyMax; ++iy){
    		fprintf (stderr, "\t\t%d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix)
      			{DrawPoint(period, RepresentationFunctionType,   ix, iy, A);}	//  
  		}
  	fprintf (stderr, "\n");	//info 

  return 0;
}



int IsInsideCircle (const int x, const int y, const int xcenter, const int ycenter, const int r){

	
  double dx = x- xcenter;
  double dy = y - ycenter;
  double d = sqrt(dx*dx+dy*dy);
  if (d<r) {    return 1;}
  return 0;
	  

} 


int PlotBigPoint(const complex double c, const unsigned char iColor, const double p_size, unsigned char A[]){

	
  int ix_seed = (creal(c)-xMin)/PixelWidth;
  int iy_seed = (yMax - cimag(c))/PixelHeight;
  int i;
	
	
  /* mark seed point by big pixel */
  int iSide =p_size*iWidth/2000.0 ; /* half of width or height of big pixel */
  int iY;
  int iX;
  for(iY=iy_seed-iSide;iY<=iy_seed+iSide;++iY){ 
    for(iX=ix_seed-iSide;iX<=ix_seed+iSide;++iX){ 
      if (IsInsideCircle(iX, iY, ix_seed, iy_seed, iSide)) {
	i= Give_i(iX,iY); /* index of _data array */
	if(i>-1 && i< iMax)
		{A[i] = iColor;}
		else { fprintf(stderr, "point is out viewport \n");}
      }
      // else {printf(" bad point \n");}
	
    }}
	
	
  return 0;
	
}



// plots raster point (ix,iy) 
int iDrawPoint( int ix, int iy, unsigned char iColor, unsigned char A[]){ 

  /* i =  Give_i(ix,iy) compute index of 1D array from indices of 2D array */
	if (ix >=ixMin && ix<=ixMax && iy >=iyMin && iy<=iyMax )
		{A[Give_i(ix,iy)] = iColor;}
		else {fprintf (stdout,"iDrawPoint :   (%d; %d) is outside\n", ix,iy); return 1; }

	return 0;
}




//------------------------- draw arc -----------------------------------------------

void dDrawPoint(const complex double c, unsigned char color,  unsigned char A[]){

	int ix, iy; // screen coordinate = indices of virtual 2D array 
	ix = (creal(c)- xMin)/PixelWidth; 
  	iy = (yMax - cimag(c))/PixelHeight; // inverse Y axis 
  	iDrawPoint( ix, iy, color, A);
}





// y = b + r sin t
// sin accepts angle expressed in radians 
void dDrawCirclePoint(const complex double circle_center, const double angle, const double radius, unsigned char color, unsigned char A[] ){

	double x = creal(circle_center) + radius * cos(angle*2.0*M_PI); // x = a + r cos t, angle is in turns, is converted to radian
	double y = cimag(circle_center) + radius * sin(angle*2.0*M_PI);
	complex double c = x + y*I;
	dDrawPoint( c, color,  A);
	
}




// angles are measured in turns so t is in [0,1] range
// 0.0 <= start_angle < end_angle <= 1.0

int DrawArc(const complex double circle_center, const double start_angle, const double end_angle, const double radius, unsigned char color, unsigned char A[]){

	
	
	double dx = xMax- xMin;
	
	int iMax = 10.0*iWidth * radius/dx;
	double dt = 1.0/iMax; // angle step ; 
	double angle;
	
	
	if (radius > dx) {fprintf(stderr," to big radius = %.16f > image_width = %.16f\n", radius, dx); return 1;}
	
	angle = start_angle;
	while (angle <= end_angle){
		dDrawCirclePoint( circle_center, angle, radius, color, A);
		dDrawCirclePoint( circle_center, angle, radius + PixelWidth, color, A); // thick curve by lazy programmer
		dDrawCirclePoint( circle_center, angle, radius - PixelWidth, color, A);
		angle += dt;
	} 
	

	return 0;

}








// -----------------------------------------------------------------------------------------
// ------------------- draw segment of the line ---------------------------------------------
//-------------------------------------------------------------------------------------------







/*
  http://rosettacode.org/wiki/Bitmap/Bresenham%27s_line_algorithm
  Instead of swaps in the initialisation use error calculation for both directions x and y simultaneously:
*/
void iDrawLine( int x0, int y0, int x1, int y1, unsigned char iColor, unsigned char A[]) 
{
  int x=x0; int y=y0;
  int dx = abs(x1-x0), sx = x0<x1 ? 1 : -1;
  int dy = abs(y1-y0), sy = y0<y1 ? 1 : -1; 
  int err = (dx>dy ? dx : -dy)/2, e2;

  for(;;){
    iDrawPoint(x, y, iColor, A);
    if (x==x1 && y==y1) break;
    e2 = err;
    if (e2 >-dx) { err -= dy; x += sx; }
    if (e2 < dy) { err += dx; y += sy; }
  }
}




int dDrawLineSegment(const double complex Z0, const double complex Z1, const unsigned char color, unsigned char A[]){

  double Zx0 = creal(Z0);
  double Zy0 = cimag(Z0);
  double Zx1 = creal(Z1);
  double Zy1 = cimag(Z1);
  int ix0, iy0; // screen coordinate = indices of virtual 2D array 
  int ix1, iy1; // screen coordinate = indices of virtual 2D array

  // first step of clipping
  //if (  Zx0 < ZxMax &&  Zx0 > ZxMin && Zy0 > ZyMin && Zy0 <ZyMax 
  // && Zx1 < ZxMax &&  Zx1 > ZxMin && Zy1 > ZyMin && Zy1 <ZyMax )
   	
  ix0= (Zx0- xMin)/PixelWidth; 
  iy0 = (yMax - Zy0)/PixelHeight; // inverse Y axis 
  ix1= (Zx1- xMin)/PixelWidth; 
  iy1= (yMax - Zy1)/PixelHeight; // inverse Y axis 
   	
  // second step of clipping
  if (ix0 >=ixMin && ix0<=ixMax && ix0 >=ixMin && ix0<=ixMax && iy0 >=iyMin && iy0<=iyMax && iy1 >=iyMin && iy1<=iyMax )
    iDrawLine(ix0,iy0,ix1,iy1,color, A) ;

  return 0;
}


int dDrawThickLineSegment(const double complex Z0, const double complex Z1, const unsigned char color, unsigned char A[]){

	dDrawLineSegment(Z0, 			Z1, 			color,  A);
	dDrawLineSegment(Z0+PixelWidth, 	Z1+PixelWidth, 		color,  A); // thic line by lazy programmer
	dDrawLineSegment(Z0-PixelWidth, 	Z1-PixelWidth, 		color,  A);
	dDrawLineSegment(Z0+PixelHeight*I, 	Z1+PixelHeight*I, 	color,  A);
	dDrawLineSegment(Z0-PixelHeight*I, 	Z1-PixelHeight*I, 	color,  A);
	
  return 0;
}

int dDrawThickExtendedLineAndArcSegment(const double complex Z0, const double complex Z1, const double complex Z2, const unsigned char color, unsigned char A[]){

	
		
	/* 
		new point on the line but Z0 is a center of a circle, Z1 is a point on the circle
		circle equation in Parametric form  
		x = a + r*cos(t)
		y = b + r*sin (t)
		where t is a parametric variable in the range 0 to 2pi, interpreted geometrically as the angle]
	*/
	double radius = cabs( Z0 - Z1);
	//double val = 180.0 / M_PI;

	double angle1 = M_PI + atan2(cimag(Z0) - cimag(Z1), creal(Z0) - creal(Z1)); //  radians in range [0,2 pi]
	double angle2 = M_PI + atan2(cimag(Z1) - cimag(Z2), creal(Z1) - creal(Z2)); //  radians in range [0,2 pi]
	
	
	
	double dr = PixelWidth*iWidth/20;
	double x = creal(Z0) +  (radius + dr) *  cos(angle1);
	double y = cimag(Z0) +  (radius + dr) *  sin(angle1);
	complex double Z1new = x+y*I;
	
	// draw thick line
	dDrawThickLineSegment(Z0, Z1new, color, A);
	// Arc
	
	if (angle2 < angle1) {swap( &angle1, &angle2);} // islands are rotating
	angle1 = angle1/(2.0*M_PI); // from radians to turns
	angle2 = angle2/(2.0*M_PI); // from radians to turns
	DrawArc( Z1, angle1, angle2, dr, color, A); // show distorion angle
	
  return 0;
}





 
// *******************************************************************************************
// ********************************** save A array to pgm file ****************************
// *********************************************************************************************

int SaveImage(const unsigned char A[], const char *shortName )
{

  FILE *fp;
  const int MaxColorComponentValue = 255;	/* color component is coded from 0 to 255 ;  it is 8 bit color file */
  
  
  
  // https://programmerfish.com/create-output-file-names-using-a-variable-in-c-c/
  char fileName[512];
  const char* fileType = ".pgm";
  sprintf(fileName,"%s%s", shortName, fileType); // 
  
  
  
  char long_comment[200];
  sprintf (long_comment, "one parameter family of complex quadratic polynomial, parameter plane ");





  // save image array to the pgm file 
  fp = fopen (fileName, "wb");	// create new file,give it a name and open it in binary mode 
  fprintf (fp, "P5\n # %s\n %d %d\n %d\n", long_comment, iWidth, iHeight, MaxColorComponentValue);	// write header to the file
  size_t rSize = fwrite (A, sizeof(A[0]), iSize, fp);	// write whole array with image data bytes to the file in one step 
  fclose (fp);

  // info 
  if ( rSize == (long unsigned int) iSize) 
  	{
  		fprintf (stdout, "\tFile %s saved ", fileName);
  		if (long_comment == NULL || strlen (long_comment) == 0)
    		printf ("\n");
  			else { fprintf (stdout,". Comment = %s \n", long_comment); }
  	}
  	else {fprintf(stdout, "wrote %zu elements out of %u requested\n", rSize,  iSize);}
  	
  	
   
  
  return 0;
}


// *****************************************************************************
//;;;;;;;;;;;;;;;;;;;;;;  program setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// **************************************************************************************

// uses global var  
// local  set up for every example from plane_examples array
int local_setup(const int k){

	fprintf(stderr, "Example island = %d \n", k);
  	// read from the array with presettings
  	plane_radius = plane_examples[k][2];
  	plane_center = plane_examples[k][0] + I*plane_examples[k][1];
  	
  	yMin = cimag(plane_center) - plane_radius;	// 
  	yMax = cimag(plane_center) + plane_radius;	//
  	
  	xMin = creal(plane_center) - plane_radius*DisplayAspectRatio;
  	xMax = creal(plane_center) + plane_radius*DisplayAspectRatio;
  	
  	/* Pixel sizes of the plane */
  	PixelWidth = (xMax - xMin) / ixMax;	//  ixMax = (iWidth-1)  step between pixels in world coordinate 
	PixelHeight = (yMax - yMin) / iyMax;
  	
  	iPixelSpacingBits = -log2( PixelWidth); // 
	
	
	eps = PixelWidth/100.0; // to see detailes smaller then pixel 
	zoom = 1.0/plane_radius;
	
	
	// for cabs2
	ER2 = ER*ER;
	ER2_DEM = ER_DEM*ER_DEM;
	eps2 = eps*eps; // precision
	
	//  if iUnknown <0 it means that it was not used !
	iUnknownTopology = -1; // 
	iUnknownPeriod = -1;
	
	int period = (int) plane_examples[k][3];
	iterMax_Period = period*3000;
	
	
	/* Undefine before redefining : C does not support any additional directive to redefine an existing macro. 
	You can define the same macro any number of times. 
	However. doing so will populate your compiler warning stack. Hence it is always advisable to first undefine an existing macro then define it again. */
	#ifdef OrbitLength
	#undef OrbitLength
	#endif
	#define OrbitLength (iterMax_Period+1) // = iterMax_Period+1
	
	

  	
  	return 0;

}



// global setup = the same for all pixels 
int global_setup()
{

	fprintf (stderr, "global setup start\t");

  /* 2D array ranges */
  iWidth = iHeight * DisplayAspectRatio; 
  iSize = iWidth * iHeight;	// size = number of points in array 
  
  iyMax = iHeight - 1;		// Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].
  ixMax = iWidth - 1;

  /* 1D array ranges */
  // i1Dsize = i2Dsize; // 1D array with the same size as 2D array
  iMax = iSize - 1;		// Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].
  
  /* create dynamic 1D arrays for colors ( shades of gray ) */
  data = malloc (iSize * sizeof (unsigned char));
  edge = malloc (iSize * sizeof (unsigned char));
   
  //
 
  	
  if (data == NULL || edge == NULL ){
    fprintf (stderr, " Setup error : Could not allocate memory\n");
    return 1;
  }

	fprintf (stderr,"\tand ends  \n");
	
  return 0;

} // ;;;;;;;;;;;;;;;;;;;;;;;;; end of the setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





const char* GiveName(const int k, const RepresentationFunctionTypeT RepresentationFunctionType)
{
	static char Name[512];
	// NoOfExample_RepresentationFunctionType                   period_p_sName
    	sprintf(Name,"%d_%d", k, RepresentationFunctionType) ; //(int) plane_examples[k][3], sName);

    	return Name;
}


// *****************************************************************************
//;;;;;;;;;;;;;;;;;;;;;;  print info  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// **************************************************************************************





int PrintInfoAboutPixelSize(int pixel_spacing_bits)
{
	fprintf(stdout, "\tpixel_spacing_bits = %d \t", pixel_spacing_bits );

	if (pixel_spacing_bits <= 40) { fprintf(stdout, "use float \n"); return 0;}
	if (pixel_spacing_bits > 60)  
		{fprintf(stdout, "use arbitrary precision number, like MPFR \n");
      return 0; }
  
  if (pixel_spacing_bits > 50)   
       {fprintf(stdout, "use long double \n");
        return 0; }

  if (pixel_spacing_bits > 40)  
       fprintf(stdout, "use double \n");
	return 0;
}


void PrintInfoAboutImage(const int k){

	fprintf(stdout, " \n");
	fprintf(stdout, "This island (=  mini Mandelbrot set) is %d example from plane_examples \n",k);
	fprintf(stdout, "\tPeriod of main  (pseudo)cardioid is %d \n", (int) plane_examples[k][3]); 
	fprintf(stdout, "\tpseudocardioid_nucleus = %.16f %+0.16f\n",creal(pseudocardioid_nucleus), cimag(pseudocardioid_nucleus));
	fprintf(stdout, "\tpseudocardioid_cusp = %.16f %+0.16f\n",creal(pseudocardioid_cusp), cimag(pseudocardioid_cusp));
	fprintf(stdout, "\tpseudocardioid_root_half = %.16f %+0.16f\n",creal(pseudocardioid_root_half), cimag(pseudocardioid_root_half));
	fprintf(stdout, "\tsecond_bulb_nucleus = %.16f %+0.16f\n",creal(second_bulb_nucleus), cimag(second_bulb_nucleus));
	// 
	double distortion = findAngleBetweenLines(pseudocardioid_nucleus, pseudocardioid_cusp, pseudocardioid_cusp, pseudocardioid_root_half);
	fprintf(stdout, "\tisland distortion = %.16f degrees = %.16f turns = %.16f radians\n", distortion, distortion/ 360.0 , distortion/0.017453292519);
	double island_size = cabs(pseudocardioid_cusp - pseudocardioid_root_half);
	fprintf(stdout, "\tisland sizes: HalfRoot-Cusp  = %.16f  = %.16f of Image Width \n", island_size, island_size/ (xMax-xMin));
	fprintf(stdout, " \n");
	fprintf(stdout, " \n");
		fprintf(stdout, "plane description\n");
	fprintf(stdout, "\tplane center = %.16f %+0.16f\n",creal(plane_center), cimag(plane_center));
	fprintf(stdout, "\tzoom = %.16e\n",zoom);
	fprintf(stdout, "\tplane radius = %.16e\n",plane_radius);
	fprintf(stdout, "\txMin = %.16e \t xMax = %.16e\n",xMin, xMax); // 
	fprintf(stdout, "\tyMin = %.16e \t yMax = %.16e\n",yMin, yMax); 
	fprintf(stdout, "\tPixelWidth = %.16e   \n", PixelWidth); 
	PrintInfoAboutPixelSize(iPixelSpacingBits);
	fprintf(stdout, " \n");
	fprintf(stdout, "Important parameters:\n");
	fprintf(stdout, "\teps = %.16e\n",eps);
	fprintf(stdout, "\tEscape Radius = ER = %.16e\n", ER);
	fprintf(stdout, "\titerMax_LastIteration = %d\n", iterMax_LastIteration);
	fprintf(stdout, "\titerMax_Period = %d\n", iterMax_Period);// iterMax_Period %d\n", RepresentationFunctionType, period,iterMax_Period);
	fprintf(stdout, " \n");
	fprintf(stdout, " \n");
	
	
	printf("==========================================================================================================================\n\n\n\n");
	fprintf(stdout, " \n");
};


int PrintInfoAboutProgam(void)
{
	printf (" \n");
	printf("gcc version: %d.%d.%d\n",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__); // https://stackoverflow.com/questions/20389193/how-do-i-check-my-gcc-c-compiler-version-for-my-eclipse
	// OpenMP version is displayed in the console 
  return 0;
}


// *****************************************************************************
//;;;;;;;;;;;;;;;;;;;;;;  Make Images  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// **************************************************************************************


int MakeImagesForExample( const int k){

/* 
LSM =100, 
LCM = 101,
DEM = 102, 
Unknown = 103, 
BD = 104, 
MBD = 105, 
BET = 106,
Period = 107,
LastIteration = 108,
AtomDomains = 109,
Test = 110,
PeriodLCM = 111,
*/

	const char *Name;
	
	local_setup(k); // for every k 
	
	
	/*
	
	// ********* LastIteration **********************
	time (&start_LastIteration);
	DrawImage(k, LastIteration,   data); // all the periods  = 108
	time (&end_LastIteration);
	
	Name = GiveName(k, LastIteration);
	SaveImage(data, Name); 

	
	ComputeBoundaries(data,edge); // boundaries of LastIteration = 101
	Name = GiveName(k, LCM);
	SaveImage(edge, Name); 
	
	
	CopyBoundaries(edge, data);
	Name = GiveName(k, "LastIteration");
	SaveImage(data, Name); 
	
	
	DrawImage(k, LSM,   data);
	Name = GiveName(k, LSM, "LSM");
	SaveImage(data, Name); 
	
	
	
	DrawImage(k, AtomDomains, data); // all the periods 
	Name = GiveName(k, "AtomDomains");
	SaveImage(data, Name); 
	
	
	
	// ********* period **********************
	time (&start_Period);
	DrawImage(k, Period, data); // only 2 main components = 107
	time (&end_Period);
	
	Name = GiveName(k, Period);
	SaveImage(data, Name); 
	
	ComputeBoundaries(data,edge); // boundaries of LastIteration = 101
	Name = GiveName(k, LCM);
	SaveImage(edge, Name); 
	
	*/
	
	
	
	// ********* MultiplierMap **********************
	//time (&start_MultiplierMap);
	DrawImage(k, MultiplierMap, data); // only 2 main components = 107
	//time (&end_MultiplierMap);
	
	Name = GiveName(k, MultiplierMap);
	//SaveImage(data, Name); 
	
	
	
	
	ComputeBoundaries(data,edge); // boundaries of LastIteration = 101
	//Name = GiveName(k, LCM);
	//SaveImage(edge, Name); 
	
	CopyBoundaries(edge, data); // 111
	//Name = GiveName(k, MultiplierMapLCM);
	//SaveImage(data, Name); 
	
	// 
	pseudocardioid_cusp = give_c_from_multiplier((int) plane_examples[k][3], pseudocardioid_nucleus, 0.0, 1.0 );
	pseudocardioid_root_half = give_c_from_multiplier((int) plane_examples[k][3], pseudocardioid_nucleus, 0.5, 1.0 );
	
	PlotBigPoint(pseudocardioid_nucleus, iColorOfBoundary, 10.0, data);
	PlotBigPoint(pseudocardioid_cusp, iColorOfBoundary, 10.0, data);
	PlotBigPoint(pseudocardioid_root_half, iColorOfBoundary, 10.0, data);
	PlotBigPoint(second_bulb_nucleus, iColorOfBoundary, 10.0, data);
	
	dDrawThickExtendedLineAndArcSegment(pseudocardioid_cusp, pseudocardioid_nucleus, pseudocardioid_root_half, iColorOfBoundary, data);
	
	dDrawThickLineSegment(pseudocardioid_nucleus, pseudocardioid_root_half, iColorOfBoundary, data);
	
	
	
	Name = GiveName(k, MultiplierMapAll);
	SaveImage(data, Name); 

	
	PrintInfoAboutImage(k);
	
	return 0;

}








void end(void){


  fprintf (stderr," allways free memory (deallocate )  to avoid memory leaks \n"); // https://en.wikipedia.org/wiki/C_dynamic_memory_allocation
  free (data);
  free(edge);
  
  PrintInfoAboutProgam();
}









// ********************************************************************************************************************
/* -----------------------------------------  main   -------------------------------------------------------------*/
// ********************************************************************************************************************

int main () {
  
  	int k = 0; // example number
	
	global_setup();
	
	for (k=0; k<kMax; ++k) { // kMax
		MakeImagesForExample(k);
	}
		
  	end();
  	
  	
  	return 0;
}

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
#include <omp.h>		// OpenMP


#define kMax 12 // number of examples, see line 211 plane_examples

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
static unsigned int ixMin = 0;	// Indexes of array starts from 0 not 1
static unsigned int ixMax;	//
static unsigned int iWidth;	// horizontal dimension of array

static unsigned int iyMin = 0;	// Indexes of array starts from 0 not 1
static unsigned int iyMax;	//

static unsigned int iHeight = 4000;	//  
// The size of array has to be a positive constant integer 
static unsigned int iSize;	// = iWidth*iHeight; 


// ----------memmory 1D arrays ==================
// unsigned char = for 1 byte ( 8 bit) colors 
unsigned char *data;
unsigned char *edge;

 

// unsigned int i; // var = index of 1D array
//static unsigned int iMin = 0; // Indexes of array starts from 0 not 1
static unsigned int iMax;	// = i2Dsize-1  = 
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


complex double pseudo_cardioid_center; // nucleus
complex double pseudo_cardioid_cusp; // 
complex double pseudo_cardioid_root_half; // common point between pseudocardioid and period*2 componnet 

/* 
plane : plane_center_x,  plane_center_y,  plane_radius, period
examples containing islands = minibrots
note that plane_center does not equal to hyperbolic component  center
plane radius is adjusted to show only 2 main components of the island
https://mrob.com/pub/muency/mainsequence.html
https://mrob.com/pub/muency/largestislands.html
https://mrob.com/pub/mu-data/largest-islands.txt
*/
double plane_examples[kMax][4] = {
	{-0.4,		+0.0,		0.8,		1}, 
	{+0.2925755,	-0.0149977, 	0.00025,	32}, 
	{-1.763,  	+0.0,		0.016,		3}, 
	{-0.15842, 	+1.03335, 	0.01,		4},  
	{+0.358431,	+ 0.643507,	0.006,		5},  
	{+0.442990,	+0.373727,	0.005,		6}, 
	{+0.432259,	+0.227315,	0.003,		7}, 
	{+0.404879,	+0.146216,	0.002,		8}, 
	{+0.378631,	+0.098841,	0.001,		9}, 
	{+0.356854, 	+0.069659,	0.001,		10},
	{+0.339454,	+0.050823,	0.001,		11},
	{+0.325631,	+0.038164,	0.001,		12}
	 
};



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
int iUnknown = 0;
int iInterior = 0;
int iExterior = 0;




/* colors = shades of gray from 0 to 255 */
unsigned char iColorOfExterior 		= 255;
unsigned char iColorOfInterior 		= 150;
unsigned char iColorOfPseudoCardioid 	= 120;
unsigned char iColorOfMainBulb 		= 150;
unsigned char iColorOfBoundary 		= 0;
unsigned char iColorOfUnknown 		= 30;


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
   iUnknown +=1;  // update pixel counters
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


int GivePeriod (const double complex c ){

	// int period = 0;
	int iMax = 10000;
	int i;
	complex double orbit[10001]; // length(orbit) = iMax + 1
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
   		//printf(" i = %d z = %f+%f \n", i, creal(orbit[i]), cimag(orbit[i]));
      
    	} // for(i=0
	
	// look for similar points = attractor 
	// go from the last point of the orbit 
	//z = orbit[0];
	for(i=iMax-1; i>0; --i){
		if ( SameValue(z, orbit[i]) ) { return iMax - i;} // period
			
	
	
	
	}
	
	
	
	
	return -1; // period not found 

}


/* show only 2 xomponent of the island : main pseudocardioid and it's main bulb */
unsigned char ComputeColorOf_Period(const int period_of_pseudocardioid,  complex double c){

	unsigned char color;
	int period = GivePeriod(c);
	
	
   		
   	if (period == period_of_pseudocardioid)
   		{ color = iColorOfPseudoCardioid ; }
   		else {
   			if (period == 2*period_of_pseudocardioid) 
   				{color= iColorOfMainBulb; }
   				else { color = iColorOfExterior; }
		}
	
	return color;


}




/* -----------  array functions = drawing -------------- */

/* gives position of 2D point (ix,iy) in 1D array  ; uses also global variable iWidth */
unsigned int Give_i (const int ix, const int iy)
{
  return ix + iy * iWidth;
}


// ***********************************************************************************************
// ********************** edge detection usung Sobel filter ***************************************
// ***************************************************************************************************

// from Source to Destination
int ComputeBoundaries(const unsigned char S[], unsigned char D[])
{
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
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
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
 
 
  fprintf(stderr, "\tcopy boundaries from S array to D array \n");
  for(iY=1;iY<iyMax-1;++iY)
    for(iX=1;iX<ixMax-1;++iX)
      {i= Give_i(iX,iY); if (S[i]==0) D[i]=0;}
 
 
 
  return 0;
}


 
 
/* ==================================================================================================
 ============================= Draw functions ===============================================================
 =====================================================================================================
*/ 
unsigned char ComputeColor(const int period_of_pseudocardioid, const RepresentationFunctionTypeT RepresentationFunctionType, const complex double c ){

	unsigned char iColor= 0;
	
	
	
	switch(RepresentationFunctionType){
	
		
		
		
		case LastIteration : 	{iColor = ComputeColorOf_last_iteration(c); 			break;}
		case Period: 		{iColor = ComputeColorOf_Period(period_of_pseudocardioid,  c); 	break;}
		
		
	
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
	unsigned int i = Give_i (ix, iy);	/* compute index of 1D array from indices of 2D array */
    	A[i] = iColor ; // 
  
  	return 0;
}




// fill array 
// uses global var :  ...
// scanning complex plane 
int DrawImage (const int k, const RepresentationFunctionTypeT RepresentationFunctionType,   unsigned char A[])
{
  	unsigned int ix, iy;		// pixel coordinate 
  	
  	const int period = (int) plane_examples[k][3];
  	
  	

  	fprintf(stderr, "compute image %d RepresentationFunctionType = %d \t  \n", k, RepresentationFunctionType);
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic,1) private(ix,iy) shared(A, ixMax , iyMax)
	// #pragma omp parallel for schedule(dynamic, 1)
  	for (iy = iyMin; iy <= iyMax; ++iy){
    		fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix)
      			{DrawPoint(period, RepresentationFunctionType,   ix, iy, A);}	//  
  		}

  return 0;
}













 
// *******************************************************************************************
// ********************************** save A array to pgm file ****************************
// *********************************************************************************************

int SaveImage(const unsigned char A[], const char *shortName )
{

  FILE *fp;
  const unsigned int MaxColorComponentValue = 255;	/* color component is coded from 0 to 255 ;  it is 8 bit color file */
  
  
  
  // https://programmerfish.com/create-output-file-names-using-a-variable-in-c-c/
  char fileName[512];
  const char* fileType = ".pgm";
  sprintf(fileName,"%s%s", shortName, fileType); // 
  
  
  
  char long_comment[200];
  sprintf (long_comment, "one parameter family of complex quadratic polynomial, parameter plane ");





  // save image array to the pgm file 
  fp = fopen (fileName, "wb");	// create new file,give it a name and open it in binary mode 
  fprintf (fp, "P5\n # %s\n %u %u\n %u\n", long_comment, iWidth, iHeight, MaxColorComponentValue);	// write header to the file
  size_t rSize = fwrite (A, sizeof(A[0]), iSize, fp);	// write whole array with image data bytes to the file in one step 
  fclose (fp);

  // info 
  if ( rSize == iSize) 
  	{
  		printf ("File %s saved ", fileName);
  		if (long_comment == NULL || strlen (long_comment) == 0)
    		printf ("\n");
  			else { printf (". Comment = %s \n", long_comment); }
  	}
  	else {printf("wrote %zu elements out of %u requested\n", rSize,  iSize);}
  	
  	
   
  
  return 0;
}






complex double  GiveCenter ()
{
  	complex double local_center = 0.0;

  	
  return local_center;
}





/*

********************************************* info 

*/





// *****************************************************************************
//;;;;;;;;;;;;;;;;;;;;;;  program setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// **************************************************************************************

// uses global var  
// local  set up for every example from plane_examples array
int local_setup(const int k){

	
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
  	
  	
  	eps = PixelWidth/1000.0; // to see detailes smaller then pixel 
  	zoom = 1.0/plane_radius;
  	
  	
  	// for cabs2
  	ER2 = ER*ER;
   	ER2_DEM = ER_DEM*ER_DEM;
  	eps2 = eps*eps; // precision
  	
  	iUnknown = 0; // 
  	
  	return 0;

}



// globa; setup = the same for all pixels 
int setup()
{

	fprintf (stderr, "setup start\n");

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
    fprintf (stderr, " Setup error : Could not allocate memory");
    return 1;
  }

	fprintf (stderr,"\tend of setup \n");
	
  return 0;

} // ;;;;;;;;;;;;;;;;;;;;;;;;; end of the setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





const char* GiveName(const int k,  const char* sName)
{
	static char Name[512];
	// NoOfExample_period_p_sName
    	sprintf(Name,"%d_period_%d_%s", k, (int) plane_examples[k][3], sName);

    	return Name;
}


// *****************************************************************************
//;;;;;;;;;;;;;;;;;;;;;;  print info  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// **************************************************************************************


int PrintInfoAboutProgam(void)
{
	printf (" \n");
	printf("gcc version: %d.%d.%d\n",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__); // https://stackoverflow.com/questions/20389193/how-do-i-check-my-gcc-c-compiler-version-for-my-eclipse
	// OpenMP version is displayed in the console 
  return 0;
}



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
	//fprintf(stdout, "pseudo_cardioid_center = %.16f %+0.16f\n",creal(pseudo_cardioid_center), cimag(pseudo_cardioid_center));
	fprintf(stdout, "\t unknown pixels = %d pixels = %.16f of all pixels\n", iUnknown, ((double) iUnknown)/iSize);
	fprintf(stdout, "\t No unknown pixels ( zero) is a perfect target, but less then 1 %% of all pixels ( ratio < 0.01 ) is also OK ( acceptable)\n");
	//
	fprintf(stdout, "plane description\n");
	fprintf(stdout, "\tplane center = %.16f %+0.16f\n",creal(plane_center), cimag(plane_center));
	fprintf(stdout, "\tzoom = %.16e\n",zoom);
	fprintf(stdout, "\tplane radius = %.16e\n",plane_radius);
	fprintf(stdout, "\tPixelWidth = %.16e   \n", PixelWidth); 
	fprintf(stdout, "\txMin = %.16e \t xMax = %.16e\n",xMin, xMax); // 
	fprintf(stdout, "\tyMin = %.16e \t yMax = %.16e\n",yMin, yMax); 
	//
	fprintf(stdout, "Last Iteration and interior detection, important parameters:\n");
	fprintf(stdout, "\teps = %.16e\n",eps);
	fprintf(stdout, "\tEscape Radius = ER = %.16e\n", ER);
	fprintf(stdout, "\titerMax_LastIteration = %d\n", iterMax_LastIteration);
	fprintf(stdout, "\tPixelWidth = %.16e   \n", PixelWidth); 
	
	//
	PrintInfoAboutPixelSize(iPixelSpacingBits);
	
	printf("==========================================================================================================================\n\n\n\n");
};



// *****************************************************************************
//;;;;;;;;;;;;;;;;;;;;;;  Make Images  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// **************************************************************************************


int MakeImages( const int k){

	const char *Name;
	
	local_setup(k); // for every k 
	
	DrawImage(k, LastIteration,   data); // all the periods 
	Name = GiveName(k, "LastIteration");
	SaveImage(data, Name); 

	/*
	ComputeBoundaries(data,edge);
	Name = GiveName(k, LCM, "LastIteration_LCM");
	SaveImage(edge, Name); 
	
	
	CopyBoundaries(edge, data);
	Name = GiveName(k, "LastIteration");
	SaveImage(data, Name); 
	
	
	DrawImage(k, LSM,   data);
	Name = GiveName(k, LSM, "LSM");
	SaveImage(data, Name); 
	*/
	
	//pseudo_cardioid_center = GiveCenter((int) plane_examples[k][3], data);
	
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
  
  	int k =0; // 
  
	setup();
	
	for (k=0; k<kMax; ++k) {
		MakeImages(k);
	}
		
  	end();
  	
  	
  	return 0;
}

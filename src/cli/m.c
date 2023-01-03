/*

G. Pastor, M. Romera and F. Montoya, “An Approach to the Ordering of One-Dimensional Quadratic Maps,” Chaos, Solitons & Fractals, Vol. 7, No. 4, 1996, pp. 565-584. doi:10.1016/0960-0779(95)00071-2




https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/Parameter_plane#Plane_types






  Adam Majewski
  adammaj1 aaattt o2 dot pl  // o like oxygen not 0 like zero 
  
  
  console program in c programing language 
===============================================================





  
  ==============================================
  
  
  Structure of a program or how to analyze the program 
  
  
  
  Creating graphic:
  * memory array
  * save it to the disk as a pgm file
  * convert pgm file to png usnigng Image Magic convert
  
  * map it to the c plane: for each pixel of plane compute c or lambda using  map_parameter
  
  
  
  

   
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


#define kMax 10 // number of examples, see line 211 plane_examples

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

static unsigned int iHeight = 1000;	//  
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

double  DisplayAspectRatio = 1.0 ; // https://en.wikipedia.org/wiki/Aspect_ratio_(image)



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


/* 
plane : center_x center_y radius
examples containing islands = minibrots
note that plane_center does not equal to hyperbolic component  center

https://mrob.com/pub/muency/mainsequence.html
https://mrob.com/pub/muency/largestislands.html
https://mrob.com/pub/mu-data/largest-islands.txt
*/
double plane_examples[kMax][3] = {
	{-1.0, 		+0.0,		1.5}, // standarad view of whole M set = period 1
	{-1.76733,  	+0.0,		0.05}, // period 3 island
	{+0.2925755,	-0.0149977, 	0.0005}, // period 32 island , distorted
	{-0.15842, 	+1.03335, 	0.03048}, // period 4 
	{+0.358431,	+ 0.643507,	0.017557}, // period 5
	{+0.442990,	+0.373727,	0.011104}, // poeriod 6
	{+0.432259,	+0.227315,	0.007423}, // period 7
	{+0.404879,	+0.146216,	0.005150}, // period 8
	{+0.378631,	+0.098841,	0.003704}, // period 9
	{+0.356854, 	+0.069659,	0.002734} // period 10
};





const int iterMax_LSM = 1001;
const int iterMax_DEM = 1001;
const int iterMax_BET = 1000;

// EscapeRadius for bailout test 
double ER = 2.0;
double ER2;	
double ER_DEM = 100.0;
double ER2_DEM; 













/* colors = shades of gray from 0 to 255 */
unsigned char iColorOfExterior = 250;
unsigned char iColorOfInterior = 200;
unsigned char iColorOfBoundary = 0;
unsigned char iColorOfUnknown = 30;





/* ------------------------------------------ functions -------------------------------------------------------------*/

 double clamp(double x, double lo, double hi) {
  return fmin(fmax(x, lo), hi);
}

//------------------complex numbers -----------------------------------------------------

 double cabs2(complex double z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}



// from screen to world coordinate ; linear mapping
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



complex double fc( const double complex z , const complex double c ){

	return z*z +c;
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
 
  // printf(" find boundaries in S array using  Sobel filter\n");   
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
 
 
  fprintf(stderr, "copy boundaries from S array to D array \n");
  for(iY=1;iY<iyMax-1;++iY)
    for(iX=1;iX<ixMax-1;++iX)
      {i= Give_i(iX,iY); if (S[i]==0) D[i]=0;}
 
 
 
  return 0;
}


// ***************************************************************************************************************************
// **************************  BET = Binary Escape Time *****************************************
// ****************************************************************************************************************************

unsigned char ComputeColorOfBET( complex double c){

	int nMax = iterMax_BET;
  	
  	unsigned char iColor;
	
  	int n;
  	
  	complex double z = critical_point;

  	for (n=0; n < nMax; n++){ //forward iteration
	
    		if (cabs2(z) > ER2) break; // esacping
    	
  		
  		z = fc(z,c); //  for speed only one family here without switch 	
   		 
  	}
  
  	if (n ==nMax)
  		{iColor = 0;} // interior = non escaping set
  		else iColor = 255 - 255.0 * ((double) n)/60; // nMax or lower walues in denominator ; exterior = escaping set
  
  
  	return iColor;


}




// ***************************************************************************************************************************
// ************************** LSM*****************************************
// ****************************************************************************************************************************

unsigned char ComputeColorOfLSM( complex double c){

	int nMax = iterMax_LSM;
  	
  	unsigned char iColor;
	
  	int n;
  	
  	complex double z = critical_point;

  	for (n=0; n < nMax; n++){ //forward iteration
	
    		if (cabs2(z) > ER2) break; // esacping
    	
  		
  		z = fc(z,c); 
   		
  	}
  
  	if (n ==nMax)
  		{iColor = 0;} // interior = non escaping set
  		else iColor = 255 - 255.0 * ((double) n)/60; // nMax or lower walues in denominator ; exterior = escaping set
  
  
  	return iColor;


}


// ***************************************************************************************************************************
// ************************** DEM = exterior DE Method where DE = Distance Estimation  only for z^+c family !!!! ************
// ****************************************************************************************************************************

/*

*c = cexp(c0) + t->center;
  *dc = dc0 * cexp(c0);
}
*/
unsigned char ComputeDolorOfDE(const double complex C )
{
  int i=0; // iteration 
   
   
  double complex Z= 0.0; // initial value for iteration Z0
  double R; // =radius = cabs(Z)
  //double D; 
  double complex dC = 1.0; // derivative
  double de; // = 2 * z * log(cabs(z)) / dc;
  int iMax = iterMax_DEM;
  unsigned char iColor;
   
    
  // iteration = computing the orbit
  for(i=0;i<iMax;i++)
    { 
    	// only for c family 
      dC = 2 * Z * dC + 1.0; 
      Z  = fc(Z, C); // Z*Z+C; // https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/qpolynomials
      
            
     
      if(cabs2(Z) > ER2_DEM) break; // exterior of M set
   
      
    } // for(i=0
   
   
  if (i == iMax) 
  	{iColor = iColorOfInterior;}// interior 
    	else { // exterior and boundary 
    		R = cabs(Z);
    		
    		//cd2 = cd2;
    		
      		de = 2.0 * R * log(R) / cabs(dC)  ; //  2 * cabs(z) * log(cabs(z)) / cabs(dc);
             	
             	             	
             	// choose only ascending part of y = tanh(x)  graph y in [ 0.0,1.0] range
             	//d = clamp( d, 0.0, 1.0);
             	// gray gradient
             	double d = tanh(de/PixelWidth ); //  map to [-3,3] range
             	d = clamp( d, 0.0, 1.0);	
             	// map from floating point in [0,1] range to integer in [0.255] range 
             	iColor = ((int)(d *255.0)) ; 
   
    		}
    
  return iColor; 
}
 
 
 

 
 
/* ==================================================================================================
 ============================= Draw functions ===============================================================
 =====================================================================================================
*/ 
unsigned char ComputeColor(const RepresentationFunctionTypeT RepresentationFunctionType, const complex double c ){

	unsigned char iColor= 0;
	
	
	
	switch(RepresentationFunctionType){
	
		case LSM :{iColor = ComputeColorOfLSM(c); break;}
		
		case DEM : {iColor = ComputeDolorOfDE( c); break; } // 
		
		
	
		default: {}
	
	
	}
	
	return iColor;



}



unsigned char  GiveColor(const RepresentationFunctionTypeT RepresentationFunctionType, const int ix, const int iy){

	complex double c = Give_c(ix,iy);
	unsigned char iColor = ComputeColor(RepresentationFunctionType, c);
	
	return iColor; 

}


// plots  raster point (ix,iy) = computes it's color and save it to the array A
int DrawPoint (const RepresentationFunctionTypeT RepresentationFunctionType,  const int ix, const int iy, unsigned char A[])
{
	
	
	unsigned char iColor = GiveColor( RepresentationFunctionType,  ix, iy);
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
  	
  	

  	fprintf(stderr, "compute image %d RepresentationFunctionType = %d \t  \n", k, RepresentationFunctionType);
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic,1) private(ix,iy) shared(A, ixMax , iyMax)
	// #pragma omp parallel for schedule(dynamic, 1)
  	for (iy = iyMin; iy <= iyMax; ++iy){
    		fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix)
      			{DrawPoint(RepresentationFunctionType,   ix, iy, A);}	//  
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






/*

********************************************* info 

*/




int PrintInfoAboutProgam()
{
	
  
  // 
  printf (" \n");
  

  printf("gcc version: %d.%d.%d\n",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__); // https://stackoverflow.com/questions/20389193/how-do-i-check-my-gcc-c-compiler-version-for-my-eclipse
  // OpenMP version is displayed in the console 
  return 0;
}






	

// uses global var  
// local  set up 
int set_plane(const int k){

	
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
  	
  	
  	zoom = 1.0/plane_radius;
  	
  	  	 
  	
  	return 0;

}


// *****************************************************************************
//;;;;;;;;;;;;;;;;;;;;;;  program setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// **************************************************************************************
// globa; setup = the same for all pixels 
int setup (int k )
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
  
  
  set_plane(k);
    
	
  
  
  
  
  ER2_DEM = ER_DEM*ER_DEM;
  ER2 = ER*ER;
 
	
  
  
   	
  /* create dynamic 1D arrays for colors ( shades of gray ) */
  data = malloc (iSize * sizeof (unsigned char));
  edge = malloc (iSize * sizeof (unsigned char));
   
  //
 
  	
  if (data == NULL || edge == NULL ){
    fprintf (stderr, " Setup error : Could not allocate memory");
    return 1;
  }

  
  
  
  fprintf (stderr," end of setup \n");
	
  return 0;

} // ;;;;;;;;;;;;;;;;;;;;;;;;; end of the setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





const char* GiveName(const int k, const RepresentationFunctionTypeT RepresentationFunctionType, const char* sName)
{

	static char Name[512];
	  	
    	sprintf(Name,"%d_%d_%s", k, RepresentationFunctionType, sName);
    	
    	
    
    	return Name;
}


int PrintInfoAboutImage(){

	fprintf(stdout, "zoom = %.16e\n",zoom);
	fprintf(stdout, "plane radius = %.16e\n",plane_radius);
	fprintf(stdout, "PixelWidth = %.16e   \n", PixelWidth); 
	fprintf(stdout, "xMin = %.16e \t xMax = %.16e\n",xMin, xMax); // 
	fprintf(stdout, "yMin = %.16e \t yMax = %.16e\n",yMin, yMax); 
	
	
	
	printf("==========================================================================================================================\n\n\n\n");
	return 0;

};

int MakeImages( const int k){


	

	const char *Name;

	DrawImage(k, LSM,   data);
	Name = GiveName(k, LSM, "LSM");
	SaveImage(data, Name); 
	
	ComputeBoundaries(data,edge);
	Name = GiveName(k, LCM, "LCM");
	SaveImage(edge, Name); 
	
	//CopyBoundaries(edge, data);
	//shortName = GiveName("LSCM",  ProjectionType);
	//SaveImage(data, shortName); 
	
	
	PrintInfoAboutImage();
	
	

	return 0;

}








int end(){


  fprintf (stderr," allways free memory (deallocate )  to avoid memory leaks \n"); // https://en.wikipedia.org/wiki/C_dynamic_memory_allocation
  free (data);
  free(edge);
  
 PrintInfoAboutProgam();
  
  return 0;

}









// ********************************************************************************************************************
/* -----------------------------------------  main   -------------------------------------------------------------*/
// ********************************************************************************************************************

int main () {
  
  	
  
	
	
	
	
	for (int k=0; k<kMax; ++k) {
			setup(k);
			MakeImages(k);
		}
		
		
	
	
	
	
	
	
  	end();
  	
  	
  	return 0;
}

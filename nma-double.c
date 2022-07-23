#include<stdio.h> 
#include<stdlib.h> 
#include<string.h> 
#include<math.h>

// SSSE3 instrinsic functions. This program requires a cpu with SSSE3 support
#include<tmmintrin.h>

// Has to be changed according to your platform ( lapack header ) 
// Apples optimized version ships with Mac OS X and resides in vecLib framework
#ifdef __APPLE__
#include<vecLib/vecLib.h>
#else
void dsyevd_(char*, char*, int*, double*, int*, double*, double*,
	     int*, int*, int*, int*);
#endif
  
typedef struct {
  double*x;
  double*y;
  double*z;
  int size;
} coordonates;

typedef struct{
  double*eigenvectors;
  double*eigenvalues;
} eigenspace; 

void visualize(coordonates protein,eigenspace space,int eigenmax); 

double* generateCutOffMatrix(double* distancematrix,coordonates protein,
			    double cutoff);

coordonates getCalphaFromXYZ(FILE* xyzfilePointer);

coordonates getCalphaFromPDB(FILE* pdbfilePointer);

double** generatedistancematrix(coordonates protein, double k,double cutoff);

double* generateHessian(double** distancematrix,
			coordonates protein);

// Reads C-alpha Atom cordinates from a xyz file into memory
coordonates getCalphaFromXYZ(FILE* xyzfilePointer) {
  
  char line[90];
  char uninteresting[20];
  int linecount=0;

  coordonates protein; 
  
  while( fgets(line, 90, xyzfilePointer) ) {
    linecount++;
  }

  // allocate memory with 16 byte alignment for eventual vector code

  protein.x = (double*)malloc(sizeof(double)*(linecount+(2-linecount%2)%2));
  protein.y = (double*)malloc(sizeof(double)*(linecount+(2-linecount%2)%2));
  protein.z = (double*)malloc(sizeof(double)*(linecount+(2-linecount%2)%2));
  
  rewind(xyzfilePointer);
  linecount = 0;
  while( fgets(line,90, xyzfilePointer) ) {
    sscanf(line,"%lf %lf %lf",
	   // uninteresting,
	   protein.x+linecount,
	   protein.y+linecount,
	   protein.z+linecount);
    linecount++;
  }
  protein.size=linecount;
  rewind(xyzfilePointer);
  return protein;
};

// Reads C-alpha Atom cordinates from a pdb file into memory
coordonates getCalphaFromPDB(FILE* pdbfilePointer) {
  
  char line[90];
  char atomcheck[20];
  char typecheck[4];
  int number,linecount=0;
  double x,y,z; 
  const char atom[5]="ATOM";
  const char calpha[3]="CA";

  coordonates protein;

  linecount=0;
  
  while( fgets(line,90,pdbfilePointer) ) {
    sscanf(line,"%s", atomcheck);
    if ( strcmp(atom,atomcheck) == 0 ) {
      sscanf(line, "%*s %*i %s",typecheck);
      if ( strcmp(typecheck,calpha) == 0 ) {
	linecount++;	
      }
    }
  }
  
  // allocate memory with 16 byte alignment for eventual vector code
  
  protein.x = (double*)malloc(sizeof(double)*(linecount+(2-linecount%2)%2));
  protein.y = (double*)malloc(sizeof(double)*(linecount+(2-linecount%2)%2));
  protein.z = (double*)malloc(sizeof(double)*(linecount+(2-linecount%2)%2));

  linecount=0;

  rewind(pdbfilePointer);
  while( fgets(line,90,pdbfilePointer) ) {
    sscanf(line,"%s", atomcheck);
    if ( strcmp(atom,atomcheck) == 0 ) {
      sscanf(line, "%*s %*i %s",typecheck);
      if ( strcmp(typecheck,calpha) == 0 ) {
      	sscanf(line+30, "%8lf%8lf%8lf",
	       protein.x + linecount,
	       protein.y + linecount,
	       protein.z + linecount);
	linecount++;
      }
    }
  }
  
  protein.size=linecount; 
  
  rewind(pdbfilePointer);

  return protein;
}

// Generates the distancematrix and calcualtes Hessian elements. Uses SSSE3 
// SIMD intrinsics to speed up calculation as much as possible. A fairly good 
// SSE instruction reference can be otained on intels AVX website: 
// http://www.intel.com/software/avx/ -> intrinsics guide. 
// Attention: Due to the nature of SSE registers mallocs have to be 
// 16 byte aligned. 

double** generatedistancematrix(coordonates protein, double k,double cutoff) {

  int i,j,m; 
  int vectorend = ((protein.size+((2-protein.size%2)%2))/2);
  int scalarrest = (2-protein.size%2)%2;

  double subxscalar, subyscalar, subzscalar; 
  double ** distancematrix = (double**)malloc(sizeof(double*)*7);
  double * cutoffmatrix;

  __m128d subx,suby,subz,subxsquare,subysquare,subzsquare,addition,distance,
    subxsuby,subxsubz,subysubz;

  __m128d vectorx_i,vectory_i,vectorz_i;

  __m128d *distancevector[7] ;
  __m128d *cutoffvector;

  __m128d kvector = _mm_load1_pd(&k);
  __m128d two = _mm_set_pd((double)2.f,(double)2.f);
  __m128d none = _mm_set_pd((double)-1.f,(double)-1.f);
  __m128d ntwo = _mm_mul_pd(two,none);

  __m128d *vx = (__m128d*)protein.x;
  __m128d *vy = (__m128d*)protein.y;
  __m128d *vz = (__m128d*)protein.z;

  // Summation vector for diagonal elements and finalscalar 

  __m128d sum[7];
  double finalscalar;

  // Vector for corrections due to vector offsets
  
  __m128d corrv;
  double * corr = (double*)&corrv;

  // Vector with zeros

  __m128d zero = _mm_setzero_pd();
  
  // masks for vectorized floating point absolute
  // v = (|a1|,|a2|) = _mm_andnot_pd(signmaskvector,a)

  double signmask = (double)0.f * (double)-1.f;
  __m128d signmaskvector = _mm_load1_pd(&signmask);

  // Dinstancematrix size is correlated with vector generation.
 
  for (i=0;i<7;i++) { 
    distancematrix[i] = 
      (double*)malloc((protein.size+scalarrest)
		     *(protein.size+scalarrest)*sizeof(double));
  }
  
  printf("DistanceMatrix: Rank = %i => Hessiangeneration by %i Vectors \n",
	 protein.size,vectorend);
  
  for (i=0;i<protein.size;i++) {

    // load one atoms coordonates into SSE registers

    vectorx_i = _mm_load1_pd(protein.x+i);
    vectory_i = _mm_load1_pd(protein.y+i);
    vectorz_i = _mm_load1_pd(protein.z+i);
    
    for (j=0;j<7;j++){
      
      // Do some pointerarithmetics to get the right line of the matrix

      distancevector[j]=
	(__m128d*)(distancematrix[j]+i*(protein.size+scalarrest));
      distancevector[j][vectorend-1] = zero;
    }

    // Distances and Off Diagonal Hessian Elements

    for (j=0;j<vectorend;j++) {

      // calculation of (|x[i] - x[j]|), (|y[i] - y[j]|), (|z[i] - z[j]|)
      subx = _mm_sub_pd(vectorx_i,vx[j]);
      suby = _mm_sub_pd(vectory_i,vy[j]);
      subz = _mm_sub_pd(vectorz_i,vz[j]);
      
      // calculation of (x[i] - x[j])^2, (y[i] - y[j])^2, (z[i] - z[j])^2
      subxsquare = _mm_mul_pd(subx,subx);
      subysquare = _mm_mul_pd(suby,suby);
      subzsquare = _mm_mul_pd(subz,subz);
      
      // calculation of r=sqrt((x[i]-x[j])^2+(y[i]-y[j])^2+(z[i]-z[j])^2 
      // addition = r^2

      addition = _mm_add_pd(subxsquare,subysquare);
      addition = _mm_add_pd(subzsquare,addition);
      distance = _mm_sqrt_pd(addition); 

      // calculation of crossed (x[i]-x[j])*(y[i]-y[j]), 
      //                        (x[i]-x[j])*(z[i]-z[j]),
      //                        (y[i]-y[j])*(z[i]-z[j])
      subxsuby = _mm_mul_pd(subx,suby);
      subxsubz = _mm_mul_pd(subx,subz);
      subysubz = _mm_mul_pd(suby,subz);

      distancevector[6][j] = distance; 
      
      // Hessian elements (off diagonal)

      distancevector[0][j] = 
	_mm_div_pd( _mm_mul_pd( _mm_mul_pd(subxsquare,ntwo), kvector), 
		    addition);
      distancevector[1][j] =
	_mm_div_pd( _mm_mul_pd( _mm_mul_pd(subysquare,ntwo), kvector),
		    addition);
      distancevector[2][j] = 
	_mm_div_pd( _mm_mul_pd( _mm_mul_pd(subzsquare,ntwo), kvector), 
		    addition);
      
      distancevector[3][j] = 
	_mm_div_pd( _mm_mul_pd( _mm_mul_pd(subxsuby,ntwo), kvector), addition);
      distancevector[4][j] = 
	_mm_div_pd( _mm_mul_pd( _mm_mul_pd(subxsubz,ntwo), kvector), addition);
      distancevector[5][j] =
	_mm_div_pd( _mm_mul_pd( _mm_mul_pd(subysubz,ntwo), kvector), addition);
    }
    
    // Scalar Correction: 
    // The matrix is (2-protein.size%2)%2 larger in order to fit all
    // calculations into SSE registers. Due to the nature of the calculation 
    // this adds unwanted non zero values into this memory area, which is 
    // truncated below. 

    if ( scalarrest != 0 ) {
      
      corrv = zero;
      for(j=0;j<protein.size%2;j++) {
	corr[j] = (double)1.f;
      }
      
      for (j=0;j<7;j++){
	
	distancevector[j]=
	  (__m128d*)(distancematrix[j]+i*(protein.size+scalarrest));
	
	distancevector[j][vectorend-1] = 
	  _mm_mul_pd(corrv,distancevector[j][vectorend-1]);
      }      

    }
    
  }
  
  // Punch some holes!
  // Sets elements zero where where interatom distance > cutoff

  cutoffmatrix = generateCutOffMatrix(distancematrix[6], protein,cutoff);
  cutoffvector = (__m128d*)cutoffmatrix;
  for (i=0;i<(protein.size+scalarrest)*(protein.size+scalarrest)/2;i++) {
    for(j=0;j<6;j++) {
      distancevector[j] = (__m128d*)distancematrix[j];
      distancevector[j][i] = _mm_mul_pd(distancevector[j][i],cutoffvector[i]);
    }
  }
  
  // Compute diagonals
  // Computes diagonal elements. Uses a horitzontal add instruction,
  // replacement of this instruction by a scalar method may add compatability 
  // to SSE 2 capable processors, i.e. > Intel Pentium 4
  
  for (i=0;i<protein.size;i++) {
    
    for (j=0;j<7;j++){
      distancevector[j]=
	(__m128d*)(distancematrix[j]+i*(protein.size+scalarrest));
    }
    
    for (j = 0; j < 6; j++) {      
      distancematrix[j][i*(protein.size+scalarrest)+i] = (double)0.f;
      sum[j] = zero;
    }    
    
    for (j = 0; j < vectorend; j++) {
      for ( m = 0; m < 6; m++) {
	sum[m] = _mm_add_pd(sum[m],distancevector[m][j]);
      }
    }
    
    for (j = 0; j < 6; j++) {
      sum[j] = _mm_hadd_pd(sum[j],zero);
      // Only needed in single precision
      // sum[j] = _mm_hadd_pd(sum[j],zero);
      sum[j] = _mm_mul_pd(sum[j],none);
      _mm_store_sd(&finalscalar,sum[j]);
      distancematrix[j][i*(protein.size+scalarrest)+i] = finalscalar;
    }    
  }
    return distancematrix;
}

// Locates positions where interatom distance is larger then cutoff. Generates
// a binary matrix for.. d[ij] > cutoff => element[ij] = 0 and 1 otherwise. 
double* generateCutOffMatrix(double* distancematrix,coordonates protein,
			    double cutoff) {

  int scalarMatrixSize = (protein.size+((2-protein.size%2)%2))
    *(protein.size+((2-protein.size%2)%2)) ;
  
  int vectorMatrixSize = scalarMatrixSize/2;

  int i;
  
  double *CutOffMatrix = (double*)malloc(scalarMatrixSize*sizeof(double));
  
  __m128d* CutOffMatrixvector = (__m128d*)CutOffMatrix;
  __m128d* distancematrixvector = (__m128d*)distancematrix;
  
  double zeroscalar = (double)0.f;
  double onescalar = (double)1.f;

  __m128d cut = _mm_load1_pd(&cutoff);
  __m128d zero = _mm_load1_pd(&zeroscalar); 
  __m128d one = _mm_load1_pd(&onescalar);

  __m128d comp,compg;

  // Branching in SSE this is way faster then an if condition

  for (i=0;i<vectorMatrixSize;i++) {
    comp = _mm_cmplt_pd(distancematrixvector[i],cut);
    compg = _mm_cmpge_pd(distancematrixvector[i],cut);
    compg = _mm_and_pd(compg,zero);
    comp = _mm_and_pd(comp,one);
    comp = _mm_or_pd(comp,compg);
    CutOffMatrixvector[i] = comp;
  }  
  return CutOffMatrix; 
}

// align hessian elements into a real matrix, in fact this is more or less 
// a preparation to call the lapack routine in order to calculate the 
// eigenvalues. Implemented in scalar code due to the complexity of memory 
// access. Copying values from one location to an other should be fast anyway.

double * generateHessian(double** distancematrix,
			coordonates protein) {

  int i,j;
  int HessianMatrixSize = (9*protein.size*protein.size);
  double * HessianMatrix = (double*)malloc((HessianMatrixSize
					  +(4-HessianMatrixSize%2)%2)
					  *sizeof(double));
  int scalarrest = (2-protein.size%2)%2;
  int vectorend = (protein.size+scalarrest)/2;


  int distancematrixLocation;

  double smallmatrix[9];
  int z = 0;

  for (i=0;i<3*protein.size;i+=3) {
    for (j=0;j<3*protein.size;j+=3) {
      
      distancematrixLocation = i/3*(protein.size+scalarrest)+j/3;
      
	HessianMatrix[(i)*(3*protein.size)+(j)] = // xx
	  distancematrix[0][distancematrixLocation];
	HessianMatrix[(i)*(3*protein.size)+(j)+1] = // xy
	  distancematrix[3][distancematrixLocation];
	HessianMatrix[(i)*(3*protein.size)+(j)+2] = //xz
	  distancematrix[4][distancematrixLocation];
	
	HessianMatrix[(i+1)*(3*protein.size)+(j)] = // xy
	  distancematrix[3][distancematrixLocation];
	HessianMatrix[(i+1)*(3*protein.size)+(j)+1] = //yy
	  distancematrix[1][distancematrixLocation];
	HessianMatrix[(i+1)*(3*protein.size)+(j)+2] = //yz
	  distancematrix[5][distancematrixLocation];
	
	HessianMatrix[(i+2)*(3*protein.size)+(j)] = //xz
	  distancematrix[4][distancematrixLocation];
	HessianMatrix[(i+2)*(3*protein.size)+(j)+1] = //yz
	  distancematrix[5][distancematrixLocation];
	HessianMatrix[(i+2)*(3*protein.size)+(j)+2] = //zz
	  distancematrix[2][distancematrixLocation];
	
	/*		
	HessianMatrix[(i)*(3*protein.size)+(j)] = 0.f; // xx
	HessianMatrix[(i)*(3*protein.size)+(j)+1] = 0.f; // xy
	HessianMatrix[(i)*(3*protein.size)+(j)+2] = 0.f; //xz

	HessianMatrix[(i+1)*(3*protein.size)+(j)] = 0.f;// xy
	HessianMatrix[(i+1)*(3*protein.size)+(j)+1] = 0.f; //yy
	HessianMatrix[(i+1)*(3*protein.size)+(j)+2] = 0.f;//yz

	HessianMatrix[(i+2)*(3*protein.size)+(j)] = 0.f;//xz
	HessianMatrix[(i+2)*(3*protein.size)+(j)+1] = 0.f; //yz
	HessianMatrix[(i+2)*(3*protein.size)+(j)+2] = 0.f;//zz
	*/
    }
    
  }
  
  return HessianMatrix;  
}

// calculates the eigenvalues and eigenvectors of the hessian matrix. Calls
// the SSYEVD subroutine from lapack. 

eigenspace diagonalize(double* HessianMatrix, coordonates protein) {

  int matrixSize = 9*protein.size*protein.size;

  // variables passed to lapack routine

  char wantEigenvectors = 'V';
  char upperTriangle = 'U';
  int order = protein.size*3;
  double* IO = (double*)malloc((matrixSize+(4-matrixSize%2)%2)*sizeof(double));
  double* eigenvalues = (double*)malloc(sizeof(double)*protein.size*3);
  double* work = (double*)malloc(sizeof(double)*3);
  int* iwork = (int*)malloc(sizeof(int)*3);
  int worklength = -1;
  int iworklength = -1;
  int info;

  // end of variables passed to lapack routine

  int i,j;
  eigenspace space; 

  // copies HessianMatrix to an other memory region. 
  
  __m128d * IOV = (__m128d*) IO; 
  __m128d * HessianMatrixV = (__m128d*) HessianMatrix;

  for (i=0; i< (matrixSize+(2-matrixSize%2)%2)/2; i++) {
    IOV[i] = HessianMatrixV[i];
  }
    
  // IO = HessianMatrix; this would be pretty freaky.. 

  dsyevd_(&wantEigenvectors,
	  &upperTriangle,
	  &order,
	  IO,
	  &order,
	  eigenvalues,
	  work,
	  &worklength,
	  iwork,
	  &iworklength,
 	  &info);
  
  // Finds optimal values ( only works for proteins of a certain size )

  worklength = (int)work[0];
  iworklength = iwork[0];

  // Makes this stuff run with some havy proteins 
  // Which are not of a certain size.. ; ) 

  if( worklength < 1 + 6*order + 2*order*order ) {
    worklength = 1 + 6*order + 2*order*order;
  }
  
  if( iworklength < 5*order + 3) {
    iworklength = 5*order + 3;
  }

  free(iwork);
  free(work);
  work = (double*)malloc(sizeof(double)*worklength);
  iwork = (int*)malloc(sizeof(int)*iworklength);

  dsyevd_(&wantEigenvectors,
	  &upperTriangle,
	  &order,
	  IO,
	  &order,
	  eigenvalues,
	  work,
	  &worklength,
	  iwork,
	  &iworklength,
 	  &info);
  
  if (info != 0) {
    printf("Warning! Diagonalisation was not successful!\n");
  }

  free(iwork);
  free(work);

  space.eigenvectors = IO;
  space.eigenvalues = eigenvalues;

  return space;
}

// Outputs visualisation trajectories. 

void visualize(coordonates protein,eigenspace space,int eigenmax) {
  
  int i,j,k;
  double findex;
  
  FILE* visualisationfile;
    
  char filename[13];
    
  int scalarrest = (2-protein.size%2)%2;


  double cosine,factor=0.8f;
  double statex,statey,statez,ampli,prefact;
  
  for(j=6;j<eigenmax;j++) { 
    
    snprintf(filename,13,"%d-mode.xyz",j+1);
    visualisationfile = fopen(filename,"w");    
   
    ampli = space.eigenvalues[j];
    ampli = 1/fabs(ampli);
    ampli = factor*ampli;   
    
    i=0;
    for(findex=0;findex<50;findex++) {
      cosine = cos(M_PI/25*findex);
      prefact = cosine*ampli;
      fprintf(visualisationfile,"%i\nx,y,z\n",protein.size);
      
      for(k=0;k<(protein.size);k++) {
	statex = protein.x[k]+
	  prefact*space.eigenvectors[j*3*protein.size+k*3];
	statey = protein.y[k]+
	  prefact*space.eigenvectors[j*3*protein.size+k*3+1];
	statez = protein.z[k]+
	  prefact*space.eigenvectors[j*3*protein.size+k*3+2];
	fprintf(visualisationfile,"CA %lf %lf %lf\n",statex,statey,statez);
      }    
      i++;
    }
    fclose(visualisationfile);
  }
}
  

// Yes we add this section for the childrens
void usage() {
  printf("\nnma file type k cut nmax \n\n");
  printf("file: pdb or xyz file \n");
  printf("type: distinguishes the file type for the program: \n");
  printf("       type 'pbd' for a pdb file (without ')\n");
  printf("       type 'xyz' for a xyz file (without ')\n");
  printf("   k: spring constant in atom - atom interaction potential \n");
  printf(" cut: Cutoff of atom - atom interaction in Angstr. \n"); 
  printf("nmax: highest mode to be visualized. nmax has to satisfy \n");
  printf("       7 < nmax < 100 \n\n");
}

int main(int argc,char** argv) { 
  
  FILE* file = fopen(argv[1],"r");
  FILE* eigenfile = fopen("eigenvalues", "w");
  FILE* moviefile;
  coordonates protein;
  double** distancematrix, *cutoffmatrix, *hessian;
  int i,j,nmax; 
  eigenspace space;
  
  double k, cutoff;

  // Do some general checking
  
  if ( file == NULL || argc !=6 ) { 
    usage(); 
    return 1; 
  }
  
  if ( 1 != sscanf(argv[3],"%lf",&k) ||
       1 != sscanf(argv[4],"%lf",&cutoff) ||
       1 != sscanf(argv[5],"%i",&nmax)) {
    usage();
    return 1;
  }
  
  if ( nmax < 7 || nmax > 100 ) {
    usage();
    return 1; 
  }

  // Do some more general checking and read contents from file into memory

  if ( 0 == strcmp(argv[2],"pdb") ) {
    protein = getCalphaFromPDB(file);
  } else if ( 0 == strcmp(argv[2],"xyz")) {
    protein = getCalphaFromXYZ(file);
  } else {
    usage();
    return 1;
  }

  fclose(file);

  distancematrix = generatedistancematrix(protein,k,cutoff);
  
  hessian = generateHessian(distancematrix,protein);
  free(distancematrix);
  space = diagonalize(hessian,protein);
  free(hessian);
  
  for (i=0;i<3*protein.size;i++) {
    fprintf(eigenfile,"%lf\n",space.eigenvalues[i]);  
  }
  fclose(eigenfile);
  visualize(protein,space,nmax);
  
  return 0;
}  
  
  

  
    
  

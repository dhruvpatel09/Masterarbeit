/** @file qcd_lanczos.c
 *  @brief Some eigenvalue routines.
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec, Juan Andres Urrea Nino
 *  @date 2011-2025
 *  @copyright GNU Public License.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

#ifndef qcd_USE_MKL 
   #ifdef qcd_USE_LAPACKE
      #include <lapacke.h>
   #endif
#else
   #include <mkl_lapacke.h>
#endif

typedef struct
{
   int degree; //Degree of polynomial.
   double alpha; //a*x + beta
   double beta;
   double gamma; //gamma*T(a*x+b)
} polynomial;

/* private functions */
//------------------------------------------------------------------------------

int setThickness(int m, int nc, int nv)
{
  int k = (int)(nc + ((m - nc)*(3./5. + (double)nv/(100.*m))));
  if (k>nv) 
     k=nv;
  return k;
}

/*Function to generate the alpha coefficient in the Chebyshev transformation*/
double calculateAlpha(double lowerLim, double upperLim)
{
   double alpha;
   alpha = 2.0/(upperLim - lowerLim);
   return alpha;
}

/*Function to generate the beta coefficient in the Chebyshev transformation.*/
double calculateBeta(double alpha, double upperLim)
{
   double beta;
   beta = 1.0 - upperLim*alpha;
   return beta;
}

/*Scalar Chebyshev polynomial*/
double scalar_Chebyshev_N(int N, double x, double val, double prev)
{
   double tmp = val;
   int i;
   for(i = N; i > 0; i--)
   {
      val = 2*x*val - prev;
      prev = tmp;
      tmp = val;
   }
   return prev;
}


/*Function to generate the gamma coefficient in the Chebyshev transformation.*/
double calculateGamma(int N, double beta)
{
   double gamma;
   gamma = 1.0/(scalar_Chebyshev_N(N, beta, beta, 1.0));
   return gamma;
}

/*
* Function to calculate the polynomial parameters.
*/
void setPol(polynomial * Pol, int N, double lowerLim, double upperLim)
{
   double a = 2.0/(upperLim - lowerLim);
   double b = 1.0 - upperLim*a;
   double g = 1.0/(scalar_Chebyshev_N(N, b, b, 1.0));
   Pol->degree = N;  //Set the degree.
   Pol->alpha = a;   //Set the alpha coefficient.
   Pol->beta = b;    //Set beta coefficient.
   Pol->gamma = g;   //Set gamma coefficient.
}

/*
* Minimalistic function to apply the Chebyshev Polynomial.
*/
void loopChebyshev(polynomial Pol, qcd_spinorComponent3d * phi1, qcd_spinorComponent3d * phi2, qcd_geometry * geo, qcd_gaugeField3d * u, double mm)
{
   int i;
   //Extract polynomial information.
   int N = Pol.degree;
   double a = Pol.alpha;
   double b = Pol.beta;
   double g = Pol.gamma;
   //Create the necessary intermediate variables.
   qcd_spinorComponent3d New;
   qcd_spinorComponent3d prev;
   qcd_spinorComponent3d tmp;
   qcd_spinorComponent3d tmp2;
   //Initialize all spinors.
   qcd_initSpinorComponent3d(&New, geo);
   qcd_initSpinorComponent3d(&prev, geo);
   qcd_initSpinorComponent3d(&tmp, geo);
   qcd_initSpinorComponent3d(&tmp2, geo);
   //Assign initial values.
   //----------------------------------------------------------------------------
   //new = (a*L + b)v1.
   qcd_applyCovariantLaplaceSpinorComponent3d(&tmp2, phi1, u, mm);
   qcd_copySpinorComponent3d(&tmp, phi1);
   qcd_scaleSpinorComponent3d(&tmp, b);
   qcd_raxpySpinorComponent3d(&New, a, &tmp2, &tmp);
   //prev = v1;
   qcd_copySpinorComponent3d(&prev, phi1);
   //tmp = new;
   qcd_copySpinorComponent3d(&tmp, &New);
   //Perform loop.
   for (i = N; i > 0; i--)
   {
      //new = 2*(aL + b)new - prev = (2aL + 2b)new - prev.
      qcd_applyCovariantLaplaceSpinorComponent3d(&tmp2, &New, u, mm);
      qcd_scaleSpinorComponent3d(&New, 2*b);
      qcd_raxpySpinorComponent3d(&New, 2*a, &tmp2, &New);
      qcd_subSpinorComponent3d(&New, &New, &prev);
      //prev = tmp.
      qcd_copySpinorComponent3d(&prev, &tmp);
      //tmp = new;
      qcd_copySpinorComponent3d(&tmp, &New);
   }
   //phi2 = g*prev + g. Guarantees positive definite matrix.
   qcd_addSpinorComponent3d(&prev, &prev, phi1);
   qcd_scaleSpinorComponent3d(&prev, g);
   qcd_copySpinorComponent3d(phi2, &prev);
   //Clean up.
   qcd_destroySpinorComponent3d(&New);
   qcd_destroySpinorComponent3d(&prev);
   qcd_destroySpinorComponent3d(&tmp);
   qcd_destroySpinorComponent3d(&tmp2);
}

void calculateRitzVecs(qcd_spinorComponent3d * Q, qcd_spinorComponent3d * V, double * U, int m, int init, int final, int order)
{
   double normQ;
   int i,w;
   if( order==0 ) // assume: smallest eigenvalues come first
   {
      for(i = init; i < final; i++)
      {
         qcd_zeroSpinorComponent3d(&(Q[i]));
         for(w = 0; w < m; w++)
         {
            qcd_raxpySpinorComponent3d(&(Q[i]), U[m*w + i], &(V[w]), &(Q[i]));
         }
         normQ = qcd_normSpinorComponent3d(&(Q[i]));
         qcd_scaleSpinorComponent3d(&(Q[i]), 1.0/normQ);
      }
   }
   else   // assume: smalles eigenvalues come last (=largest of the polynomial)
   {
      for(i = init; i < final; i++)
      {
         qcd_zeroSpinorComponent3d(&(Q[i]));
         for(w = 0; w < m; w++)
         {
            qcd_raxpySpinorComponent3d(&(Q[i]), U[m*w + (m-1-i)], &(V[w]), &(Q[i]));
         }
         normQ = qcd_normSpinorComponent3d(&(Q[i]));
         qcd_scaleSpinorComponent3d(&(Q[i]), 1.0/normQ);
      }
   }
}


/*Function to find maximum of an array.*/
int needOrtho(double * vector, int N, double threshold)
{
   int flag = 0;
   int i;
   for(i = 0; i < N; i++)
   {
      flag += (vector[i] >= threshold);
   }
   return flag;
}

#ifdef qcd_USE_LAPACKE
 
/*Eigendecomposition using MRRR method for tridiagonal matrix.*/
lapack_int useDSTEMR(double * alpha, double * beta, lapack_int n, double * eval, double * evec, double tol)
{
   //Set Lapack parameters.
   double * d = malloc(n*sizeof(double));  //TODO: pre-allocate as part of alpha/beta
   double * e = malloc(n*sizeof(double));
   lapack_int m = n; //Number of eigenpairs found.
   lapack_int * isuppz = malloc(2*m*sizeof(int));
   lapack_logical tryrac = 1;
   if((d==NULL) || (e==NULL) || (isuppz==NULL))
   {
      fprintf(stderr,"out of memory in useDSTEMR\n");
      return(-1);
   }
   //Fill the diagonal and off diagonals.
   memcpy(d,alpha,m*sizeof(double));
   memcpy(e,beta, m*sizeof(double));
   //Query workspaces.
   //LAPACKE_dstemr(LAPACK_ROW_MAJOR, 'V', 'A', n, d, e, 0., 0., 0, n, &m, eval, evec, n, -1, isuppz, &tryrac);
   //Perform calculation.
   lapack_int info = LAPACKE_dstemr(LAPACK_ROW_MAJOR, 'V', 'A', n, d, e, 0., 0., 0, n, &m, eval, evec, n, n, isuppz, &tryrac);
   free(d);
   free(e);
   free(isuppz);
   return info;
}

lapack_int useDSYEVR(double * T, double * eval, double * evec, lapack_int n, double abstol)
{
   lapack_int m = n;
   lapack_int * isuppz = malloc(2*m*sizeof(int));
   if(isuppz==NULL)
   {
      fprintf(stderr,"out of memory in useDSYEVR\n");
      return(-1);
   }
   lapack_int info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'A', 'L', n, T, n, 0., 0., 0, n, abstol, &m, eval, evec, n, isuppz);
   free(isuppz);
   return info;
}

#endif


/* public functions */

int qcd_fsign(double f)
{
   //return (f>0)?1:((f<0)?-1:0);
    return (f>=0)?1:-1;
}//end qcd_fsign

void qcd_printLanczosStatus(qcd_lanczosStatus stat)
{
   if(stat.solverType == 1)
   {
      if(stat.usedLapack)
         printf("qcd_lanczosReortho with lapack was used\n");
      else
         printf("qcd_lanczosReortho without lapack was used\n");
    }else if (stat.solverType == 2)
   {
      if(stat.usedLapack)
         printf("qcd_lanczosCheby with lapack was used\n");
      else
         printf("qcd_lanczosCheby without lapack was used\n");
   } else
   {     
      printf("no eigenvalue problem was solved\n");
      return;
   }
   printf("   restarts %i, # of small systems solved %i, largest %i, max residuum %e\n",stat.restarts, stat.nRed, stat.maxRed, stat.maxRes);
}


void qcd_bubblesort(int *Ind, double *a, int n)
{
   /* take vector a of length n
    * return index-vector Ind such that
    * a[Ind[0]] <= a[Ind[1]] <= a[Ind[2]] <= ... <= a[Ind[n-1]]
    *
    * the rather primitive bubble-sort should be replaced
    * by something efficient in the future
    */

   int i,j;
   int tmp;

   for(i=0; i<n; i++)
      Ind[i]=i;

   for(i=1; i<n; i++)
   {
      j=i;
      while( j>0 && (a[Ind[j]] < a[Ind[j-1]]) )
      {
         //swap j and j-1
         tmp = Ind[j-1];
         Ind[j-1] = Ind[j];
         Ind[j] = tmp;
         j--;
      }
   }
}//end qcd_bubblesort



void qcd_givens(double *c, double *s, double a, double b)
{
   /*  computes c and s such that
    * [ c  s ]T  [ a ]   [ r ]
    * [-s  c ]   [ b ] = [ 0 ]
    *
    * Givens-rotation as in Golub & van Loan p. 216
    */
    double tau;
    if( fabs(b) > fabs(a) )
    {
       tau = -a/b;
       *s  = 1.0/sqrt(1.0+tau*tau);
       *c  = *s * tau;
    }else
    {
       tau = -b/a;
       *c  = 1.0/sqrt(1.0+tau*tau);
       *s  = *c * tau;
    }
    return;
}//end qcd_givens


void qcd_implicitSymmetricQRStepWithShift(double *alpha, double *beta, double *Q, int N, int start, int stop)
{
   /* implicit Symmetric QR Step with Wilkinson shift
    * as in Golub & van Loan algorithm 8.3.2
    *
    * alpha, beta define tridiagonal N x N matrix  T=  [alpha[0] beta[0]  0...            ]
    *                                                  [beta[0]  alpha[1] beta[1]  0 ...  ]
    *                                                  [0        .        .        .      ]
    *                                                  [0               beta[n-1] alpha[n]]
    * Q is a full N x N matrix, Q[i*N+j] ith entry of eigenvector j
    * start, stop defines sub matrix
    * t = T(start:stop, start:stop)
    * and the algorithm overwrites t with Z^T t Z with Z being givens rotations with property
    * Z^T(t-mu) is upper triangular and mu is eigenvalue of t's trailing 2x2 principal submatrix
    * Q is updated: Q <- Q diag(1,Z,1)
    */

   double d, mu, x, z, tt, c, s;
   int k,i;
   double tmp1, tmp2, tmp3;
   double *tmp5, *tmp6;

   if(start>=stop || stop >= N)
      return;

   tmp5 = (double*) malloc(N*sizeof(double));
   tmp6 = (double*) malloc(N*sizeof(double));
   d  = (alpha[stop-1] - alpha[stop]) * 0.5;
   tt = beta[stop-1];
   tt *= tt;
   mu = alpha[stop] - tt/(d + qcd_fsign(d)*sqrt(d*d + tt));
   x = alpha[start] - mu;
   z = beta[start];

   for(k=start; k<stop; k++)
   {
      qcd_givens(&c, &s, x, z);
      //apply Givens matrix to submatrix
      tmp1 =  c*c*alpha[k]    - 2*c*s*beta[k] +s*s*alpha[k+1];
      tmp2 =  c*s*(alpha[k]-alpha[k+1]) +(c*c-s*s)*beta[k];
      tmp3 =  s*s*alpha[k] +    2*c*s*beta[k] +c*c*alpha[k+1];
      if(k>0)
      {
         beta[k-1] = c*beta[k-1] - s*z;
      }
      z    = -s*beta[k+1];                  // the new T_k+2,k
      beta[k+1] *=  c;
      alpha[k]   = tmp1;
      beta[k]    = tmp2;
      alpha[k+1] = tmp3;
      //apply Givens matrix to Q:  Q<- QG
      for(i=0; i<N; i++)
      {
         tmp5[i] = Q[i*N+k]*c - Q[i*N+k+1]*s;
         tmp6[i] = Q[i*N+k]*s + Q[i*N+k+1]*c;
      }
      for(i=0; i<N; i++)
      {
         Q[i*N+k]   = tmp5[i];
         Q[i*N+k+1] = tmp6[i];
      }
      x = beta[k];
   }
   free(tmp5);
   free(tmp6);
   return;
}//end qcd_implicitSymmetricQRStepWithShift


//void qcd_symmetricQR(double *alpha, double *beta, double **Q, int N)
void qcd_symmetricQR(double *al, double *be, int N, double *eval, double *Q, double tol)
{
   /* symmetric QR algorithm as in Golub & van Loan alg. 8.3.3 
    * al, be: diagonal & sub/super diagonal of a symmetric matrix. Both length N
    * Q: N eigenvectors. Q[i*N+j] ith entry of eigenvector j
    */

   int start=0, stop=N-1;
   int i,j;
   int *index;
   double *alpha, *beta;
   
   if (N<1)
      return;
   
   if (N==1)
   {
      eval[0] = al[0];
      Q[0] = 1.0;
      return;
   }
   
   // get workspace 
   alpha = (double*) malloc(N*sizeof(double));
   beta  = (double*) malloc(N*sizeof(double));
   index = (int*) malloc(N*sizeof(int));
   
   if(alpha==NULL || beta==NULL || index==NULL)
   {
      fprintf(stderr,"Error in qcd_symmetricQR: out of memory\n");
      return;
   }
   
   memcpy(alpha, al, N*sizeof(double));
   memcpy(beta, be, N*sizeof(double));
   
   
   
   //initialize eigenvectors
   memset(Q,0,N*N*sizeof(double));
   for(i=0; i<N; i++)
   {
      Q[i*N+i] = 1.0;
   }

   beta[N-1]=0;
   //start loop
   while(1)
   {
      for(i=0; i<N-1; i++)
         if(fabs(beta[i]) <tol*(fabs(alpha[i])+fabs(alpha[i+1])))
            beta[i] = 0;
      //find start and stop
      stop=N-1;
      while(beta[stop-1]==0 && stop>0)
         stop--;

      start=stop;
      while(beta[start-1] != 0 && start>0)
         start--;

      if(stop==0)
         break;
      qcd_implicitSymmetricQRStepWithShift(alpha, beta, Q, N, start, stop);
      //printf("start: %i, stop %i\n",start,stop);
      //printf("%0.16lf %0.16lf %0.16lf\n",beta[0],beta[1],beta[2]);
   }
   
   // sort:

   qcd_bubblesort(index,alpha,N);
   //printf("sorting %i eigenvalues and eigenvectors\n",N);
   //reorder vectors
   for(j=0; j<N; j++) // components
   {
       //use eval as tmp space, reorder j component of every eigenvector
       for (i=0; i<N; i++)
           eval[i] = Q[j*N+index[i]];
       for (i=0; i<N; i++)
           Q[j*N+i] = eval[i];
   }
   //reorder eigenvalues
   for(i=0; i<N; i++)
   {
      //printf("eval_%i = %f at index %i\n",i, alpha[index[i]],index[i]); fflush(stdout);
      eval[i] = alpha[index[i]];
   }
   
   free(alpha);
   free(beta);
   free(index);
}//end qcd_symmetricQR

void qcd_defaultLanczosParams(qcd_lanczosParams *par)
{
   par->res=1e-10;
   par->redTol=1e-12;
   par->nLanczosFactor=6;
   par->deg=4;
   par->upper=10.0;
   par->tolReortho=1e-10;
   par->maxRestarts=30;
   par->redFreqFactor=1;
}


int qcd_lanczosCheby(double *eval, qcd_spinorComponent3d *evec, qcd_lanczosStatus *stat, qcd_gaugeField3d *u, int nvec, qcd_lanczosParams *par)
{
   polynomial pol;
   qcd_geometry *geo;
   double n0;
   qcd_spinorComponent3d res;
   qcd_spinorComponent3d *V; // Lanczos Vectors
   //qcd_spinorComponent3d *Q; // Ritz Vectors
   int reorthoCount = 0;
   double normR, resErr;            //Norms for quality check
   int qcdOut = 0;                         //Error status
   qcd_complex_16 overlap;                 //Overlap
   qcd_complex_16 tmpC;
   double *alpha;     //Alpha coefficients
   double *beta;      //Beta coefficients
   double *O1,*O2,*O3;//Orthogonality estimators
   double *T;         //Reduced Lanczos matrix
   double *evecT;     //Eigenvectors
   double *evalT;     //Eigenvalues
   double la;
   double phase;
   int i,j;           //loop variables
   int m=nvec*par->nLanczosFactor;
   int k;             //Thickness (lanczos)
   int counter;       //counter for lanczos iterations
   double maxRes=0.0; //return value
   int restarts=0;
   int nRed=0;
   int maxRed=0;
   int myid = u->geo->myid;
   int root3d = u->geo->root3d[u->geo->Pos[0]*u->geo->lL[0]];
   
  
   //test initialization of u, evec
   
   geo = u->geo;
   if (stat!=NULL)
   {
      stat->restarts=0;
      stat->nRed=0;
      stat->maxRed=0;
      stat->maxRes=100;
      stat->solverType=2;
#ifdef qcd_USE_LAPACKE
         stat->usedLapack=1;
#else
         stat->usedLapack=0;
#endif
   }
   setPol(&pol, par->deg, par->upper, 15.0);
   
   //if(geo->myid==0){ printf("qcd_lanczosCheby: polynomial set\n"); fflush(stdout);}

   //----------------------------------------------------------------------------
   //allocate memory
   // TODO: test whether allignment speeds things up
   qcdOut += qcd_initSpinorComponent3d(&res, geo);  //residual vector
   alpha = malloc((m + 1)*sizeof(double));  //Alpha coefficients.
   qcdOut += (alpha==NULL);
   beta = malloc((m + 1)*sizeof(double));   //Beta coefficients.
   qcdOut += (beta==NULL);
   O1 = malloc((m + 1)*sizeof(double));     //Orthogonality estimators.
   qcdOut += (O1==NULL);
   O2 = malloc((m + 1)*sizeof(double));
   qcdOut += (O2==NULL);
   O3 = malloc((m + 1)*sizeof(double));
   qcdOut += (O3==NULL);
   T = malloc(m*m*sizeof(double));          //Reduced Lanczos matrix.
   qcdOut += (T==NULL);
   evecT = malloc(m*m*sizeof(double));      //Eigenvectors.
   qcdOut += (evecT==NULL);
   evalT = malloc(m*sizeof(double));        //Eigenvalues
   qcdOut += (eval==NULL);
   V = (qcd_spinorComponent3d*) malloc((m+2)*sizeof(qcd_spinorComponent3d));
   qcdOut += (V==NULL);
   //Q = (qcd_spinorComponent3d*) malloc(m*sizeof(qcd_spinorComponent3d));
   //qcdOut += (Q==NULL);
   if(qcdOut > 0)
   {
      fprintf(stderr,"Process %i: Out of memory in qcd_lanczosCheby\n",geo->myid);
      return 1;
   }
   for(i = 0; i < m + 2; i++)
   {
      qcdOut += qcd_initSpinorComponent3d(&(V[i]), geo);
   }
   if(qcdOut > 0)
   {
      fprintf(stderr,"Process %i: Out of memory in qcd_lanczosCheby\n",geo->myid);
      return 1;
   }
  
   counter = 1;
   //Initialize error estimators to 0.
   memset(O1, 0, (m+1)*sizeof(double));
   memset(O2, 0, (m+1)*sizeof(double));
   memset(O3, 0, (m+1)*sizeof(double));
   
   restarts=0;
   nRed=0;

   
   // main lanczos iteration
   while(counter <= par->maxRestarts)
   {
      if(counter > 1) //restart
      {
         restarts++;
         //Set thickness.
         //k = setThickness(m, nc, nvec);
         k = nvec; //fixed restart width
         //if(geo->myid==0){ printf("qcd_lanczosCheby: restart with thickness %i\n",k); fflush(stdout);}
         
         //Set initial k Lanczos vectors with Ritz vectors.
         for(i = 0; i < k; i++)
         {
            qcd_copySpinorComponent3d(&(V[i+1]), &(evec[i]));
         }
         //Set k+1-th Lanczos vector as the previous residual vector.
         qcd_copySpinorComponent3d(&(V[k+1]), &(V[m+1]));
      }
      else //first run
      {
         k = 0;  //thickness.
         //Create initial constant vector with norm 1.
         qcd_setSpinorComponent3d(&(V[1]), (qcd_complex_16) {1.0,0.0});
         n0 = qcd_normSpinorComponent3d(&(V[1]));
         qcd_scaleSpinorComponent3d(&(V[1]), 1.0/n0);
         //if(geo->myid==0){ printf("qcd_lanczosCheby: first vector initialized\n"); fflush(stdout);}
      }
      //r = Ap*V(:,k+1);
      loopChebyshev(pol, &(V[k+1]), &(V[0]), geo, u, 0.0);
      //alpha(k+1) = V(:,k+1)'*r;
      alpha[k+1] = qcd_mulAdjointSpinorComponent3dSpinorComponent3dR(&(V[k+1]), &(V[0]));
      //r = r - alpha(k+1).*V(:,k+1);
      qcd_raxpySpinorComponent3d(&(V[0]), -alpha[k+1], &(V[k+1]), &(V[0]));
      //Orthogonalize restarting vector against all previous ones.
      for(i = 1; i <= k; i++)
      {
         beta[i] = qcd_mulAdjointSpinorComponent3dSpinorComponent3dR(&(V[i]), &(V[0]));
         qcd_raxpySpinorComponent3d(&(V[0]), -beta[i], &(V[i]), &(V[0]));
      }
      beta[k+1] = qcd_normSpinorComponent3d(&(V[0]));
      //V(k+2) = r/beta(k+1).
      qcd_copySpinorComponent3d(&(V[k+2]), &(V[0]));
      qcd_scaleSpinorComponent3d(&(V[k+2]), 1/beta[k+1]);
      //Loop to calculate vectors k+3,...,m+1
      for(i = k + 2; i < m + 1; i++)
      {
         //if(geo->myid==0){ printf("qcd_lanczosCheby: i=%i\n",i); fflush(stdout);}
         //r = Ap*V(:,i);
         loopChebyshev(pol, &(V[i]), &(V[0]), geo, u, 0.0);
         //alpha(i) = V(:,i)'*r;
         alpha[i] = qcd_mulAdjointSpinorComponent3dSpinorComponent3dR(&(V[i]),&(V[0]));
         //r = r - alpha(i).*V(:,i) - beta(i-1).*V(:,i-1)
         //1. r = -alpha(i)*V(i) + r.
         qcd_raxpySpinorComponent3d(&(V[0]), -alpha[i], &(V[i]), &(V[0]));
         //2. r = -beta(i-1)*V(i-1) + r.
         qcd_raxpySpinorComponent3d(&(V[0]), -beta[i-1], &(V[i-1]), &(V[0]));
         //beta(i) = norm(r);
         beta[i] = qcd_normSpinorComponent3d(&(V[0]));
         /*
         * Reorthogonalization schemes:
         * Initial Run: Partial.
         * Restarts: Full.
         */
         //For restarts.
         if (counter > 1)
         {
            //Full reorthogonalization.
            for(j = 1; j < i + 1; j++)
            {
               overlap = qcd_mulAdjointSpinorComponent3dSpinorComponent3d(&(V[j]),&(V[0]));
               overlap = (qcd_complex_16) {-overlap.re, -overlap.im};
               qcd_axpySpinorComponent3d(&(V[0]), overlap, &(V[j]), &(V[0]));
            }
         }
         else //initial run
         {
            //Periodic Reorthogonalization.
            //Calculate orthogonality estimates.
            O2[i] = 1.0;
            O2[i-1] = 2.22e-16*beta[2]*0.5/beta[i-1];
            O3[1] = (alpha[1] - alpha[i])*O2[1] + beta[1]*O2[2] - beta[i-1]*O1[1];
            O3[1] = O3[i]/beta[i];
            for(j = 2; j <= i-2; j++)
            {
               O3[j] = (alpha[j] - alpha[i])*O2[j] + beta[j]*O2[j+1] - beta[i-1]*O1[j] + beta[j-1]*O2[j-1] + 2.22e-16*0.5*(beta[j+1] + beta[i]);
               O3[j] /= beta[i];
            }
            
            //If needed, perform reorthogonalization.
            if(needOrtho(O3, i-1, par->tolReortho) > 0)
            {
               reorthoCount += 1;
               //Orthogonalize vectors i and i+1 against all previous ones.
               for(j = 1; j < i; j++)
               {
                  //Vector i.
                  overlap = qcd_mulAdjointSpinorComponent3dSpinorComponent3d(&(V[j]),&(V[i]));
                  overlap = (qcd_complex_16) {-overlap.re, -overlap.im};
                  qcd_axpySpinorComponent3d(&(V[i]), overlap, &(V[j]), &(V[i]));
                  O2[j] = 0.0;
                  //Vector i+1
                  overlap = qcd_mulAdjointSpinorComponent3dSpinorComponent3d(&(V[j]), &(V[0]));
                  overlap = (qcd_complex_16) {-overlap.re, -overlap.im};
                  qcd_axpySpinorComponent3d(&(V[0]), overlap, &(V[j]), &(V[0]));
                  O3[j] = 0.0;
               }
               //Orthogonalize i+1 with i.
               overlap = qcd_mulAdjointSpinorComponent3dSpinorComponent3d(&(V[i]), &(V[0]));
               overlap = (qcd_complex_16) {-overlap.re, -overlap.im};
               qcd_axpySpinorComponent3d(&(V[0]), overlap, &(V[i]), &(V[0]));
            }
            //Update the orthogonality estimates for next iteration.
            memcpy(O1,O2,(i-1)*sizeof(double)); //Copy O2 to O1.
            memcpy(O2,O3,(i-1)*sizeof(double)); //Copy O3 to O2.
         }
         //At this point orthogonalization has been performed.
         //beta(i) = norm(r);
         beta[i] = qcd_normSpinorComponent3d(&(V[0]));
         //V(:,i+1) = r./beta(i);
         qcd_copySpinorComponent3d(&(V[i+1]), &(V[0]));
         qcd_scaleSpinorComponent3d(&(V[i+1]), 1.0/beta[i]);
      }//m Lanczos iterations done
      
      //Reduced eigenvalue problem.
      //For restarts.
      if(counter > 1)
      {
         //Populate Tm with zeros.
         memset(T, 0, m*m*sizeof(double));
         /*First k alpha and beta coefficients.*/
         for(i = 1; i <= k; i++)
         {
            alpha[i] = evalT[m - i];
            T[m*(i-1) + i - 1] = alpha[i];
            T[m*k + i - 1] = beta[i];
            T[m*(i-1) + k] = beta[i];
         }
         /*Input last m-k alpha and beta coefficients.*/
         for(i = k; i < m - 1; i++)
         {
            T[m*i + i + 1] = beta[i+1];
            T[m*(i+1) + i] = beta[i+1];
            T[m*i + i] = alpha[i+1];
         }
         T[m*(m - 1) + m - 1] = alpha[m];
#ifdef qcd_USE_LAPACKE
         /*Use MR^3.*/
         if(useDSYEVR(T, evalT, evecT, m, par->redTol) !=0)
         {
            fprintf(stderr,"Error in useDSYEVR, qcd_lanczos.c\n");
            return 3;
         }
         maxRed=m;
         nRed++;
#else
         fprintf(stderr, "No restarts possible (yet) without lapack!\n");
         stat->usedLapack=0;
         goto CLEANUP;
         //TODO: implement jacobi
#endif
         //if(geo->myid==0){ printf("qcd_lanczosCheby: reduced problem solved (DSYEVR)\n"); fflush(stdout);}
      }
      else //initial run
      {
         //Here we use dstemr (Robust Representations)
#ifdef qcd_USE_LAPACKE
         if( useDSTEMR(&(alpha[1]), &(beta[1]), m, evalT, evecT,par->redTol) != 0)
         {
            fprintf(stderr,"Error in useDSTEMR, qcd_lanczos.c\n");
            return 4 ;
         }
         maxRed=m;
         nRed++;
#else
         qcd_symmetricQR(&(alpha[1]), &(beta[1]), m, evalT, evecT,par->redTol);
         maxRed=m;
         nRed++;
#endif
      }

      calculateRitzVecs(evec, &(V[1]), evecT, m, 0, nvec, 1);
      //if(geo->myid==0){ printf("qcd_lanczosCheby: eigenvectors computed\n"); fflush(stdout);}
      
      maxRes=0;
      for(i = nvec-1; i >= 0; i--)
      {
         //Calculate eigenvalue of Laplacian.
         //(evalT are P(lambda)
         qcd_applyCovariantLaplaceSpinorComponent3d(&res, &(evec[i]), u, 0.0);
         la = qcd_mulAdjointSpinorComponent3dSpinorComponent3dR(&(evec[i]), &res);
         eval[i] = la;
         qcd_raxpySpinorComponent3d(&res, -la, &(evec[i]), &res);
         normR = qcd_normSpinorComponent3d(&res);
         resErr = normR/fabs(la);
         maxRes = fmax(maxRes,resErr);
         if(maxRes>par->res && counter<par->maxRestarts)
         {
            break; // next iteration
         }
      }
      counter += 1;
      if (maxRes<par->res)
          break; //while loop
   }//while-loop (thick restarts)

   //--------------------------------------------------------------------------
   //Fix phase convention --> First entry of each eigenvector is real and positive.
   //Iterate over all relevant eigenvectors.
   for(i = 0; i < nvec; i++)
   {
      tmpC = evec[i].D[0][0];     //Grab the first entry.
      phase = qcd_CARG(tmpC);    //Extract its phase.
      MPI_Bcast(&phase, 1, MPI_DOUBLE, 0, geo->comm3d);
      //Multiply the entire eigenvector by the opposite phase.
      qcd_mulSpinorComponent3dC(&(evec[i]), (qcd_complex_16) {cos(phase), -sin(phase)});
   }
  
  
   //cleanup
#ifndef qcd_USE_LAPACKE
   CLEANUP:
#endif
   if (stat != NULL)
   {
      stat->restarts=restarts;
      stat->nRed=nRed;
      stat->maxRed=maxRed;
      stat->maxRes=maxRes;
   }
   
   qcd_destroySpinorComponent3d(&res);
   free(alpha);
   free(beta);
   free(O1);
   free(O2);
   free(O3);
   free(T);
   free(evecT);
   free(evalT);
   for(i = 0; i < m + 2; i++)
   {
      qcd_destroySpinorComponent3d(&(V[i]));
   }
   free(V);
   if( maxRes>par->res )
   {
      if(myid==root3d)
         fprintf(stderr,"Warning! Lanczos did not converge! Largest residuum %e, goal %e\n",maxRes,par->res);
      return 2;
   }
   if( eval[nvec-1] >= par->upper )
   {
      if(myid==root3d)
         fprintf(stderr,"Warning! largest eigenvalue not smaller than par.upper, %e >= %e\n",eval[nvec-1],par->upper);
      return 5;
   }
   return 0;
}//end qcd_lanczosCheby


int qcd_lanczosReortho(double *eval, qcd_spinorComponent3d *evec, qcd_lanczosStatus *stat, qcd_gaugeField3d *u, int nvec, qcd_lanczosParams *par)
{
   /* calculates the lowest n eigenvalues eval and eigenvectors evec of
    * the 3d covariant laplacian with mass-shift 0.0
    * u is the gauge-field of the slice
    *
    * The first par.nLanczosFactor*n Lanczos iterations with complete re-orthogonalization
    * are carried out (memory!)
    *
    * the quality of an eigenvalue/eigenvector pair is
    *  q = (A*evec - eval*evec)
    * If the biggest q out of the first n is larger than the precision goal, an error is returned
    */

   qcd_spinorComponent3d *q;
   qcd_spinorComponent3d tmp, tmp2;
   int i,j,k,l;
   int nLanczos = nvec*par->nLanczosFactor;
   double *beta;
   double *alpha;
   double *smallEvec;
   double *smallEval;
   int myid = evec[0].geo->myid;
   qcd_complex_16 overlap, tmpC;
   double nrm, phase;
   double quality=0;
   qcd_geometry *geo = u->geo;
   int root3d = u->geo->root3d[geo->Pos[0]*geo->lL[0]];
  int nRed=0, maxRed=0;
   
   if (stat != NULL)
   {
      stat->restarts=0;
      stat->nRed=0;
      stat->maxRed=0;
      stat->maxRes=100;
      stat->solverType=1;
#ifdef qcd_USE_LAPACKE
         stat->usedLapack=1;
#else
         stat->usedLapack=0;
#endif
   }
   beta      = (double*) malloc((nLanczos+1)*sizeof(double));
   alpha     = (double*) malloc((nLanczos+1)*sizeof(double));
   smallEvec = (double*) malloc(nLanczos*nLanczos*sizeof(double));
   smallEval = (double*) malloc(nLanczos*sizeof(double));
   j=0;
   if(smallEvec==NULL || smallEval==NULL || alpha==NULL || beta==NULL)
   {
      if(myid==root3d) fprintf(stderr,"Out of memory in lanczos\n");
      j++;
   }

   q = (qcd_spinorComponent3d*) malloc((nLanczos+2)*sizeof(qcd_spinorComponent3d));
   if( q==NULL )
   {
      if(myid==root3d) fprintf(stderr,"Out of memory in lanczos\n");
      j++;
   }
   for(i=0; i<nLanczos+2; i++)
      j+= qcd_initSpinorComponent3d(&(q[i]),geo);

   j += qcd_initSpinorComponent3d(&tmp, geo);
   j += qcd_initSpinorComponent3d(&tmp2, geo);
   MPI_Bcast(&j, 1, MPI_INT, 0, geo->comm3d);
   if(j>0)
   {
      if(myid==root3d) fprintf(stderr,"Out of memory in lanczos\n");      
      return 1;
   }

   qcd_zeroSpinorComponent3d(&(q[0]));
   // random starting vector
   qcd_rand(&(q[1].D[0][0].re), geo->lV3*6);
   q[1].BMup2date = 0;
   q[1].BPup2date = 0;
   qcd_scaleSpinorComponent3d(&q[1], 1.0/qcd_normSpinorComponent3d(&q[1]));
   alpha[0] = 0.0;
   beta[0] = 0.0;
   k=1;

   while(k<nLanczos) //lanczos iteration with complete reorthogonalization
   {
      qcd_applyCovariantLaplaceSpinorComponent3d(&tmp, &(q[k]), u, 0.0);
      qcd_addSpinorComponent3d(&(q[0]), &(q[0]), &tmp);
      alpha[k] = qcd_mulAdjointSpinorComponent3dSpinorComponent3dR(&(q[k]), &tmp);
      qcd_raxpySpinorComponent3d(&(q[0]), -alpha[k], &(q[k]), &q[0]);
      //reorthogonalize q[0] against q[1]...q[k]
      for(j=1; j<=k; j++)
      {
         overlap = qcd_mulAdjointSpinorComponent3dSpinorComponent3d(&(q[j]),&(q[0]));
         overlap = (qcd_complex_16) {-overlap.re, -overlap.im};
         qcd_axpySpinorComponent3d(&(q[0]),overlap, &(q[j]), &(q[0]));
      }
      beta[k] = qcd_normSpinorComponent3d(&(q[0]));
      qcd_copySpinorComponent3d(&(q[k+1]), &(q[0]));
      qcd_scaleSpinorComponent3d(&(q[k+1]),1.0/beta[k]);
      qcd_copySpinorComponent3d(&(q[0]), &(q[k]));
      qcd_scaleSpinorComponent3d(&(q[0]),-beta[k]);
      //k++;

      if((k % (par->redFreqFactor*nvec)) == 0 || k==nLanczos-1) // construct vectors and see whether they're any good
      {
         nRed++;
         maxRed=k;
         //diagonalize real, symmetric tri-diagonal matrix diag(beta,alpha,beta)
#ifdef qcd_USE_LAPACKE
         useDSTEMR(&(alpha[1]), &(beta[1]), k, smallEval,  smallEvec, par->redTol);
#else
         qcd_symmetricQR(&(alpha[1]), &(beta[1]), k, smallEval, smallEvec, par->redTol);
#endif
         calculateRitzVecs(evec , &(q[1]), smallEvec, k, 0, nvec, 0);
         quality=0;
             
         for(l=nvec-1; l>=0; l--)
         {
            //test how good the eigenvector/eigenvalue pair is
            qcd_applyCovariantLaplaceSpinorComponent3d(&tmp, &(evec[l]), u, 0.0);
            //nrm = qcd_mulAdjointSpinorComponent3dSpinorComponent3dR(&(evec[l]),&tmp);
            qcd_copySpinorComponent3d(&tmp2,&(evec[l]));
            qcd_scaleSpinorComponent3d(&tmp2,smallEval[l]);
            qcd_subSpinorComponent3d(&tmp, &tmp, &tmp2);
            nrm = qcd_normSpinorComponent3d(&tmp)/fabs(smallEval[l]);
            if( nrm>quality ) quality=nrm;
            if (quality> par->res) break;
         }
         if(quality < par->res) k=nLanczos+1; // end while loop
      }
      k++;
   }

   // normalize eigenvectors; first entry positive,real
   for(i = 0; i < nvec; i++)
   {
      //nrm = qcd_normSpinorComponent3d(&evec[i]);
      //qcd_scaleSpinorComponent3d(&evec[i],1.0/nrm);
      tmpC = evec[i].D[0][0];    //Grab the first entry.
      phase = qcd_CARG(tmpC);    //Extract its phase.
      MPI_Bcast(&phase, 1, MPI_DOUBLE, 0, geo->comm3d);
      //Multiply the entire eigenvector by the opposite phase.
      qcd_mulSpinorComponent3dC(&(evec[i]), (qcd_complex_16) {cos(phase), -sin(phase)});
   }

   memcpy(eval,smallEval,nvec*sizeof(double));
   
   if (stat != NULL)
   {
      stat->restarts=0;
      stat->nRed=nRed;
      stat->maxRed=maxRed;
      stat->maxRes=quality;
   }

   //clean up
   for(i=0; i<nLanczos; i++)
      qcd_destroySpinorComponent3d(&(q[i]));
   qcd_destroySpinorComponent3d(&tmp);
   qcd_destroySpinorComponent3d(&tmp2);
   free(alpha);
   free(beta);
   free(smallEvec);
   free(smallEval);
   free(q);
   if( quality<par->res )
   {
      if(myid==root3d)
         fprintf(stderr,"Warning! Lanczos did not converge! Largest residuum %e, goal %e\n",quality,par->res);
      return 2;
   }
   return 0;
}//end qcd_lanczosReortho

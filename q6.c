#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <lapacke.h>

void read_input_data(double *Lx, double *Ly, int *nx, int *ny, double *tf, double *td, double *v, double *u0x, double *u0y);
void read_coeff_data(long A, double *K);

struct band_mat{
	long ncol;        /* Number of columns in band matrix            */
	long nbrows;      /* Number of rows (bands in original matrix)   */
	long nbands_up;   /* Number of bands above diagonal           */
	long nbands_low;  /* Number of bands below diagonal           */
	double *array;    /* Storage for the matrix in banded format  */
	/* Internal temporary storage for solving inverse problem */
	long nbrows_inv;  /* Number of rows of inverse matrix   */
	double *array_inv;/* Store the inverse if this is generated */
	int *ipiv;        /* Additional inverse information         */
};

typedef struct band_mat band_mat;

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters.  */ 
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
	bmat->nbrows 		= nbands_lower + nbands_upper + 1;
	bmat->ncol   		= n_columns;
	bmat->nbands_up 	= nbands_upper;
	bmat->nbands_low 	= nbands_lower;
	bmat->array      	= (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
	bmat->nbrows_inv 	= bmat->nbands_up*2 + bmat->nbands_low + 1;
	bmat->array_inv  	= (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
	bmat->ipiv       	= (int *) malloc(sizeof(int)*bmat->ncol);
	
	if (bmat->array==NULL||bmat->array_inv==NULL) {
		return 0;
	}  
	/* Initialise array to zero */
	long i;
	for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
		bmat->array[i] = 0.0;
	}
  return 1;
}

/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double *getp(band_mat *bmat, long row, long column) {
	int bandno = bmat->nbands_up + row - column;
	if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
		printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
		exit(1);
	}
	return &bmat->array[bmat->nbrows*column + bandno];
}

/* Retrun the value of a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double getv(band_mat *bmat, long row, long column) {
	return *getp(bmat,row,column);
}

void setv(band_mat *bmat, long row, long column, double val) {
	*getp(bmat,row,column) = val;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays                                          */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
	/* Copy bmat array into the temporary store */
	int i,bandno;
	for(i=0;i<bmat->ncol;i++) { 
		for (bandno=0;bandno<bmat->nbrows;bandno++) {
			bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] 
			= bmat->array[bmat->nbrows*i+bandno];
		}
    x[i] = b[i];
	}

	long nrhs = 1;
	long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
	int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, 
		bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
	return info;
}

int printmat(band_mat *bmat) {
	long i,j;
	for(i=0; i<bmat->ncol;i++) {
		for(j=0; j<bmat->nbrows; j++) {
			printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
		}
	}
	return 0;
}

/*Check that a grid point has valid coordinates */
int is_valid(long j, long p, long J, long P) {
	return (j>=0)&&(j<J)&&(p>=0)&&(p<P);
}

/* Return the 1D element index corresponding to a particular grid point.
   We can rewrite this function without changing the rest of the code if
   we want to change the grid numbering scheme!
   Output: long integer with the index of the point
   Input:
   long j:  The X grid point index
   long k:  The Y grid point index
   long P:  The number of Y points.
*/
long indx( long j, long p, long P, long J) {

	// include the boundary conditions
	if(j==-1){j=J-1;}
	if(j==J){j=0;}
	if(p==-1){p=P-1;}
	if(p==P){p=0;}

	return j*P + p;

}

/* Return the 2D point corresponding to a particular 1D grid index */
void gridp(long indx, long P, long *j, long *p) {
	*j = indx%P;
	*p = indx - (*j)*P;
}


long folding(long i, long N) {
	long x;

	if (i < 0.5*N){
		x = 2*i;
	}
	else {
		x = 2*(N-i)-1;
	}
	return x;
}

long unfolding(long i, long N) {
	long x;

	if (i%2 == 0){
		x = i/2.0;
	}
	else {
		x = N-(i+1)/2.0;
	}
	return x;
}

void arrayfolded(double *folded, double *unfoldedarray, long A, int nx, int ny) {
	int unknown_indx;

	for(int j=0; j<nx; j++) {
		for(int p=0; p<ny; p++) {
			unknown_indx = indx(j,p,ny,nx);

			folded[unknown_indx] = unfoldedarray[folding(unknown_indx, A)];

		}
	}
}

void arrayunfolded(double *unfolded, double *foldedarray, long A, int nx, int ny) {
	int unknown_indx;

	for(int j=0; j<nx; j++) {
		for(int p=0; p<ny; p++) {
			unknown_indx = indx(j,p,ny,nx);

			unfolded[unknown_indx] = foldedarray[unfolding(unknown_indx, A)];

		}
	}
}

void freeforall(band_mat *bmat, double* w, double * w_next, double* phi, double * foldedw, double* unfoldedphi, double *ux, double *uy, double *K) {
	free(bmat->array); free(bmat->array_inv); free(bmat->ipiv); free(w); free(w_next); free(phi); free(ux); free(uy); free(K); free(foldedw); free(unfoldedphi);
}

//====================================================================================================

int main() {
	FILE * output; //create A reference to the output file
	output = fopen("output.txt", "w"); //create an output file if not already existing

	//Define pointers for all of the variables, to pass them to the input function:
	double 	Lx;			double 	*pLx		= &Lx; 	// right x boundary
	double 	Ly;			double 	*pLy		= &Ly; 	// right y boundary
	int 	nx;			int 	*pnx 		= &nx; 	// no. of x grid points
	int 	ny;			int 	*pny 		= &ny; 	// no. of y grid points
	double 	tf;			double 	*ptf		= &tf; 	// simulation time
	double 	td;			double 	*ptd		= &td; 	// diagnostic timestep
	double 	v;			double 	*pv			= &v; 	// viscosity
	double 	u0x;		double 	*pu0x		= &u0x; // av. velocity (x)
	double 	u0y;		double 	*pu0y		= &u0y; // av velocity (y)

	//Read the required parameters in from input.txt
	read_input_data(pLx,pLy,pnx,pny,ptf,ptd,pv,pu0x,pu0y);

	long A = nx*ny; /* Total size of problem is number of grid points on 2D plane */

	//Allocate sensible memory
	double *K = malloc(A*sizeof(double));
	//Read the coefficients from coefficients.txt
	read_coeff_data(A,K);
	//Compute the grid spacing that will be used in the finite difference method for solving the PDEs
	double dx = Lx/(nx);
	double dy = Ly/(ny);
	double dt = 0.01*dx*dy;

	int nsteps_per_output = td/dt + 1; 
  	dt =  td/nsteps_per_output;
	
	// Prepare grid storage
	double *phi; 	phi = malloc(A*sizeof(double)); 	// flow field at current vorticity
	double *ux; 	ux 	= malloc(A*sizeof(double));		// current velocity (x)
	double *uy; 	uy 	= malloc(A*sizeof(double)); 	// current velocity (y)
	double *w;		w 	= malloc(A*sizeof(double)); 	// current vorticity and initialisation
	double *w_next 		= malloc(A*sizeof(double)); 	// current vorticity
	double *foldedw 	= malloc(A*sizeof(double)); 	// folded arrray for omega
	double *unfoldedphi	= malloc(A*sizeof(double)); 	// folded arrray for omega
	
	for (int i=0; i < A; i++) {
		w[i] = 0;
	}

	band_mat bmat;
  	long nbands_low = A; /* The matrix is block tridiagonal: one block on either side of the
                               blocks of diagonals */  
  	long nbands_up  = nbands_low;
  	init_band_mat(&bmat, nbands_low, nbands_up, A);

  	double ctime = 0.0;
	long j,p;

  	for(p=0; p<ny; p++) {
		for(j=0; j<nx; j++) {
			//double x = j*dx;
			//double y = p*dy;
			long ii = indx(j,p,ny,nx);
            fprintf(output,"%lf\n", w[ii]);
            //printf("%lf\n",w[ii]);
        }
    }

    long nstep = 0;

  	// iterate over time
  	while (ctime<tf) {
	  	 /* Loop over the equation number and set the matrix
	     values equal to the coefficients of the grid values 
	     note boundaries treated with special cases       */
		for(j=0; j<nx; j++) {
			for(p=0; p<ny; p++) {   
				long unknown_indx = indx(j,p,ny,nx);
				
				setv(&bmat,folding(unknown_indx, A),folding(indx(j-1,p,ny,nx), A),1.0/(dx*dx));	      		
				setv(&bmat,folding(unknown_indx, A),folding(indx(j+1,p,ny,nx), A),1.0/(dx*dx));	      		
				setv(&bmat,folding(unknown_indx, A),folding(indx(j,p-1,ny,nx), A),1.0/(dy*dy));	      		
				setv(&bmat,folding(unknown_indx, A),folding(indx(j,p+1,ny,nx), A),1.0/(dy*dy));
	      		setv(&bmat,folding(unknown_indx, A),folding(unknown_indx, A),-2.0/(dx*dx) - 2.0/(dy*dy));
			}
	  	}

	  	// Fold omega
	  	arrayfolded(foldedw, w, A, nx, ny);

	  	// now we can solve for phi
		solve_Ax_eq_b(&bmat, phi, foldedw);

		/*
		// unfold omega
		arrayunfolded(w, folded, A);
		*/

		// unfold phi
		arrayunfolded(unfoldedphi, phi, A, nx, ny);

		for(j=0; j<nx; j++) {
			for(p=0; p<ny; p++) {
				long unknown_indx = indx(j,p,ny,nx);
				phi[unknown_indx] = unfoldedphi[unknown_indx];
			}
		}

		
		// now that we know phi, we can find u
		for(j=0; j<nx; j++) {
			for(p=0; p<ny; p++) {
				long ii = indx(j,p,ny,nx);
				ux[ii] = ((phi[indx(j,p+1,ny,nx)]-phi[indx(j,p-1,ny,nx)])/(2*dx)) + u0x;
				uy[ii] = -((phi[indx(j+1,p,ny,nx)]-phi[indx(j-1,p,ny,nx)])/(2*dy)) + u0y;
			}
		}
		
		for(j=0; j<nx; j++) {
			for(p=0; p<ny; p++) {

				long ii = indx(j,p,ny,nx);
				// first define all of the derivatives in the equation for dw/dt:

				double dK_dx 	= (K[indx(j+1,p,ny,nx)] 	- K[indx(j-1,p,ny,nx)])/(2*dx);
				double duy_dx 	= (uy[indx(j+1,p,ny,nx)] 	- uy[indx(j-1,p,ny,nx)])/(2*dx);
				double dK_dy 	= (K[indx(j,p+1,ny,nx)] 	- K[indx(j,p-1,ny,nx)])/(2*dy);
				double dux_dy 	= (ux[indx(j,p+1,ny,nx)] 	- ux[indx(j,p-1,ny,nx)])/(2*dy);

				double d2w_dx2 	= (w[indx(j+1,p,ny,nx)] + w[indx(j-1,p,ny,nx)] - 2*w[ii])/(dx*dx);
				double d2w_dy2 	= (w[indx(j,p+1,ny,nx)] + w[indx(j,p-1,ny,nx)] - 2*w[ii])/(dy*dy);

				double dw_dx 	= (w[indx(j+1,p,ny,nx)] 	- w[indx(j-1,p,ny,nx)])/(2*dx);
				double dw_dy 	= (w[indx(j,p+1,ny,nx)] 	- w[indx(j,p-1,ny,nx)])/(2*dy);

				// plug into the full formula for w_next:

				w_next[ii] = dt*(- (dK_dx*uy[ii]) - (duy_dx*K[ii]) + (dK_dy*ux[ii]) 
					+ (dux_dy*K[ii]) + v*(d2w_dx2 + d2w_dy2) - ((dw_dx*ux[ii])+(dw_dy*uy[ii]))) + w[ii];
			}
		}
		
		//Output the results:
    	if (nstep%nsteps_per_output==0) {
			for(p=0; p<ny; p++) {
				for(j=0; j<nx; j++) {
					//double x = j*dx;
					//double y = p*dy;

					long ii = indx(j,p,ny,nx);
	                fprintf(output,"%lf\n", w_next[ii]);
	            }
	        }
		}
		nstep++;
		
		// swap the values in preparation for the next timestep iteration
		for (int jj=0;jj<A;jj++) {
			w[jj] = w_next[jj];
		}
		// move to next timestep
		ctime+=dt;
	}
	fclose(output);
	//
	
	freeforall(&bmat, w, w_next, phi, foldedw, unfoldedphi, ux, uy, K);
	return 1;
}


// Function to read the values for each of the parameters in input.txt
void read_input_data(double *Lx, double *Ly, int *nx, int *ny, double *tf, double *td, double *v, double *u0x, double *u0y) {
	FILE *input; //input.txt
	if(!(input=fopen("input.txt","r"))) {
		printf("Error opening parameters in input.txt\n");
		exit(1);
	}
	if(9!=fscanf(input,"%lf %lf %i %i %lf %lf %lf %lf %lf",Lx,Ly,nx,ny,tf,td,v,u0x,u0y)) {
       printf("Error reading parameters from input.txt\n");
       exit(1);
	}
	fclose(input);
}
void read_coeff_data(long A, double *K) {
	FILE *coeff; //coefficients.txt
	if(!(coeff=fopen("coefficients.txt","r"))) {
		printf("Error opening parameters in coefficients.txt\n");
		exit(1);
	}
	int j;
	for(j=0;j<A;j++) {
   		if (1!=fscanf(coeff,"%lf",K+j)) {
       		printf("Error reading parameters from coefficients.txt\n");
       		exit(1);
		}
	}
fclose(coeff);	
}

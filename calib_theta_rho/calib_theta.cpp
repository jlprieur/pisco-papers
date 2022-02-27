/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Program calib_theta 
* To compute the theta zero for speckle interferometry
* using images of stars which moves because of the Earth rotation 
*
* Assumes that the stars moves roughtly along the X axis
*
* JLP 
* Version 11-01-2017
-------------------------------------------------------------------*/
#define DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
// #include <jlp_ftoc.h>
#include "jlp_numeric.h"
#include "jlp_fitsio.h"

#define NCOEFF 2

/* Defined here: */
static int process_list_for_calib_theta(char *list_fname);
static int calib_theta_for_one_file(char *in_name, double *mean_theta0, 
                                    double *sigma_theta0);
static int find_all_centers(double *dble_image1, double *xcent, double *ycent,
                            INT4 nx, INT4 ny, INT4 idim);
static int compute_centers(double *dble_image1, double *xcent, double *ycent,
                           int *good_col, INT4 nx, INT4 ny, INT4 idim, 
                           double *phi, int *npts, int hwidth);
static int compute_matrices(double *AA, double *BB, double *xcent, 
                            double *ycent, int npts, int nc);
static int compute_residuals(double *xcent, double *ycent, int npts, 
                             double *phi, double *rms_resid);
static int reject_bad_columns(double *xcent, double *ycent, int nx, 
                              int *good_col, double *phi, double threshold);
static int fit_centers(double *AA, double *BB, double *xcent, 
                       double *ycent, double *phi, double *theta0,
                       double *rms_resid, int npts, int ncoeff);

int main(int argc, char *argv[])
{
int status;
char list_fname[128];

printf(" Program calib_theta to compute theta zero \n");
printf(" JLP Version 06-12-2017 \n");

/* One parameters only is allowed to run the program: */
/* Carefull: 7 parameters always, using JLP "runs" */
if(argc == 7 && *argv[3]) argc = 4;
if(argc == 7 && *argv[2]) argc = 3;
if(argc == 7 && *argv[1]) argc = 2;
if(argc != 2)
  {
  printf(" Syntax: calib_theta input_list\n"); 
  printf(" input_list: list of long integration image fits files\n");
  printf(" Fatal: Syntax error: argc=%d\n",argc);
  exit(-1);
  }

/* Input of parameters with the command line: */
 strcpy(list_fname, argv[1]);

 status = process_list_for_calib_theta(list_fname);

return(status);
}
/*************************************************************************
*
*************************************************************************/
static int process_list_for_calib_theta(char *list_fname)
{
int status, k, nfiles;
char in_fname[128];
double mean_theta0, sigma_theta0, mean_theta1, sigma_theta1;
double ssum, ssumsq, sumweights, ww;
FILE *fp;

if((fp = fopen(list_fname,"r")) == NULL) {
  printf("process_list_for_calib_theta/Error reading list of files >%s<\n",
          list_fname);
  return(-1);
  }

/* First go to count the files: */
k = 0;
ssum = 0.;
ssumsq = 0.;
sumweights = 0.;
while(!feof(fp)) {
  if(fscanf(fp, "%s\n", in_fname)) {
   k++;
   status  = calib_theta_for_one_file(in_fname, &mean_theta0, &sigma_theta0);
// DEBUG
   sigma_theta0 = 1.;
   ww = mean_theta0 / sigma_theta0;
   sumweights += 1. / sigma_theta0;
   ssum += ww;
   ssumsq += ww * ww;
   }
}
fclose(fp);
nfiles = k;

mean_theta1 = ssum / sumweights;
printf("mean value: %f (nfiles=%d)\n", mean_theta1, nfiles);
sigma_theta1 = sqrt(ssumsq / sumweights - mean_theta1 * mean_theta1);
printf("sigma: %f\n", sigma_theta1);

return(status);
}
/*************************************************************************
*
*************************************************************************/
static int calib_theta_for_one_file(char *in_name, double *mean_theta0, 
                                    double *sigma_theta0) 
{
double *AA, *BB, phi[NCOEFF];
double *xcent, *ycent; 
double theta0, rms_resid, old_rms_resid, threshold; 
int *good_col;
double *dble_image1;
int status, ncoeff, npts, hwidth, nx1, ny1, idim;
char comments1[81];
int i, imax;

AA = NULL; BB= NULL; ncoeff = NCOEFF; 
*mean_theta0 = 0.;
*sigma_theta0 = 0.;

printf("calib_theta_for_one_file: processing %s\n", in_name);

/**********************************************************/
  status = JLP_LoadFITSImage(in_name, comments1, &dble_image1, &nx1, &ny1);
  if(status != 0) {
     printf("calib_theta_for_one_file/Fatal error in JLP_LoadFITSImage status=%d\n", status); 
     exit(-1);
     }

/* Problem with the ICCD when nx = 384: last two columns are bad ! */
  idim = nx1;
  if(nx1 == 384) nx1 = 382;

  xcent = (double *)malloc(nx1 * sizeof(double));
  ycent = (double *)malloc(nx1 * sizeof(double));
/* Flags for good columns */
  good_col = (int *)malloc(nx1 * sizeof(double));
  for(i = 0; i < nx1; i++) good_col[i] = 1;

/* Matrices for least-squares: */
AA = (double *)malloc(nx1 * ncoeff * sizeof(double));
BB = (double *)malloc(nx1 * sizeof(double));

/* First iteration with all columns: */
  find_all_centers(dble_image1, xcent, ycent, nx1, ny1, nx1);
  status = fit_centers(AA, BB, xcent, ycent, phi, &theta0, &rms_resid,
                       nx1, ncoeff);

/* 6 more iterations: */
  if(!status) {

  imax = 6;
  for(i = 0; i < imax && !status; i++) {
   old_rms_resid = rms_resid;
/* Compute the center of gravity of each column within
* a sigment of size (2 * half_width): */
   hwidth = MAXI(3. * rms_resid, 6);
   compute_centers(dble_image1, xcent, ycent, good_col, nx1, ny1, idim, phi, 
                   &npts, hwidth);

/* Solve problem again */
   status = fit_centers(AA, BB, xcent, ycent, phi, &theta0, &rms_resid,
                        npts, ncoeff);

/* Select columns by removing all columns for which the residuals 
* are larger than 2.5 sigma (with a minimum value of 1 pixel): */
   threshold = MAXI(2.5 * rms_resid, 1.);
   reject_bad_columns(xcent, ycent, nx1, good_col, phi, threshold);
   }
  }

*mean_theta0 = theta0;
*sigma_theta0 = rms_resid;

  free(xcent);
  free(ycent);
if(AA != NULL) free(AA);
if(BB != NULL) free(BB);

return(0);
}
/******************************************************************
* find_all_centers
* Determine the center of the trace in all columns (first guess)
*
* INPUT:
* dble_image1[nx, ny]: input image
*
* OUTPUT:
* xcent, ycent: location of the maxima (along a column)
*******************************************************************/
static int find_all_centers(double *dble_image1, double *xcent, double *ycent,
                            INT4 nx, INT4 ny, INT4 idim)
{
double sum, sumw;
double wmax;
int j1, j2, jmax, hwidth;
register int i, j;

/* Select a window of 2*hwidth+1 pixels: */
   hwidth = 8;

/* Do not take the whole range since bad pixels close to the edges */
j1 = ny / 3 ; j2 = (ny * 2) / 3;
/* For each column look for maximum inside the range determined by the previous
* iteration (to be more robust and do not be trapped on gamma ray events) 
*/
/* Loop on the columns */
for(i = 0; i < nx; i++) {

/* Localize the maximum along the column in the selected range: */
   wmax = dble_image1[i + j1 * idim];
   jmax = j1;
   for(j = j1 + 1; j < j2; j++) {
     if(dble_image1[i + j * idim] > wmax) {
       wmax = dble_image1[i + j * idim];
       jmax = j;
       }
     }

/* Compute more accurately the location of the maximum: */
   j1 = MAXI(jmax - hwidth, 0);
   j2 = MINI(jmax + hwidth + 1, ny);
   sum = 0; sumw = 0.; 
   for(j = j1; j < j2; j++) {
     sum += dble_image1[i + j * idim] * (double)j; 
     sumw += dble_image1[i + j * idim];
     } 
   ycent[i] = sum / sumw;
   xcent[i] = (double)i;
/* DEBUG
printf(" xcent[%d] = %.1f ycent = %.1f (j1=%d j2=%d jmax=%d wmax=%f)\n", 
         i, xcent[i], ycent[i], j1, j2, jmax, wmax);
*/
}
   
return(0);
}
/******************************************************************
* compute_centers
* Determine the center of gravity of the trace along the columns using the
* result of the coefficient fit as a first guess 
*
*
* INPUT:
* dble_image1[nx, ny]: input image
* phi: coefficients of the polynomial 
* hwidth: half width (select a window of 2*hwidth+1 pixels)
* good_col: flags set to one for good columns
*
* OUTPUT:
* xcent, ycent: location of the maxima (along a column)
* npts: number of good columns to be used for the next fit
*******************************************************************/
static int compute_centers(double *dble_image1, double *xcent, double *ycent,
                           int *good_col, INT4 nx, INT4 ny, INT4 idim, 
                           double *phi, int *npts, int hwidth)
{
double sum, sumw, yyc, backg;
double ww, wmax;
int j1, j2, jmax, jjc, jback1, jback2, nn;
register int i, j, k;

/* Loop on the columns */
k = 0;
for(i = 0; i < nx; i++) {

  if(good_col[i]) {
/* Location of the center using the polynomial */
   yyc = phi[0] - phi[1] * xcent[i];
   jjc = (int)(yyc + 0.5);

/* Select a window of 2*hwidth+1 pixels for computing the gravity center: */
   j1 = MAXI(jjc - hwidth, 0);
   j2 = MINI(jjc + hwidth + 1, ny);

/* Localize the maximum along the column in the selected range: */
   wmax = dble_image1[i + j1 * idim];
   jmax = j1;
   for(j = j1 + 1; j < j2; j++) {
     ww = dble_image1[i + j * idim];
     if(ww > wmax) {
       wmax = ww;
       jmax = j;
       }
     }
/* New window for the location of the maximum: */
   j1 = MAXI(jmax - hwidth, 0);
   j2 = MINI(jmax + hwidth + 1, ny);
/* Select a window of 10*hwidth+1 pixels for computing the backgound: */
   jback1 = MAXI(jmax - 5*hwidth, 0);
   jback2 = MINI(jmax + 5*hwidth + 1, ny);

/*
if(i < 10) 
  printf("jjc = %d jmax=%d wmax=%f backg=%f\n", jjc, jmax, wmax, backg);
*/

/* Compute the backgound along the column: */
   backg = 0.;
   nn = 0;
   for(j = jback1; j < j1; j++) {
     backg = dble_image1[i + j * idim];
     nn++;
     }
   for(j = 2; j < jback2; j++) {
     backg = dble_image1[i + j * idim];
     nn++;
     }
   backg /= (double)nn;

/* Compute more accurately the location of the maximum: */
   sum = 0; sumw = 0.; 
   for(j = j1; j < j2; j++) {
     ww = dble_image1[i + j * idim] - backg;
     sum += ww * (double)j; 
     sumw += ww;
     } 
   ycent[k] = sum / sumw;
   xcent[k] = (double)i;
   k++;
/* DEBUG
printf(" k=%d Column %d ycent = %.1f (j1=%d j2=%d jmax=%d wmax=%f)\n", 
         k, i, ycent[k], j1, j2, jmax, wmax);
*/
 } /* EOF test on good_col[i] */
}
   
*npts = k;
printf("compute_centers/Successful selection of %d columns (nx=%d)\n",
        *npts, nx);
return(0);
}
/*********************************************************************
* To compute the elements of matrices AA and BB 
* which correspond to the least-square problem:
*  AA^* AA phi = AA^* BB 
*
* with phi: polynomial of 1st order to be multiplied with x
*  phi(x,y) = a x  + c
*********************************************************************/
static int compute_matrices(double *AA, double *BB, double *xcent, 
                            double *ycent, int npts, int nc)
{
register int i;

/* Loop on the data points */
   for(i = 0; i < npts; i++) {
     AA[0 + i * nc] = 1;
     AA[1 + i * nc] = xcent[i];
     BB[i] = ycent[i];
   }

return(0);
}
/*********************************************************************
* Compute the residuals 
*********************************************************************/
static int compute_residuals(double *xcent, double *ycent, int npts, 
                             double *phi, double *rms_resid)
{
double sum, sumsq, ww;
register int i;

  sum = 0.;
  sumsq = 0;
  for(i = 0; i < npts; i++) {
   ww = ycent[i] - phi[0] - phi[1] * xcent[i];
   sum += ww;
   sumsq += ww * ww;
   }
sum /= (double)npts;
sumsq = sumsq/(double)npts - sum * sum;
sumsq = sqrt(sumsq);
/*
printf("compute_residuals/mean=%f sigma=%f\n", sum, sumsq);
*/
*rms_resid = sumsq;
return(0);
}  
/*********************************************************************
* Reject bad columns 
*********************************************************************/
static int reject_bad_columns(double *xcent, double *ycent, int nx, 
                              int *good_col, double *phi, double threshold)
{
double ww;
register int i, k;

  k = 0;
  for(i = 0; i < nx; i++) {
   ww = ABS(ycent[i] - phi[0] - phi[1] * xcent[i]);
   if(ww >= threshold && good_col[i]) {
       good_col[i] = 0;
       k++;
       } 
   }
printf("reject_bad_columns/removing %d columns\n", k);

return(0);
}  
/********************************************************************
*
********************************************************************/
static int fit_centers(double *AA, double *BB, double *xcent, 
                       double *ycent, double *phi, double *theta0,
                       double *rms_resid, int npts, int ncoeff)
{
int ifail, status;
register int i;

  *theta0 = 0.;
  *rms_resid = 0.;

  compute_matrices(AA, BB, xcent, ycent, npts, ncoeff);

/* Initial guess: */
  for(i = 0; i < ncoeff; i++) phi[i] = 0.;
/* Solve problem
    AA^* AA phi = AA^* BB
 with conjugate gradients:
*/
  status = JLP_CGRAD(AA, BB, phi, &ncoeff, &npts, &ifail);
  if(ifail != 0 || status) {
     printf("From JLP_CGRAD/status=%d ifail=%d\n", status, ifail);
     status = -2;
  } 
  else {
/* Output coefficients: b + ax */
   printf("  ycent = %f + %f * xcent\n", phi[0], phi[1]);
/* Compute residuals */
   compute_residuals(xcent, ycent, npts, phi, rms_resid);
   *theta0 =  atan(phi[1])*180./PI;
   printf("theta=%f degrees (rms_residuals=%f)\n", *theta0,
           *rms_resid);
  }

return(ifail);
}

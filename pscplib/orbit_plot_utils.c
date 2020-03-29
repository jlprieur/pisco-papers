/*************************************************************************
* orbit_plot_utils.c
* To plot the curves rho(epoch), theta(epoch) and XY orbits in the sky plane
*
*
* JLP
* Version 13/07/2018
**************************************************************************/
#include "orbit_plot_utils.h"
#include "residuals_utils.h"
#include "jlp_splot_idv_prototypes.h" // JLP_DEVICE_CURVE 
#include "jlp_string.h" // jlp_trim_string

static int orbit_plot_dx_or_dy(char *dxy_infile, char *comments, int npts_max, 
                               char *plotfile, int plot_dx, 
                               int smoothed_values);
static int cleanup_negative_values(float *xx, float *yy, int *npts);
static int orbit_rescale_theta(float *xplot, float *yplot, int npts_max, 
                               int *npts, int ncurves);

/*************************************************************************
* orbit_plot_rho
* To plot rho versus epoch
*
* INPUT:
* measures_infile: name of the input file
* comments: comments to be written in the plot
* orbit_infile: file with orbital elements
* npts_max: number of points (estimated from the number of lines in the input file)
* iformat: 1=raw (from WDS)\n");
*          2=corrected measures (from 1BIN.FOR of orbit_weight.c)\n");
*          3=measures and ephemerids (from residuals_2.c)\n");
*          4=raw or corrected measures and orbit
* iformat = 1: epoch, rho_O, theta_O, n_nights, author, aperture, instrument
* iformat = 2: epoch, rho_O, theta_O, n_nights, author, aperture, weight 
* iformat = 3: epoch, rho_O, rho_C, Drho_O-C, theta_O, theta_C, Dtheta, author
* iformat = 4: epoch, rho_O, theta_O, n_nights, author, aperture, ... 
* iformat = 5: epoch, Dx, Dy, n_nights, author, aperture, weight 
*
* plotfile: output file name 
*************************************************************************/
int orbit_plot_rho(char *measures_infile, char *comments, char *orbit_infile, 
                   int nber_of_orbits, int npts_max, int iformat, 
                   char *plotfile)
{
float *xplot, *yplot;
int grid, jlp_axes, icol_x, icol_y, status, ncurves, isize; 
int npts[2], orbit_format;
char xlabel[41], ylabel[41], title[81], plotdev[128];
double Omega_node[3], omega_peri[3], i_incl[3], e_eccent[3], T_periastron[3]; 
double Period[3], a_smaxis[3], mean_motion[3], orbit_equinox[3];
double c_tolerance[3], rho_c, theta_c, epoch, epoch1, epoch2;
register int i, k;

if(iformat < 1 || iformat > 4) {
  fprintf(stderr, "orbit_plot_rho/Error: iformat=%d not alllowed here!\n", 
          iformat);
  return(-1);
  }

isize = npts_max * sizeof(float);
if(iformat == 3 || iformat == 4) isize = 2 * isize;
xplot = (float *)malloc(isize);
yplot = (float *)malloc(isize);
ncurves = 1;

/* 
* iformat = 1: epoch, rho_O, theta_O, n_nights, author, aperture, instrument
* iformat = 2: epoch, rho_O, theta_O, n_nights, author, aperture, weight 
* iformat = 3: epoch, rho_O, rho_C, Drho_O-C, theta_O, theta_C, Dtheta, author
* iformat = 4: epoch, rho_O, theta_O, n_nights, author, aperture, ... 
* iformat = 5: epoch, Dx, Dy, n_nights, author, aperture, weight 
*/
icol_x = 1;
icol_y = 2;

/* Read epoch and rho_O (iformat from 1 to 4) */
status = orbit_read_data_from_file(measures_infile, icol_x, icol_y, npts_max,
                                   xplot, yplot, &npts[0]);

/* Negative values for rho correspond to flagged data (not for residuals!)*/
cleanup_negative_values(xplot, yplot, &npts[0]);

if(status || !npts[0]) {
  fprintf(stderr, "orbit_plot_rho/Error reading rho and epoch in >%s<\n", 
          measures_infile);
  return(-1);
  }

/* Read epoch and rho_C */
if(iformat == 3) {
  icol_x = 1;
  icol_y = 3;
  status = orbit_read_data_from_file(measures_infile, icol_x, icol_y, npts_max,
                                     &xplot[npts_max], &yplot[npts_max], 
                                     &npts[1]);
  if(status) ncurves = 1;
  else ncurves = 2;
/* Compute rho_C from orbital elements: */
  } else if(iformat == 4) {
/* Marco's format: */
   orbit_format = 1;
   status = read_orbital_elements_from_file(orbit_infile, orbit_format,
                           Omega_node, omega_peri, i_incl, e_eccent, 
                           T_periastron, Period, a_smaxis, mean_motion, 
                           orbit_equinox, nber_of_orbits);
/*  c_tolerance = smallest increment allowed in the iterative process
*                used for solving Kepler's equation
*/
   for(k = 0; k < nber_of_orbits; k++) 
   c_tolerance[k] = ABS(1.5E-5 * cos(i_incl[k])
                      / sqrt((1.0 + e_eccent[k])/(1.0 - e_eccent[k])));

  epoch1 = xplot[0];
  epoch2 = xplot[npts[0]-1];
  for(i = 0; i < npts_max; i++) {
/* Compute rho_C from orbital elements at equally spaced
* epochs along the range of epochs: */
    epoch = epoch1 + i * (epoch2 - epoch1)/(double)(npts_max-1);
    compute_ephemerid_of_multiple_system(nber_of_orbits, Omega_node, 
                      omega_peri, i_incl, e_eccent, T_periastron,
                      Period, a_smaxis, mean_motion, epoch, c_tolerance,
                      &theta_c, &rho_c);
    xplot[npts_max + i] = epoch; 
    yplot[npts_max + i] = rho_c; 
    }
    npts[1] = npts_max;
    ncurves = 2;
  }
strcpy(xlabel, "epoch [year]");
strcpy(ylabel, "rho [arcsec]");
strcpy(title, "");
grid = 0;
jlp_axes = 1;

 sprintf(plotdev,"landscape/%s_rho.ps", plotfile);
 sprintf(plotdev,"square/%s_rho.ps", plotfile);
 jlp_trim_string(plotdev,60);
 
orbit_plot_curves(xplot, yplot, npts_max, npts, ncurves, jlp_axes, grid,
                  xlabel, ylabel, title, measures_infile, comments, plotdev);

free(xplot);
free(yplot);
return(0);
}
/*************************************************************************
* orbit_plot_theta
* To plot theta versus epoch
*
* INPUT:
* measures_infile: name of the input file
* comments: comments to be written in the plot
* orbit_infile: file with orbital elements
* npts_max: number of points (estimated from the number of lines in the input file)
* iformat: 1=raw (from WDS)\n");
*          2=corrected measures (from 1BIN.FOR of orbit_weight.c)\n");
*          3=measures and ephemerids (from residuals_2.c)\n");
*          4=raw or corrected measures and orbit
* plotfile: output file name 
*************************************************************************/
int orbit_plot_theta(char *measures_infile, char *comments, char *orbit_infile, 
                     int nber_of_orbits, int npts_max, int iformat, char *plotfile)
{
float *xplot, *yplot;
int grid, jlp_axes, icol_x, icol_y, status, isize; 
int ncurves, npts[2], orbit_format;
char xlabel[41], ylabel[41], title[81], plotdev[128], outfile[128];
double Omega_node[3], omega_peri[3], i_incl[3], e_eccent[3], T_periastron[3];
double Period[3], a_smaxis[3], mean_motion[3], orbit_equinox[3];
double c_tolerance[3], rho_c, theta_c, epoch, epoch1, epoch2;
register int i, k;

if(iformat < 1 || iformat > 4) {
  fprintf(stderr, "orbit_plot_theta/Error: iformat=%d not alllowed here!\n", 
          iformat);
  return(-1);
  }

isize = npts_max * sizeof(float);
if(iformat == 3 || iformat == 4) isize = 2 * isize;
xplot = (float *)malloc(isize);
yplot = (float *)malloc(isize);
ncurves = 1;

/* 
* iformat = 1: epoch, rho_O, theta_O, n_nights, author, aperture, instrument
* iformat = 2: epoch, rho_O, theta_O, n_nights, author, aperture, weight 
* iformat = 3: epoch, rho_O, rho_C, Drho_O-C, theta_O, theta_C, Dtheta, author
* iformat = 4: epoch, rho_O, theta_O, n_nights, author, aperture, ... 
* iformat = 5: epoch, Dx, Dy, n_nights, author, aperture, weight 
*/
if(iformat != 3) {
  icol_x = 1;
  icol_y = 3;
  } else {
  icol_x = 1;
  icol_y = 5;
  }
/* Read epoch and theta_O */
status = orbit_read_data_from_file(measures_infile, icol_x, icol_y, npts_max,
                                   xplot, yplot, &npts[0]);
if(status || !npts[0]) {
  fprintf(stderr, "orbit_plot_theta/Error reading theta and epoch in >%s<\n", 
          measures_infile);
  return(-1);
  }
/* Read epoch and theta_C */
if(iformat == 3) {
  icol_x = 1;
  icol_y = 6;
  status = orbit_read_data_from_file(measures_infile, icol_x, icol_y, npts_max,
                                     &xplot[npts_max], &yplot[npts_max], 
                                     &npts[1]);
  if(status) ncurves = 1;
  else ncurves = 2;
/* Compute theta_C from orbital elements: */
  } else if(iformat == 4) {
/* Marco's format: */
   orbit_format = 1;
   status = read_orbital_elements_from_file(orbit_infile, orbit_format,
                           Omega_node, omega_peri, i_incl, e_eccent, 
                           T_periastron, Period, a_smaxis, mean_motion, 
                           orbit_equinox, nber_of_orbits);
/*  c_tolerance = smallest increment allowed in the iterative process
*                used for solving Kepler's equation
*/
   for(k = 0; k < nber_of_orbits; k++) 
   c_tolerance[k] = ABS(1.5E-5 * cos(i_incl[k])
                      / sqrt((1.0 + e_eccent[k])/(1.0 - e_eccent[k])));

  epoch1 = xplot[0];
  epoch2 = xplot[npts[0]-1];
  for(i = 0; i < npts_max; i++) {
/* Compute rho_C from orbital elements at equally spaced
* epochs along the range of epochs: */
    epoch = epoch1 + i * (epoch2 - epoch1)/(double)(npts_max-1);
/* Compute the ephemerids corresponding to the observation epoch: */
    compute_ephemerid_of_multiple_system(nber_of_orbits, Omega_node, 
                      omega_peri, i_incl, e_eccent, T_periastron,
                      Period, a_smaxis, mean_motion, epoch, c_tolerance,
                      &theta_c, &rho_c);
    xplot[npts_max + i] = epoch;
    yplot[npts_max + i] = theta_c;
    }
    npts[1] = npts_max;
    ncurves = 2;
  }

strcpy(xlabel, "epoch [year]");
strcpy(ylabel, "theta [degrees]");
strcpy(title, "");
grid = 0;
jlp_axes = 1;

/* Subtract 360 degrees if needed: */
 orbit_rescale_theta(xplot, yplot, npts_max, npts, ncurves);

/*
sprintf(plotdev,"landscape/%s_theta.ps", plotfile);
*/
sprintf(plotdev,"square/%s_theta.ps", plotfile);
jlp_trim_string(plotdev,60);
 
orbit_plot_curves(xplot, yplot, npts_max, npts, ncurves, jlp_axes, grid,
                  xlabel, ylabel, title, measures_infile, comments, plotdev);

sprintf(outfile,"%s_theta.dat", plotfile);
orbit_output_curves(xplot, yplot, npts_max, npts, ncurves, outfile);

free(xplot);
free(yplot);
return(0);
}
/*************************************************************************
* orbit_plot_Dx
* To plot Dx versus epoch
*
* INPUT:
* measures_infile: name of the input file
* comments: comments to be written in the plot
* npts_max: number of points (estimated from the number of lines in the input file)
* iformat = 5: epoch, Dx, Dy, n_nights, author, aperture, weight 
*
* plotfile: output file name 
*************************************************************************/
int orbit_plot_Dx(char *dxy_infile, char *comments, int npts_max, 
                  char *plotfile, int smoothed_values)
{
int plot_dx = 1, status;
char plotdev[128];

sprintf(plotdev,"landscape/%s_dx.ps", plotfile);
jlp_trim_string(plotdev,60);
 
printf("Plotdev= %s\n", plotdev);
status = orbit_plot_dx_or_dy(dxy_infile, comments, npts_max, plotdev, plot_dx,
                             smoothed_values);
return(status);
}

/*************************************************************************
* orbit_plot_Dy
* To plot Dy versus epoch
*
* INPUT:
* measures_infile: name of the input file
* comments: comments to be written in the plot
* npts_max: number of points (estimated from the number of lines in the input file)
* iformat = 5: epoch, Dx, Dy, n_nights, author, aperture, weight 
*
* plotfile: output file name 
*************************************************************************/
int orbit_plot_Dy(char *dxy_infile, char *comments, int npts_max, 
                  char *plotfile, int smoothed_values)
{
int plot_dx = 0, status;
char plotdev[128];

sprintf(plotdev,"landscape/%s_dy.ps", plotfile);
jlp_trim_string(plotdev,60);
 
printf("orbitplot_Dy/plotdev= %s\n", plotdev);
status = orbit_plot_dx_or_dy(dxy_infile, comments, npts_max, plotdev, plot_dx, 
                             smoothed_values);
return(status);
}
/*************************************************************************
*
*************************************************************************/
static int orbit_plot_dx_or_dy(char *dxy_infile, char *comments, int npts_max, 
                               char *plotfile, int plot_dx, int smoothed_values)
{
float *xplot, *yplot;
int remove_negative_values;
int grid, jlp_axes, icol_x, icol_y, status, ncurves, isize, npts[2];
char xlabel[41], ylabel[41], title[81], plotdev[128];

isize = npts_max * sizeof(float);
xplot = (float *)malloc(isize);
yplot = (float *)malloc(isize);
ncurves = 1;

/* 
* iformat = 5: epoch, Dx, Dy, n_nights, author, aperture, weight 
*/
/* Read epoch and Dx or Dy */
if(plot_dx) {
  sprintf(plotdev,"landscape/%s_dx.ps", plotfile);
  jlp_trim_string(plotdev,60);
  icol_x = 1;
  icol_y = 2;
  } else {
  sprintf(plotdev,"landscape/%s_dy.ps", plotfile);
  jlp_trim_string(plotdev,60);
  icol_x = 1;
  icol_y = 3;
  }

/* Read epoch and Dx or Dy */
if(smoothed_values) {
printf("WARNING: smoothed values! \n");
remove_negative_values = 0;
status = orbit_read_smoothed_data_from_file(dxy_infile, icol_x, icol_y, 
                                            npts_max, xplot, yplot, &npts[0], 
                                            remove_negative_values);
} else {
status = orbit_read_data_from_file(dxy_infile, icol_x, icol_y, npts_max,
                                   xplot, yplot, &npts[0]);
}

if(status || !npts[0]) {
  fprintf(stderr, "orbit_plot_dx_or_dy/Error reading dx/dy and epoch in >%s<\n", 
          dxy_infile);
  return(-1);
  }

strcpy(xlabel, "epoch [year]");
if(plot_dx)
  strcpy(ylabel, "Dx [arcsec]");
else
  strcpy(ylabel, "Dy [arcsec]");

strcpy(title, "");
grid = 0;
jlp_axes = 1;

printf("orbitplot_dx_or_dy/plotdev= %s\n", plotdev);
orbit_plot_curves(xplot, yplot, npts_max, npts, ncurves, jlp_axes, grid,
                  xlabel, ylabel, title, dxy_infile, comments, plotdev);

free(xplot);
free(yplot);
return(0);
}
/**************************************************************************
* Routine to plot one or two curves (yplot versus xplot)
*
* INPUT:
* xplot[ncurves*npts_max], yplot[ncurves*npts_max]: arrays to be plotted
* npts_max: maximum size of curve to be plotted
* npts[]: number of points of the curves contained in arrays xplot, yplot
* ncurves: number of curves
* jlp_axes: flag set to 1 if jlp_axes instead of mongo axes 
* grid: flag set to 1 if a grid is wanted for the plot
* xlabel, ylabel, title: caption of figure
* measures_infile, comments: information to be printed as a comment in small fonts
* plotdev: plotting device
**************************************************************************/
int orbit_plot_curves(float *xplot, float *yplot, int npts_max, int *npts,
                      int ncurves, int jlp_axes, int grid, char *xlabel, 
                      char *ylabel, char *title, char *measures_infile, 
                      char *comments, char *plotdev)
{
float errx[1], erry[1], xout[20], yout[20];
float xmin, xmax, ymin, ymax, expand;
float xmin_user, xmax_user, ymin_user, ymax_user;
int nout, error_bars, ticks_in, auto_scale, iplan;
int full_caption, xgrid, ygrid, xlog10, ylog10, idv;
int y_is_reversed = 0;
int status, nout_max = 20;
char nchar[8], pcolor[30], out_filename[64];

orbit_curves_min_max(xplot, yplot, npts_max, npts, ncurves, &xmin, &ymin, 
                     &xmax, &ymax);

/* Initialize plotting device: 
* plan: flag set to 1 if same scale in X and Y
*/
 iplan = 0;
 strcpy(out_filename, "tmp.ps");
 status = JLP_DEVICE_CURVE(plotdev, out_filename, &xmin, &xmax, &ymin, &ymax, 
                           &iplan, title, &idv);
if(status) {
   printf(" Fatal error opening graphic device: %s \n", plotdev);
   exit(-1);
   }

/* Display the curve: */
strcpy(pcolor,"Default");
error_bars = 0;
/* 8 = empty circles
* 2 = empty triangles 
* 3 = filled triangles 
* 4 = + 
* 5 = x 
* 510 seems OK
*/
strcpy(&nchar[0],"510");
strcpy(&nchar[4],"L");
full_caption = 0;
/* Same as NEWPLOT2 but with a better scale and the possibility of grids: */
xlog10 = 0; ylog10 = 0;
if(grid) {
 xgrid = 1; ygrid = 1; 
 } else {
 xgrid = 0; ygrid = 0; 
 }
expand = 1.2;
ticks_in = 0;
auto_scale = 0;
iplan = 1;
xmin_user = 0.;
xmax_user = 1.;
ymin_user = 0.;
ymax_user = 1.;
newplot210(xplot, yplot, errx, erry, npts, npts_max, ncurves,
          xmin_user, xmax_user, ymin_user, ymax_user, auto_scale, iplan,
          y_is_reversed,
          xlabel, ylabel, title, nchar, pcolor, xout, yout, &nout, nout_max,
          error_bars, measures_infile, comments, full_caption, jlp_axes, 
          xgrid, ygrid, xlog10, ylog10, expand, ticks_in, idv);

/* Close display device and free idv number: */
JLP_SPCLOSE(&idv);

return(0);
}
/***************************************************************************
* Parameters of the frame
* 
* INPUT:
* xx[npts_max], yy[npts_max]: arrays to be plotted
* npts: number of points of the curves contained in xx, yy arrays
* ncurves: number of curves
*
* OUTPUT:
* xmin,ymin,xmax,ymax: boundaries for the plot
***************************************************************************/
int orbit_curves_min_max(float *xx, float *yy, int npts_max, int *npts, 
                         int ncurves, float *xmin, float *ymin, float *xmax, 
                         float *ymax)
{
float w1;
register int i, k;

*xmin = xx[0]; *xmax = *xmin;
*ymin = yy[0]; *ymax = *ymin;
for(k = 0; k < ncurves; k++) {
  for(i = 1; i < npts[k]; i++) {
    *xmin = MINI(*xmin, xx[i + k * npts_max]);
    *xmax = MAXI(*xmax, xx[i + k * npts_max]);
    *ymin = MINI(*ymin, yy[i + k * npts_max]);
    *ymax = MAXI(*ymax, yy[i + k * npts_max]);
    }
  }

if(*xmin == *xmax) *xmax += 1.;
if(*ymin == *ymax) *ymax += 1.;

/* Plus or minus 5%: */
w1 = (*xmax - *xmin)/20;
*xmin -= w1;
*xmax += w1;

w1 = (*ymax - *ymin)/20;
*ymin -= w1;
*ymax += w1;

return(0);
}
/*************************************************************************
* Load (xplot, yplot) from columns of the input file
*
* INPUT:
* measures_infile: name of the file
* icol_x: column number of X values
* icol_y: column number of X values
* npts_max: number of points corresponding to the size of xplot, yplot 
*
* OUTPUT:
* xplot, yplot: arrays filled with input X, Y values
* npts: number of points read from input file 
*************************************************************************/
int orbit_read_data_from_file(char *measures_infile, int icol_x, int icol_y, 
                               int npts_max, float *xplot, float *yplot, 
                               int *npts)
{
float ww1, ww2;
int status1, status2;
char in_line[80];
FILE *fp_in;

if((fp_in = fopen(measures_infile, "r")) == NULL) {
  fprintf(stderr, "orbit_read_data_from_file/Fatal error opening input file >%s<\n",
          measures_infile);
  return(-1);
  }

*npts = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 80, fp_in)) {
/* Check if it is not a comment: */ 
    if(in_line[0] != '%') {

status1 = orbit_read_column_from_line(in_line, icol_x, &ww1);
status2 = orbit_read_column_from_line(in_line, icol_y, &ww2);
    if(status1 || status2) {
       fprintf(stderr, "orbit_read_data_from_file/Error reading line %d (end of file?)\n",
               (*npts) + 1);
       break;
     } else {
       if(*npts > npts_max -1) {
          fprintf(stderr, "orbit_read_data_from_file/Fatal error: npts_max=%d npts=%d\n",
                  npts_max, *npts);
          exit(-1);
          }
       xplot[*npts] = ww1; 
       yplot[*npts] = ww2;
       (*npts)++;
     }
    } /* EOF in_line != % */
  } /* EOF fgets */
 } /* EOF while */

fclose(fp_in);

return(0);
}
/************************************************************************
*
*
************************************************************************/
int orbit_read_measures_from_file(char *measures_infile, float *xplot, 
                                  float *yplot, int *npts0, int iformat,
                                  int npts_max)
{
int icol_x, icol_y, remove_negative_values, status;
/* 
* iformat = 1: epoch, rho_O, theta_O, n_nights, author, aperture, instrument
* iformat = 2: epoch, rho_O, theta_O, n_nights, author, aperture, weight 
* iformat = 3: epoch, rho_O, rho_C, Drho_O-C, theta_O, theta_C, Dtheta, author
* iformat = 4: epoch, rho_O, theta_O, n_nights, author, aperture, ... 
* iformat = 5: epoch, Dx, Dy, n_nights, author, aperture, weight 
*/
if(iformat != 3) {
 icol_x = 2;
 icol_y = 3;
 } else {
 icol_x = 2;
 icol_y = 5;
 }
/* Read Dx and Dy if iformat == 5 */
if(iformat == 5) {
  remove_negative_values = 0;
  status = orbit_read_smoothed_data_from_file(measures_infile, icol_x, icol_y, 
                                              npts_max, xplot, yplot, npts0,
                                              remove_negative_values);
/* Read rho_O and theta_O if iformat != 5 */
  } else {
  status = orbit_read_data_from_file(measures_infile, icol_x, icol_y, npts_max,
                                     xplot, yplot, npts0);
  }

/* DEBUG:
printf("DEBUG/orbit_plot_skyplane/x0=%f y0=%f npts=%d npts_max=%d\n", 
        xplot[0], yplot[0], npts[0], npts_max);
*/
/* Handling error: */
if(status || *npts0 == 0) {
  fprintf(stderr, "orbit_plot_skyplane/Error reading rho_O and theta_O in >%s<\n", 
          measures_infile);
  status = -1; 
  }

return(status);
}
/*************************************************************************
* Load (xplot, yplot) from columns of the input file
*
* INPUT:
* measures_infile: name of the file
* icol_x: column number of X values
* icol_y: column number of X values
* npts_max: number of points corresponding to the size of xplot, yplot 
*
* OUTPUT:
* xplot, yplot: arrays filled with input X, Y values
* npts: number of points read from input file 
*************************************************************************/
int orbit_read_smoothed_data_from_file(char *measures_infile, int icol_x, 
                                       int icol_y, int npts_max, float *xplot, 
                                       float *yplot, int *npts, 
                                       int remove_negative_values)
{
float ww1, ww2, sum1, sum2, epoch, epoch0, epoch_bin_width;
int status1, status2, nvalues;
char in_line[80];
FILE *fp_in;

if((fp_in = fopen(measures_infile, "r")) == NULL) {
  fprintf(stderr, "orbit_read_smoothed_data_from_file/Fatal error opening input file >%s<\n",
          measures_infile);
  return(-1);
  }

/* Two years for averaging the measurements: */
epoch_bin_width = 2.;

*npts = 0;
epoch0 = 0.;
sum1 = 0.;
sum2 = 0.;
nvalues = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 80, fp_in)) {
/* Check if it is not a comment: */ 
    if(in_line[0] != '%') {

orbit_read_column_from_line(in_line, 1, &epoch);
/* Compute mean of values within one year: */
if(epoch > epoch0 + epoch_bin_width) {
  if(nvalues > 0) {
       xplot[*npts] = sum1 /(float)nvalues; 
       yplot[*npts] = sum2 /(float)nvalues; 
       (*npts)++;
     }
  sum1 = 0.;
  sum2 = 0.;
  nvalues = 0;
  epoch0 = epoch;
  }
status1 = orbit_read_column_from_line(in_line, icol_x, &ww1);
status2 = orbit_read_column_from_line(in_line, icol_y, &ww2);
    if(status1 || status2) {
       fprintf(stderr, "orbit_read_smoothed_data_from_file/Error reading line %d (end of file?)\n",
               (*npts) + 1);
       break;
     } else {
       if(*npts > npts_max -1) {
          fprintf(stderr, "orbit_read_smoothed_data_from_file/Fatal error: npts_max=%d npts=%d\n",
                  npts_max, *npts);
          exit(-1);
          }
     }
     if((sum1 < 0 || sum2 < 0) && remove_negative_values) {
     } else {
     sum1 += ww1;
     sum2 += ww2;
     nvalues++;
     }
    } /* EOF in_line != % */
  } /* EOF fgets */
 } /* EOF while */
fclose(fp_in);

return(0);
}
/*************************************************************************
* Estimate the maximum value of the number of points contained in a file
* (necessary for further memory allocation)
*
* INPUT:
* measures_infile: name of the file
*
* OUTPUT:
* npts_max: number of lines (discarding commented lines)
*************************************************************************/
int orbit_get_npts_max(char *measures_infile, int *npts_max)
{
char in_line[80];
FILE *fp_in;

if((fp_in = fopen(measures_infile, "r")) == NULL) {
  fprintf(stderr, "orbit_get_npts_max/Fatal error opening input file >%s<\n",
          measures_infile);
  return(-1);
  }

*npts_max = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 80, fp_in)) {
/* Check if it is not a comment: */ 
    if(in_line[0] != '%') {
    (*npts_max)++;
    } /* EOF in_line != % */
  } else {
  break;
  } /* EOF fgets */
 } /* EOF while */
fclose(fp_in);

if(*npts_max == 0) {
  fprintf(stderr, "orbit_get_npts_max/empty file: npts_max=0 !\n");
  return(-1);
  } else {
  printf("orbit_get_npts_max/npts_max=%d\n", *npts_max);
  }

return(0);
}
/**************************************************************************
* Curve xx=epoch yy=rho
* Negative values for rho correspond to bad flagged data 
*
**************************************************************************/
static int cleanup_negative_values(float *xx, float *yy, int *npts)
{
int old_npts, k;
register int i;

old_npts = *npts;
k = 0;
for(i = 0; i < old_npts; i++) {
  if(yy[i] > 0) {
     xx[k] = xx[i];
     yy[k] = yy[i];
     k++;
     }
  } 
*npts = k;

if(*npts == 0) {
  fprintf(stderr, "cleanup_negative_values/Fatal error: npts = 0!\n");
  exit(-1);
  } else {
  printf("DDEBUG: after cleanup: npts=%d (before: npts=%d)\n", *npts, old_npts);
  }
return(0);
}
/*************************************************************************
* To read float value in a given column,
* assuming that the line is a series of float values separated with blanks.
*
* INPUT:
* inline: character string to be read
* icol: column number
*
* OUTPUT:
* value: value of column #icol
**************************************************************************/
int orbit_read_column_from_line(char *in_line, int icol, float *value)
{
int status, nval;

  switch (icol) {
    case 1:
      nval = sscanf(in_line, "%f", value);
      break;
    case 2:
      nval = sscanf(in_line, "%*f %f", value);
      break;
    case 3:
      nval = sscanf(in_line, "%*f %*f %f", value);
      break;
    case 4:
      nval = sscanf(in_line, "%*f %*f %*f %f", value);
      break;
    case 5:
      nval = sscanf(in_line, "%*f %*f %*f %*f %f", value);
      break;
    case 6:
      nval = sscanf(in_line, "%*f %*f %*f %*f %*f %f", value);
      break;
    case 7:
      nval = sscanf(in_line, "%*f %*f %*f %*f %*f %*f %f", value);
      break;
    case 8:
      nval = sscanf(in_line, "%*f %*f %*f %*f %*f %*f %*f %f", value);
      break;
    default:
      fprintf(stderr, "orbit_read_column_from_line/Fatal error : icol=%d is not allowed!\n", icol);
      exit(-1);
      break;
   }

if(nval == 1) status = 0;
 else status = -1;

return(status);
}
/*************************************************************************
* To read float value in a given column,
* assuming that the line is a series of float values separated with blanks.
*
* INPUT:
* inline: character string to be read
* icol: column number
*
* OUTPUT:
* value: value of column #icol
**************************************************************************/
int dble_read_column_from_line(char *in_line, int icol, double *value)
{
int status, nval;

  switch (icol) {
    case 1:
      nval = sscanf(in_line, "%lf", value);
      break;
    case 2:
      nval = sscanf(in_line, "%*lf %lf", value);
      break;
    case 3:
      nval = sscanf(in_line, "%*lf %*lf %lf", value);
      break;
    case 4:
      nval = sscanf(in_line, "%*lf %*lf %*lf %lf", value);
      break;
    case 5:
      nval = sscanf(in_line, "%*lf %*lf %*lf %*lf %lf", value);
      break;
    case 6:
      nval = sscanf(in_line, "%*lf %*lf %*lf %*lf %*lf %lf", value);
      break;
    case 7:
      nval = sscanf(in_line, "%*lf %*lf %*lf %*lf %*lf %*lf %lf", value);
      break;
    case 8:
      nval = sscanf(in_line, "%*lf %*lf %*lf %*lf %*lf %*lf %*lf %lf", value);
      break;
    default:
      fprintf(stderr, "read_column_from_line/Fatal error : icol=%d is not allowed!\n", icol);
      exit(-1);
      break;
   }

if(nval == 1) status = 0;
 else status = -1;

return(status);
}
/**************************************************************************
* Routine to output one or two curves (yplot versus xplot) to file
*
* INPUT:
* xplot[ncurves*npts_max], yplot[ncurves*npts_max]: arrays to be plotted
* npts_max: maximum size of curve to be plotted
* npts[]: number of points of the curves contained in arrays xplot, yplot
* ncurves: number of curves
* outfile: output file name
*
**************************************************************************/
int orbit_output_curves(float *xplot, float *yplot, int npts_max, int *npts, 
                        int ncurves, char *outfile)
{
FILE *fp_out;
register int i;

if((fp_out = fopen(outfile, "w")) == NULL) {
  fprintf(stderr, "orbit_output_curves/Fatal error opening output file >%s<\n",
          outfile);
  return(-1);
  }

 fprintf(fp_out,"%%ncurves=%d npts1=%d, npts2=%d\n", ncurves, npts[0], npts[1]);

/* First curve: */
 fprintf(fp_out,"%% Curve #1\n");
 fprintf(fp_out,"%d \n", npts[0]);
for(i = 0; i < npts[0]; i++) 
   fprintf(fp_out, "%f %f\n", xplot[i], yplot[i]);

/* Second curve: */
if(ncurves > 1) {
  fprintf(fp_out,"%% Curve #2\n");
  fprintf(fp_out,"%d \n", npts[1]);
  for(i = 0; i < npts[1]; i++) 
   fprintf(fp_out, "%f %f\n", xplot[npts_max + i], yplot[npts_max + i]);
  }

fclose(fp_out);
return(0);
}
/**************************************************************************
* Routine to rescale theta values to avoid "steps" of 360 degrees
*
* INPUT:
* xplot[ncurves*npts_max], yplot[ncurves*npts_max]: arrays to be plotted
* npts_max: maximum size of curve to be plotted
* npts[]: number of points of the curves contained in arrays xplot, yplot
* ncurves: number of curves
*
**************************************************************************/
static int orbit_rescale_theta(float *xplot, float *yplot, int npts_max, 
                               int *npts, int ncurves)
{
register int i, j, k;

printf("\n orbit_rescale_theta: ncurves=%d\n", ncurves);

/* First iteration with increasing X values: */ 
 for(k = 0; k < ncurves; k++) {
  for(i = 0; i < npts[k]; i++) {
/* Makes a diagnostic on each curve: */
    if(yplot[i + k * npts_max] - yplot[i-1 + k * npts_max] > 300) {
printf("Iteration 1: Discontinuity in theta found for curve #%d at pixel #%d\n", k, i);
/* Then solves the problem: */
     for(j = i; j < npts[k]; j++) yplot[j + k * npts_max] -= 360.; 
     i = npts[k];
     }
   }
 }

/* Second iteration with decreasing X values: */ 
 for(k = 0; k < ncurves; k++) {
  for(i = npts[k] - 1; i > 1; i--) {
/* Makes a diagnostic on each curve: */
    if(yplot[i - 1 + k * npts_max] - yplot[i + k * npts_max] > 300) {
printf("Iteration 2: Discontinuity in theta found for curve #%d at pixel #%d\n", k, i);
/* Then solves the problem: */
     for(j = i-1; j >= 0; j--) yplot[j + k * npts_max] -= 360.; 
     i = 0;
     }
   }
 }

return(0);
}

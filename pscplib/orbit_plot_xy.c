/*************************************************************************
* orbit_plot_xy.c
* To plot the curves rho(epoch), theta(epoch) and XY orbits in the sky plane
*
*
* JLP
* Version 13/01/2018
**************************************************************************/
#include "orbit_plot_utils.h"
#include "residuals_utils.h"
#include "jlp_splot_idv_prototypes.h" // JLP_DEVICE_CURVE 
#include "jlp_string.h"       // jlp_trim_string 

static int compute_resid_vectors(char *measures_infile, float *xstart, 
                                 float *ystart, float *xend, float *yend, 
                                 int npts_max, int *nresid);
static int load_pisco_measures(char *measures_infile, float *xplot, 
                                float *yplot, int npts_max, int *npts_pisco); 
static int draw_resid_vectors(float *xstart, float *ystart, float *xend, 
                              float *yend, int nresid, int idv);
static int draw_sense_of_motion(int ix1, int iy1, double rad0, 
                               int north_to_east, int idv);
static int draw_cross_and_north_east_label(float data_range, int north_to_east,
                                           char *plot_title, int idv);
static int draw_line_of_apsids(char *orbit_infile, int *north_to_east, 
                               int draw_apsids, int idv);

/*************************************************************************
* orbit_plot_skyplane
* Plot the orbit in the plane of the sky
* Full orbit if iformat==4, or part of the orbit only in other cases 
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
* iplot:  1=plot rho vs epoch and theta vs epoch
*         2=plot dx vs epoch and dy vs epoch
*         3=plot XY_full_orbit 
*         4=plot XY_part_of_orbit
* resid_vectors: flag set to 1 if residual vectors are to be plotted
* plotfile: output file name 
* plot_title: title of the plot 
*************************************************************************/
int orbit_plot_skyplane(char *measures_infile, char *comments, 
                        char *orbit_infile, int nber_of_orbits, int npts_max, 
                        int iformat, int iplot, int resid_vectors, 
                        char *plotfile, char *plot_title)
{
float *xplot, *yplot, rho, theta, epoch;
float *rho_resid, *theta_resid, *xstart, *ystart, *xend, *yend;
int icol_x, icol_y, status, k, npts[3], ncurves, isize, npts0, nresid;
int orbit_format, ncur_max, draw_apsids, npts_pisco;
char xlabel[41], ylabel[41], title[81], plotdev[128];
double Omega_node[3], omega_peri[3], i_incl[3], e_eccent[3], T_periastron[3];
double Period[3], a_smaxis[3], mean_motion[3], orbit_equinox[3];
double c_tolerance[3], rho_c, theta_c;
register int i, ic;

isize = npts_max * sizeof(float);
ncur_max = 1;
// 3 curves: 1=data 2=orbit 3=pisco_data
if((iplot == 3) || (iplot == 4)) ncur_max = 3; 
isize *= ncur_max;
xplot = (float *)malloc(isize);
yplot = (float *)malloc(isize);

status = orbit_read_measures_from_file(measures_infile, xplot, yplot, &npts0, 
                                       iformat, npts_max);
if(status != 0) return(-1);
  
npts[0] = npts0;
ncurves = 1;

/* Read rho_C and theta_C (part of the orbit only if iformat==3)*/
  if((iformat == 3) && (iplot == 4)) {
   icol_x = 3;
   icol_y = 6;
   status = orbit_read_data_from_file(measures_infile, icol_x, icol_y, npts_max,
                                      &xplot[npts_max], &yplot[npts_max], 
                                      &npts[1]);
   if(status || !npts[1]) {
     fprintf(stderr, "orbit_plot_skyplane/Error reading rho_C and theta_C in >%s<\n", 
            measures_infile);
     return(-1);
     }

   ncurves++;
/* Compute rho_C and theta_C from orbital elements at equally spaced
* epochs along the period (full orbit): */
  } else if(iplot == 3) {
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

   for(i = 0; i < npts_max; i++) {
/* Compute rho_C from orbital elements at equally spaced
* epochs along the period: */
    epoch = T_periastron[0] + i * Period[0]/(double)(npts_max-1);
/* Compute the ephemerids corresponding to the observation epoch: */
    compute_ephemerid_of_multiple_system(nber_of_orbits, Omega_node, 
                      omega_peri, i_incl, e_eccent, T_periastron,
                      Period, a_smaxis, mean_motion, epoch, c_tolerance,
                      &theta_c, &rho_c);
    xplot[i + npts_max] = rho_c;
    yplot[i + npts_max] = theta_c;
    }
   npts[1] = npts_max;
   ncurves++;
}

// Load PISCO measures:
  if((iplot == 3) || (iplot == 4)) {
    status = load_pisco_measures(measures_infile, &xplot[npts_max * 2], 
                                 &yplot[npts_max * 2], npts_max, &npts_pisco); 
    if((status == 0) && (npts_pisco > 0)){ 
      npts[2] = npts_pisco;
      ncurves++;
      }
    }

// Compute the residual vectors:
nresid = 0;
if(resid_vectors && (iformat == 3)) {
  xstart = (float *)malloc(npts_max * sizeof(float));
  xend = (float *)malloc(npts_max * sizeof(float));
  ystart = (float *)malloc(npts_max * sizeof(float));
  yend = (float *)malloc(npts_max * sizeof(float));
  status = compute_resid_vectors(measures_infile, xstart, ystart, xend, yend, 
                                 npts_max, &nresid);
  if(status != 0) {
    fprintf(stderr, "compute_resid_vectors/Error cannot plot the residuals: status=%d ! \n", 
            status);
    nresid = 0;
    }
  } else if(resid_vectors && (iformat != 3)) {
    fprintf(stderr, "Error cannot plot the residuals if iformat!=3 ! \n");
  }

/* Transformation with a change of origin for theta (zero at the bottom of the 
* plot), except for iformat==5 (residuals with Drho, Dtheta)
*/
if(iformat != 5) {
  for(ic = 0; ic < ncurves; ic++) { 
    k = 0;
    for(i = 0; i < npts[ic]; i++) { 
      rho = xplot[i + ic * npts_max]; 
      theta = yplot[i + ic * npts_max]; 
      theta -= 90.;
/* Conversion to radians: */
      theta *= (PI/180.);
/* Discard all negative rho's */
      if(rho > 0) {
// Convert [rho, theta] to [ rho cos(theta), rho sin(theta) ] for ic=0,1 :
        xplot[k + ic * npts_max] = rho * cos(theta);
        yplot[k + ic * npts_max] = rho * sin(theta);
        k++;
      }
    } /* EOF loop on i */
    npts[ic] = k;
  } /* EOF loop on ic */
} /* EOF case iformat != 5 */

strcpy(xlabel, "X (arcsec)");
strcpy(ylabel, "Y (arcsec)");
strcpy(title, "");

sprintf(plotdev,"square/%s_orbit.ps", plotfile);
jlp_trim_string(plotdev,60);

/* Draw line of apsids (Epoch of periastron, and opposite at +Period/2) */
if((iplot == 3) && (orbit_infile[0] != '\0')) {
 draw_apsids = 1;
 } else {
 draw_apsids = 0;
 }
 
/* DEBUG
printf("orbit_plot_skyplane/calling orbit_plot_XY with x0=%f y0=%f npts=%d npts_max=%d ncurves=%d icurve=0\n", 
        xplot[0], yplot[0], npts[0], npts_max, ncurves);
printf("orbit_plot_skyplane/calling orbit_plot_XY with x0=%f y0=%f npts=%d npts_max=%d ncurves=%d icurve=1\n", 
        xplot[1], yplot[1], npts[1], npts_max, ncurves);
if(ncurves > 2) {
 printf("orbit_plot_skyplane/calling orbit_plot_XY with x0=%f y0=%f npts=%d npts_max=%d ncurves=%d icurve=2\n", 
        xplot[2], yplot[2], npts[2], npts_max, ncurves);
 }
*/
orbit_plot_XY(xplot, yplot, npts_max, npts, ncurves, xlabel, ylabel, title,
              measures_infile, comments, orbit_infile, plotdev,
              draw_apsids, plot_title, xstart, ystart, 
              xend, yend, nresid);

free(xplot);
free(yplot);
return(0);
}
/*************************************************************************
*
* Assume here that iformat == 3
*
*************************************************************************/
static int compute_resid_vectors(char *measures_infile, float *xstart, 
                                 float *ystart, float *xend, float *yend, 
                                 int npts_max, int *nresid)
{
int i, k, icol_x, icol_y, npts_start, npts_end; 
int status_start, status_end, status = -1;
double rho, theta;

// iformat = 3: epoch, rho_O, rho_C, Drho_O-C, theta_O, theta_C, Dtheta_O-C, author
// Read rho_O, theta_O (to compute later xend, yend)
  icol_x = 2;
  icol_y = 5;
  status_end = orbit_read_data_from_file(measures_infile, icol_x, icol_y, 
                                         npts_max, xend, yend, &npts_end);
// iformat = 3: epoch, rho_O, rho_C, Drho_O-C, theta_O, theta_C, Dtheta_O-C, author
// Read rho_C, theta_C (to compute later xstart, ystart)
  icol_x = 3;
  icol_y = 6;
  status_start = orbit_read_data_from_file(measures_infile, icol_x, icol_y, 
                                           npts_max, xstart, ystart, 
                                           &npts_start);
  if((status_start == 0) && (status_end == 0) && (npts_start == npts_end)) {
    *nresid = npts_start;
    status = 0;
    } else {
    fprintf(stderr, "compute_resid_vectors/Error reading the residuals, npts_start=%d npts_end=%d : will not plot them !\n",
            npts_start, npts_end);
    *nresid = 0;
    }

/* DEBUG
 printf("BEFORE/compute_resid_vectors/nresid=%d\n", *nresid);
 for(i = 0; i < 10; i++) {
   printf("i=%d C:xstart=%f ystart=%f O:xend=%f yend=%f\n", 
           i, xstart[i], ystart[i], xend[i], yend[i]);
   }
*/

// Compute xstart, ystart, xend, yend:
  k = 0;
  for(i = 0; i < *nresid; i++) {
// Valid data only when rho > 0:
    if((xstart[i] > 0) && (xend[i] > 0)) {
// rho_C, theta_C:
/* Transformation with a change of origin for theta (zero at the bottom of the 
plot) */
      rho = xstart[i];
      theta = ystart[i] - 90.;
/* Conversion to radians: */
      theta *= (PI/180.);
      xstart[k] = rho * cos(theta);
      ystart[k] = rho * sin(theta);
// rho_O, theta_O:
      rho = xend[i];
      theta = yend[i] - 90.;
/* Conversion to radians: */
      theta *= (PI/180.);
      xend[k] = rho * cos(theta);
      yend[k] = rho * sin(theta);
      k++;
    } // EOF if rho > 0
  } // EOF loop on i

*nresid = k;

printf("ZZZ/compute_nresid: nresid = %d\n", *nresid);

return(status);
}
/**************************************************************************
* Assume here that iformat == 3
*
*************************************************************************/
static int load_pisco_measures(char *measures_infile, float *xplot, 
                                float *yplot, int npts_max, int *npts)
{
float ww1, ww2;
int status1, status2, icol_x, icol_y;
char in_line[80];
FILE *fp_in;

if((fp_in = fopen(measures_infile, "r")) == NULL) {
  fprintf(stderr, "load_pisco_measures/Fatal error opening input file >%s<\n",
          measures_infile);
  return(-1);
  }

// Assume here that iformat == 3
// iformat = 3: epoch, rho_O, rho_C, Drho_O-C, theta_O, theta_C, Dtheta_O-C, author
// Read rho_O, theta_O 
  icol_x = 2;
  icol_y = 5;

*npts = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 80, fp_in)) {
/* Check if it is not a comment: */
    if(in_line[0] != '%') {

// Read observer in cols 63, 64, 65
   if(!strncmp(&in_line[63], "Sca", 3) 
     || !strncmp(&in_line[63], "Pru", 3)) { 
printf("ZZZZ: inline=%s", in_line);
    status1 = orbit_read_column_from_line(in_line, icol_x, &ww1);
    status2 = orbit_read_column_from_line(in_line, icol_y, &ww2);
    if(status1 || status2) {
       fprintf(stderr, "load_pisco_measures/Error reading line %d (end of file?)\n",
               (*npts) + 1);
       break;
     } else {
       if(*npts > npts_max -1) {
          fprintf(stderr, "load_pisco_measures/Fatal error: npts_max=%d npts=%d\n",
                  npts_max, *npts);
          exit(-1);
          }
       xplot[*npts] = ww1;
       yplot[*npts] = ww2;
       (*npts)++;
     }
    } // EOF "Sca" or "Pru" 
    } /* EOF in_line != % */
  } /* EOF fgets */
 } /* EOF while */

fclose(fp_in);

printf("load_pisco_measures/Nber of PISCO observations found in file %s npts=%d\n", 
       measures_infile, *npts);
return(0);
}
/**************************************************************************
* Routine to plot an orbit in the sky plane
*
* INPUT:
* xplot[ncurves*npts_max], yplot[ncurves*npts_max]: arrays to be plotted
* npts[]: number of points of the curves contained in arrays xplot, yplot
* xlabel, ylabel, title: caption of figure
* measures_infile, comments: information to be printed as a comment in small fonts
* plotdev: plotting device
**************************************************************************/
int orbit_plot_XY(float *xplot, float *yplot, int npts_max, int *npts, 
                  int ncurves, char *xlabel, char *ylabel, char *title,
                  char *measures_infile, char *comments, char *orbit_infile,
                  char *plotdev, int draw_apsids, char *plot_title,
                  float *xstart, float *ystart, float *xend, float *yend,
                  int nresid)
{
float errx[1], erry[1], xout[20], yout[20];
float xmin, xmax, ymin, ymax, cross_width, expand;
float xmin_user, xmax_user, ymin_user, ymax_user;
float xmin_data, xmax_data, ymin_data, ymax_data, xdata_range, ydata_range;
float scale_fact, data_range;
int nout, error_bars, plan, ticks_in;
int auto_scale, iplan, full_caption, idv;
int status, nout_max = 20, north_to_east;
int y_is_reversed = 0;
char nchar[16], pcolor[4*32], out_filename[64];

orbit_curves_min_max(xplot, yplot, npts_max, npts, ncurves, 
                     &xmin, &ymin, &xmax, &ymax);

/* Initialize plotting device: 
* plan: flag set to 1 if same scale in X and Y
*/
 plan = 1;
 strcpy(out_filename, "tmp.ps");
 status = JLP_DEVICE_CURVE(plotdev, out_filename, 
                           &xmin, &xmax, &ymin, &ymax, &plan,
                           title, &idv);
if(status) {
   printf(" Fatal error opening graphic device: %s \n", plotdev);
   exit(-1);
   }

/* Display the curve: */
error_bars = 0;
/* 
* 2 = empty triangles 
* 3 = filled triangles 
* 4 = + 
* 5 = x 
* 510 seems OK
* 8 = empty circles
* 9 = filled circles
*/
// icurve= 0 : data points
// strcpy(&nchar[0],"510");
strcpy(&nchar[0],"510");

// icurve= 1 : orbit 
// strcpy(&nchar[4],"L1");
strcpy(&nchar[4],"L0");

// icurve= 2 : pisco measures 
strcpy(&nchar[8],"920");

// Colors: Red, Green, Black, etc
strcpy(&pcolor[0],"Green");
strcpy(&pcolor[32],"Black");
strcpy(&pcolor[32*2],"Red");

full_caption = 0;
expand = 0.8;
ticks_in = 1;
auto_scale = 0;
iplan = 1;

jlp_splot_min_max_for_curves(xplot, yplot, errx, erry, npts, npts_max, ncurves,
                       &xmin_data, &xmax_data, &ymin_data, &ymax_data, 
                       error_bars, iplan);
xdata_range = xmax_data - xmin_data;
ydata_range = ymax_data - ymin_data;
data_range = MAXI(xdata_range, ydata_range);

// Margin around the orbit as a fraction of the maximum width of the orbit
// scale_fact : fraction of the plot used for the orbit:
scale_fact = 0.15;

xmin_user = xmin_data - 1. * data_range * scale_fact;
xmax_user = xmax_data + 0.9 * data_range * scale_fact;
ymin_user = ymin_data - 0.7 * data_range * scale_fact;
ymax_user = ymax_data + 1.6 * data_range * scale_fact;

/* DEBUG
printf("DEBUG/plot_XY/x0=%f y0=%f npts=%d npts_max=%d\n", 
        xplot[0], yplot[0], npts[0], npts_max);
*/
NEWPLOT2(xplot, yplot, errx, erry, npts, &npts_max, &ncurves,
         &xmin_user, &xmax_user, &ymin_user, &ymax_user, &auto_scale, &iplan,
         &y_is_reversed,
         xlabel, ylabel, title, nchar, pcolor, xout, yout, &nout, &nout_max,
         &error_bars, measures_infile, comments, &full_caption, &expand,
         &ticks_in, &idv);

/* Draw line of apsids (Epoch of periastron, and opposite at +Period/2) */
if(orbit_infile[0] != '\0') {
   draw_line_of_apsids(orbit_infile, &north_to_east, draw_apsids, idv);
/* Draw central cross and North-East label: */
   draw_cross_and_north_east_label(data_range, north_to_east, plot_title, idv);
  }

// Draw residual vectors
draw_resid_vectors(xstart, ystart, xend, yend, nresid, idv);

/* Close display device and free idv number: */
JLP_SPCLOSE(&idv);

return(0);
}
/***********************************************************************
* Draw residual vectors
*
***********************************************************************/
static int draw_resid_vectors(float *xstart, float *ystart, float *xend, 
                              float *yend, int nresid, int idv)
{
float x1, x2, y1, y2;
int i, rr, gg, bb;

// Return directly if nresid is negative or null:
if(nresid <= 0) return(-1);

printf("DEBUG/draw_resid_vectors/ nresid=%d\n", nresid);

// Dark green:
   rr = 0;
   bb = 0;
   gg = 180;
   JLP_SETCOLOR(&rr, &gg, &bb, &idv);

// Draw residual vectors:
  for(i = 0; i < nresid; i++) {
    x1 = xstart[i];
    y1 = ystart[i];
    x2 = xend[i];
    y2 = yend[i];
    JLP_LINE1(&x1, &y1, &x2, &y2, &idv);
   }

// Black:
   rr = 0;
   bb = 0;
   gg = 0;
   JLP_SETCOLOR(&rr, &gg, &bb, &idv);

return(0);
}
/**********************************************************************
*
* Draw central cross and North-East label: 
*
* INPUT:
* cross_width (in arcseconds)
************************************************************************/
static int draw_cross_and_north_east_label(float data_range, int north_to_east,
                                           char *plot_title, int idv)
{
int ix1, ix2, iy1, iy2, max_length, idrawit;
float x1, x2, y1, y2, angle, expand, length, cross_width, ne_width;
double rad0;
char xlabel[20];

/* Draw central cross (to do the same as "2binpl.for"): */
/* User coordinates (arcseconds with orbit center at 0,0) */
cross_width = data_range * 0.1;
x1 = -cross_width/2.; x2 = cross_width/2.; y1 = 0.; y2 = 0.;
JLP_LINE1(&x1, &y1, &x2, &y2, &idv);
x1 = 0.; x2 = 0.; y1 = -cross_width/2.; y2 = cross_width/2.;
JLP_LINE1(&x1, &y1, &x2, &y2, &idv);

/* MGO coordinates */
expand = 0.6;
if(expand <= 0.8) {
  ne_width = 1600;
  } else {
  ne_width = 1800;
  }
// bottom right:
// ix1 = 30000 - ne_width; iy1 = 9000 + ne_width; 
// top left:
ix1 = 8600; iy1 = 27600 - ne_width; 
JLP_RELOC(&ix1, &iy1, &idv);
ix2 = ix1 + ne_width;
iy2 = iy1;
JLP_DRAW(&ix2, &iy2, &idv);
/*
void JLP_SPLABEL(char *xlabel, int *max_length, int *ix, int *iy,
                 float *angle, float *expand, int *idrawit, float *length,
                 int *idv1)
*/
max_length = 1;
strcpy(xlabel, "E");
ix2 += 50;
iy2 -= 300;
angle = 0.;
idrawit = 1;
JLP_SPLABEL(xlabel, &max_length, &ix2, &iy2, &angle, &expand, &idrawit,
            &length, &idv);


// Draw title at nearly the same height
ix2 = ix1 + 5000;
iy2 = iy1 - 200;
JLP_SPLABEL(plot_title, &max_length, &ix2, &iy2, &angle, &expand, &idrawit,
            &length, &idv);

JLP_RELOC(&ix1, &iy1, &idv);
ix2 = ix1;
iy2 = iy1 - ne_width;
JLP_DRAW(&ix2, &iy2, &idv);
if(expand <= 0.8) {
  ix2 -= 500;
  iy2 -= 800;
  } else {
  ix2 -= 700;
  iy2 -= 1000;
  }
strcpy(xlabel, "N");
JLP_SPLABEL(xlabel, &max_length, &ix2, &iy2, &angle, &expand, &idrawit,
            &length, &idv);

// Draw an arc of circle:
rad0 = ne_width * 1.1 + expand * 600.;
draw_sense_of_motion(ix1, iy1, rad0, north_to_east, idv);

return(0);
}
/******************************************************************
*
******************************************************************/
static int draw_sense_of_motion(int ix1, int iy1, double rad0, 
                               int north_to_east, int idv)
{
int i, ix0, iy0, ixstart, iystart;;
double angle0, angle1;

 angle0 = - (double)PI / 2. + PI / 8.; 
 for(i = 0; i < 100; i++) {
  angle1 = angle0 + (double)i * PI / (2. * 180.);
  ix0 = (int)((double)ix1 + rad0 * cos(angle1));
  iy0 = (int)((double)iy1 + rad0 * sin(angle1));
  if(i > 0)  JLP_DRAW(&ix0, &iy0, &idv);
  JLP_RELOC(&ix0, &iy0, &idv);
  }

// Draw arrow from north to east (at the top):
 if(north_to_east == 1) {
  ixstart = ix0;
  iystart = iy0;
  JLP_RELOC(&ixstart, &iystart, &idv);
  ix0 -= 600; 
  iy0 -= 400;
  JLP_DRAW(&ix0, &iy0, &idv);
  JLP_RELOC(&ixstart, &iystart, &idv);
  ix0 += 800; 
  iy0 -= 200;
  JLP_DRAW(&ix0, &iy0, &idv);
// Draw arrow from east to north (at the bottom:
  } else {
  ixstart = (int)((double)ix1 + rad0 * cos(angle0));
  iystart = (int)((double)iy1 + rad0 * sin(angle0));
  JLP_RELOC(&ixstart, &iystart, &idv);
  ix0 = ixstart;
  iy0 = iystart;
  ix0 += 300; 
  iy0 += 600;
  JLP_DRAW(&ix0, &iy0, &idv);
  JLP_RELOC(&ixstart, &iystart, &idv);
  ix0 += 400; 
  iy0 -= 600;
  JLP_DRAW(&ix0, &iy0, &idv);
  }

return(0);
}
/******************************************************************
* Draw line of apsids 
*
* INPUT:
* orbit_infile
********************************************************************/
static int draw_line_of_apsids(char *orbit_infile, int *north_to_east, 
                               int draw_apsids, int idv)
{
int status, orbit_format, nber_of_orbits, rr, gg, bb;
int lwidth, ltype;
float x1, x2, y1, y2, epoch;
double Omega_node[3], omega_peri[3], i_incl[3], e_eccent[3], T_periastron[3];
double Period[3], a_smaxis[3], mean_motion[3], orbit_equinox[3];
double c_tolerance[3], rho_c1, theta_c1, rho_c2, theta_c2;
double rho_c1b, theta_c1b;

orbit_format = 1;

// Plot only the line of apsids of the first orbit:
nber_of_orbits = 1;
status = read_orbital_elements_from_file(orbit_infile, orbit_format,
                           Omega_node, omega_peri, i_incl, e_eccent, 
                           T_periastron, Period, a_smaxis, mean_motion, 
                           orbit_equinox, nber_of_orbits);

/* c_tolerance = smallest increment allowed in the iterative process
*                used for solving Kepler's equation
*/
 c_tolerance[0] = ABS(1.5E-5 * cos(i_incl[0])
                      / sqrt((1.0 + e_eccent[0])/(1.0 - e_eccent[0])));

/* Look for the position of the companion at Periastron: 
*/
 epoch = T_periastron[0];
 compute_ephemerid_of_multiple_system(nber_of_orbits, Omega_node, 
                      omega_peri, i_incl, e_eccent, T_periastron,
                      Period, a_smaxis, mean_motion, epoch, c_tolerance,
                      &theta_c1, &rho_c1);
// Small difference of epochs to have a small difference in the angles too..
 epoch = T_periastron[0] + 0.1;
 compute_ephemerid_of_multiple_system(nber_of_orbits, Omega_node, 
                      omega_peri, i_incl, e_eccent, T_periastron,
                      Period, a_smaxis, mean_motion, epoch, c_tolerance,
                      &theta_c1b, &rho_c1b);
// Look for the sense of motion of the companion:
 if(theta_c1b - theta_c1 > 0) 
   *north_to_east = 1;
 else 
   *north_to_east = 0;

 epoch = T_periastron[0] + Period[0]/2.;
 compute_ephemerid_of_multiple_system(nber_of_orbits, Omega_node, 
                      omega_peri, i_incl, e_eccent, T_periastron,
                      Period, a_smaxis, mean_motion, epoch, c_tolerance,
                      &theta_c2, &rho_c2);

/* Transformation with a change of origin for theta (zero at the bottom of the 
* plot)
*/
 theta_c1 -= 90.;
 theta_c2 -= 90.;


/* Draw line of apsids (to do the same as "4binpl.for"): */
/* User coordinates (arcseconds) */
 if(draw_apsids) {
   x1 = rho_c1 * cos(theta_c1 * PI / 180.); 
   y1 = rho_c1 * sin(theta_c1 * PI / 180.); 
   x2 = rho_c2 * cos(theta_c2 * PI / 180.); 
   y2 = rho_c2 * sin(theta_c2 * PI / 180.); 
// Medium gray (128/256 = 0.5):
   rr = 128;
   bb = 128;
   gg = 128;
   JLP_SETCOLOR(&rr, &gg, &bb, &idv);
   lwidth = 1;
   ltype = 1;
   JLP_SETLINEPARAM(&lwidth, &ltype, &idv);
   JLP_LINE1(&x1, &y1, &x2, &y2, &idv);
   rr = 0;
   bb = 0;
   gg = 0;
   JLP_SETCOLOR(&rr, &gg, &bb, &idv);
   lwidth = 0;
   ltype = 0;
   JLP_SETLINEPARAM(&lwidth, &ltype, &idv);
   }

return(0);
}

/*************************************************************************
* orbit_plot_utils.h
* Prototypes of routines included in "orbit_plot_utils.c"
* To plot the curves rho(epoch), theta(epoch) and XY orbits in the sky plane
*
* JLP
* Version 09/01/2012
**************************************************************************/
#ifndef _orbit_plot_utils_h /* BEOF sentry */
#define _orbit_plot_utils_h 
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <jlp_catalog_utils.h>

#ifndef MAXI
#define MAXI(a,b) ((a) < (b)) ? (b) : (a)
#endif
#ifndef MINI
#define MINI(a,b) ((a) > (b)) ? (b) : (a)
#endif

/* Declaring linkage specification to have "correct names"
* that can be linked with C programs */

#ifdef __cplusplus
extern "C" {
#endif

int orbit_plot_rho(char *measures_infile, char *comments, 
                   char *orbit_file, int nber_of_orbits, int npts_max, 
                   int iformat, char *plotfile);
int orbit_plot_theta(char *measures_infile, char *comments, 
                     char *orbit_file, int nber_of_orbits, int npts_max, 
                     int iformat, char *plotfile);
int orbit_plot_curves(float *xplot, float *yplot, int npts_max, int *npts, 
                      int ncurves, int jlp_axes, int grid, char *xlabel, 
                      char *ylabel, char *title, char *measures_infile, 
                      char *comments, char *plotdev);
int orbit_plot_Dx(char *dxy_infile, char *comments, int npts_max, 
                  char *plotfile, int smoothed_values);
int orbit_plot_Dy(char *dxy_infile, char *comments, int npts_max, 
                  char *plotfile, int smoothed_values);
int orbit_plot_XY(float *xplot, float *yplot, int npts_max, int *npts, 
                  int ncurves, char *xlabel, char *ylabel, char *title,
                  char *measures_infile, char *comments, char *orbit_file,
                  char *plotdev, int draw_absids, char *plot_title,
                  float *xstart, float *ystart, float *xend, float *yend,
                  int nresid);
int orbit_get_npts_max(char *measures_infile, int *npts_max);
int orbit_output_curves(float *xplot, float *yplot, int npts_max, int *npts, 
                        int ncurves, char *outfile);
int orbit_curves_min_max(float *xx, float *yy, int npts_max, int *npts,
                         int ncurves, float *xmin, float *ymin, float *xmax,
                         float *ymax);
int orbit_read_column_from_line(char *in_line, int icol, float *value);
int orbit_read_measures_from_file(char *measures_infile, float *xplot,
                                  float *yplot, int *npts0, int iformat,
                                  int npts_max);
int orbit_read_data_from_file(char *measures_infile, int icol_x, int icol_y,
                              int npts_max, float *xplot, float *yplot,
                              int *npts);
int orbit_read_smoothed_data_from_file(char *measures_infile, int icol_x,
                                       int icol_y, int npts_max, float *xplot,
                                       float *yplot, int *npts,
                                       int remove_negative_values);

// orbit_plot_xy
int orbit_plot_skyplane(char *measures_infile, char *comments, 
                        char *orbit_file, int nber_of_orbits, int npts_max, 
                        int iformat, int iplot, int resid_vectors, 
                        char *plotfile, char *plot_title);

#ifdef __cplusplus
}
#endif

#endif

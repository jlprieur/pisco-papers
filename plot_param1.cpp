/************************************************************************
* "plot_param1.c"
*
* To plot with a parameter file, for plain files and latex tables 
*
* JLP 
* Version 21/02/2019
*************************************************************************/
//#include "jlp_catalog_utils.h"
#include "jlp_plot1.h"      // JLP_Plot1 class 
#include "jlp_splot_idv.h"  // newplot110() 
#include "astrom_utils1.h"  // latex_read_fvalue()

#define DEBUG

int read_latex_table_file(char *in_latex_fname, int icol_x, int icol_y,
                          int icol_xerr, int icol_yerr, float **x_plot,
                          float **y_plot, float **xerr_plot, float **yerr_plot,
                          int *npts, int *error_bars);
int plot_latex_data(float *xplot, float *yplot, float *errx, 
                    float *erry, int *npoints, int nmax, int error_bars);

int main(int argc, char *argv[])
{
char in_latex_fname[128], in_plot_param_fname[128];
int icol_x, icol_y, icol_xerr, icol_yerr;
float *x_plot, *y_plot, *xerr_plot, *yerr_plot;
int npoints[4], npts, error_bars, nmax;
FILE *fp_in_param;
JLP_Plot1 *jlp_plt0;

if(argc != 2) {
  printf("Syntax: plot_param1 param_file\n");
  return(-1);
}
strcpy(in_plot_param_fname, argv[1]);

printf("OK: in_plot_param_fname=%s\n", in_plot_param_fname);

jlp_plt0 = new JLP_Plot1(in_plot_param_fname);

/* Scan the input latex table
*/
icol_x = 1;
icol_y = 2;
icol_xerr = 3;
icol_yerr = 4;
read_latex_table_file(in_latex_fname, icol_x, icol_y, icol_xerr, 
                      icol_yerr, &x_plot, &y_plot, &xerr_plot, &yerr_plot,
                      &npts, &error_bars); 
nmax = npts;
npoints[0] = npts;

/*
plot_latex_data(x_plot, y_plot, xerr_plot, yerr_plot, npoints, 
                nmax, error_bars);
*/
return(0);
}
/************************************************************
*
****************************************************************/
/*************************************************************
* Draw the curve (called by UpdateBackupDC)
*************************************************************/
int plot_latex_data(float *xplot, float *yplot, float *errx, 
                    float *erry, int *npoints, int nmax, int error_bars)
{
int jlp_axes_are_wanted, xgrid_is_wanted, ygrid_is_wanted;
int x_is_log10, y_is_log10, ticks_in;
float expand;
int ncurves, iplan, y_is_reversed;
char plotdev[128], xlabel[40], ylabel[40], title[40];
char nchar[4*4], pcolor[30*4];
int i, k;

strcpy(plotdev,"square");
strcpy(title,"Title");
strcpy(xlabel,"X axis");
strcpy(ylabel,"Y axis");

ncurves = 1;

// Set default color and line type for the curves:
for(k = 0; k < ncurves; k++) {
  sprintf(&nchar[k*4], "L%.1d", k);
  strcpy(&pcolor[k*30],"Default");
}

/* in jlp_splot_idv/jlp_newplot.cpp : 
int jlp_newplot110(float *xplot, float *yplot, float *errx, float *erry,
                   int *npoints, int nmax, int ncurves, int iplan, 
                   int y_is_reversed, char *xlabel, char *ylabel, 
                   char *title, char *nchar, char *pcolor, int error_bars, 
                   char *plotdev, int jlp_axes_are_wanted,
                   int xgrid_is_wanted, int ygrid_is_wanted,
                   int x_is_log10, int y_is_log10, float expand,
                   int ticks_in)
*/
 iplan = 0;
 y_is_reversed = 0;
 error_bars = 0;
 jlp_axes_are_wanted = 0;
 xgrid_is_wanted = 0;
 ygrid_is_wanted = 0;
 x_is_log10 = 0;
 y_is_log10 = 0;
 expand = 1.6;
 ticks_in = 0;
 jlp_newplot110(xplot, yplot, errx, erry, npoints, nmax, ncurves, iplan, 
                y_is_reversed, xlabel, ylabel, title, nchar, pcolor, 
                error_bars, plotdev, jlp_axes_are_wanted,
                xgrid_is_wanted, ygrid_is_wanted, x_is_log10, y_is_log10, 
                expand, ticks_in);

return(0);
}
/************************************************************************
* Scan the latex table 
*
* INPUT:
*  in_latex_fname: filename of the input file containing the latex table 
*  icol_x : number of the column containing x_plot
*  icol_y : number of the column containing y_plot
*  icol_xerr : number of the column containing xerr_plot
*  icol_yerr : number of the column containing yerr_plot
*  x_plot: array of the X data to be plotted
*  y_plot: array of the Y data to be plotted
*  xerr_plot: array of the X error data to be plotted
*  yerr_plot: array of the Y error data to be plotted
*
* OUTPUT:
*  npts : number of points (curve x_plot, y_plot)
*  error_bars : flag set to one if error_bars
*
*************************************************************************/
int read_latex_table_file(char *in_latex_fname, int icol_x, int icol_y,
                          int icol_xerr, int icol_yerr, float **x_plot,
                          float **y_plot, float **xerr_plot, float **yerr_plot,
                          int *npts, int *error_bars) 
{
char in_line[256];
int iline, ipts, stat_x, stat_y, stat_xerr, stat_yerr, nlines, iverbose = 1;
double ww_x, ww_y, ww_xerr, ww_yerr;
FILE *fp_in_latex;

// First reading of the file to determine for the total number of lines 
// (i.e. maximum size of output arrays)
/* Open the input latex table: */
  if((fp_in_latex = fopen(in_latex_fname, "r")) == NULL) {
    fprintf(stderr, "plot_latex/Fatal error opening in latex file: %s\n",
           in_latex_fname);
    return(-1);
   }

iline = 0;
while(!feof(fp_in_latex)) {
  if(fgets(in_line, 256, fp_in_latex)) {
    iline++;
  }
}
nlines = iline;
#ifdef DEBUG
 printf("First step: latex file succesfully read: nlines = ", nlines);
#endif

fclose(fp_in_latex);
*x_plot = new float[nlines];
*y_plot = new float[nlines];
*xerr_plot = new float[nlines];
*yerr_plot = new float[nlines];
*error_bars = 0;

// Initialize error values to zero:
for(ipts = 0; ipts < nlines; ipts++) {
  (*xerr_plot)[ipts] = 0.;
  (*yerr_plot)[ipts] = 0.;
 }

// Second reading of the file to read the data: 
/* Re-open the input latex table: */
  if((fp_in_latex = fopen(in_latex_fname, "r")) == NULL) {
    fprintf(stderr, "plot_latex/Fatal error opening in latex file: %s\n",
           in_latex_fname);
    return(-1);
   }
iline = 0;
ipts = 0;
while(!feof(fp_in_latex)) {
  if(fgets(in_line, 256, fp_in_latex)) {
// Read the x and y values in columns icol_x and icol_y:
      stat_x = latex_read_fvalue(in_line, &ww_x, icol_x, iverbose);
      stat_y = latex_read_fvalue(in_line, &ww_y, icol_y, iverbose);
// Read the x and y error values in columns icol_xerr and icol_yerr:
      stat_xerr = latex_read_fvalue(in_line, &ww_xerr, icol_xerr, iverbose);
      stat_yerr = latex_read_fvalue(in_line, &ww_yerr, icol_yerr, iverbose);
      if((stat_x == 0) && (stat_y == 0)){
        (*x_plot)[ipts] = (float)ww_x;
        (*y_plot)[ipts] = (float)ww_y;
        if(stat_xerr == 0) {
           *error_bars = 1;
           (*xerr_plot)[ipts] = (float)ww_xerr;
           }
        if(stat_yerr == 0) {
           *error_bars = 1;
           (*yerr_plot)[ipts] = (float)ww_yerr;
           }
        ipts++;
      } else {
        fprintf(stderr,"Warning: error reading value: %s (line=%d, icol_x=%d)\n", 
                in_line, iline, icol_x);
      } // If status 
    iline++;
   } // If fgets 
 } /* EOF while ... */
#ifdef DEBUG
 printf("Second step: from input_latex, %d lines sucessfully read\n", iline);
#endif

*npts = ipts;

fclose(fp_in_latex);
return(0);
}

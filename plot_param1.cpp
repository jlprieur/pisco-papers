/************************************************************************
* "plot_param1.c"
*
* To plot with a parameter file, for plain files and latex tables 
*
* JLP 
* Version 31/03/2020
*************************************************************************/
//#include "jlp_catalog_utils.h"
#include "jlp_plot1.h"      // JLP_Plot1 class 
#include "jlp_splot_idv.h"  // newplot110() 
#include "jlp_read_latex_tables.h"  // jlp_read_latex_table_to_float 

#define DEBUG

int plot_latex_data(float *xplot, float *yplot, float *errx, 
                    float *erry, int *npoints, int nmax, int error_bars);

int main(int argc, char *argv[])
{
char in_latex_fname[128], in_plot_param_fname[128];
int icol_x, icol_y, icol_xerr, icol_yerr;
float *x_plot, *y_plot, *xerr_plot, *yerr_plot;
int npoints[4], npts, npts_size, error_bars, nmax;
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
strcpy(in_latex_fname, "test_latex.tex");
jlp_read_latex_table_to_float(in_latex_fname, icol_x, icol_y, icol_xerr, 
                      icol_yerr, &x_plot, &y_plot, &xerr_plot, &yerr_plot,
                      &npts_size, &npts, &error_bars); 
nmax = npts;
npoints[0] = npts;

/*
plot_latex_data(x_plot, y_plot, xerr_plot, yerr_plot, npoints, 
                nmax, error_bars);
*/
return(0);
}
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

/*************************************************************************
* Program latex_add_magnitude
* 
* To process and update Latex table of red dwarfs with WDS and GAIA data 
*
\hline 
00014+3937 & HLD60 & K0V+K1V & 9.09 & 9.77 & \nodata & 20.42 & 1.91 & 19.336 & 0.020 & 5.16 & & & \\
*
* JLP
* Version 22/01/2022
*************************************************************************/
#include <stdio.h>
#include <stdlib.h> // exit(-1)
#include <string.h>
#include <math.h>  // log10()
#include <ctype.h>  // isalpha(), isdigit()
#include "latex_utils.h" // latex_get_column_item(), latex_remove_column()
#include "jlp_string.h"  // (in jlplib/jlp_fits/ ) jlp_cleanup_..

#include "latex_utils.h"  // latex_get_column_item()
//#include "astrom_utils1.h" 
//#include "astrom_utils2.h" 

#define MAX_LENGTH 256 

#define DEBUG

static int latex_read_parallax(char *in_line, int in_absp_col, int in_plx1_col, 
                               int in_plx2_col, int in_err_plx1_col, 
                               int in_err_plx2_col, double *plx0, 
                               double *err_plx0, double *absp0); 
static int latex_add_magnitude(FILE *fp_in, FILE *fp_out,
                               int in_mag_col, int in_absp_col,
                               int in_plx1_col, int in_plx2_col,
                               int in_err_plx1_col, int in_err_plx2_col,
                               int out_absmag_col);
static int compute_abs_mag(double gmag, double absp, double plx, double *Gmag);

int main(int argc, char *argv[])
{
int i, status, in_mag_col, nval;
int in_plx1_col, in_plx2_col, out_abspmag_col; 
int in_err_plx1_col, in_err_plx2_col, in_absp_col; 
char filein[128], fileout[128];
FILE *fp_in, *fp_out;

/* If command line with "runs" */
if(argc == 7){
 if(*argv[5]) argc = 6;
 else if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
 else argc = 1;
 }

if(argc != 4)
  {
  printf(" Syntax: latex_add_magnitude in_latex_file out_latex_file in_mag_col,in_absp_col,in_plx1_col,in_err_plx1_col,in_plx2_col,in_err_plx2_col,out_abspmag_col\n");
  exit(-1);
  }
else
  {
  strcpy(filein,argv[1]);
  strcpy(fileout,argv[2]);
  nval = sscanf(argv[3], "%d,%d,%d,%d,%d,%d,%d", &in_mag_col,
                 &in_absp_col, &in_plx1_col, &in_err_plx1_col, &in_plx2_col, 
                 &in_err_plx2_col, &out_abspmag_col);
  if(nval != 7) {
    fprintf(stderr, "Error, should have 7 column numbers !\n");
    exit(-1);
    }
  }

printf(" OK: filein=%s fileout=%s \n", filein, fileout);
printf("in_mag_col=%d in_absp_col=%d \n", in_mag_col, in_absp_col);
printf("in_plx1_col=%d in_plx2_col=%d out_abspmag_col=%d \n", 
        in_plx1_col, in_plx2_col, out_abspmag_col);
printf("in_err_plx1_col=%d in_err_plx2_col=%d \n", 
        in_err_plx1_col, in_err_plx2_col);

if((fp_in = fopen(filein,"r")) == NULL) {
  printf(" Fatal error opening input file %s \n",filein);
  exit(-1);
  }

if((fp_out = fopen(fileout,"w")) == NULL) {
  printf(" Fatal error opening output file %s \n",fileout);
  fclose(fp_in);
  exit(-1);
  }

latex_add_magnitude(fp_in, fp_out, in_mag_col, in_absp_col,
                   in_plx1_col, in_plx2_col, in_err_plx1_col,
                   in_err_plx2_col, out_abspmag_col);

fclose(fp_in);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input latex file 
************************************************************************/
static int latex_read_parallax(char *in_line, int in_absp_col, int in_plx1_col, 
                               int in_plx2_col, int in_err_plx1_col, 
                               int in_err_plx2_col, double *plx0, 
                               double *err_plx0, double *absp0) 
{
int status, iline, i, verbose_if_error;
double plx1_value, plx2_value, err_plx1_value, err_plx2_value;
char buffer[256];

*absp0 = 0.;
*plx0 = 0.;
*err_plx0 = 0.;
verbose_if_error = 0;

// absp_mag in in_absp_col:
   status = latex_get_column_item(in_line, buffer, in_absp_col, 
                                  verbose_if_error);
   if(status == 0) {
     sscanf(buffer, "%lf", absp0);
     }

// parallax in in_plx1_col:
   status = latex_get_column_item(in_line, buffer, in_plx1_col, 
                                  verbose_if_error);
   if(status == 0) {
     sscanf(buffer, "%lf", &plx1_value);
     }
// parallax in in_plx2_col:
   status = latex_get_column_item(in_line, buffer, in_plx2_col, 
                                  verbose_if_error);
   if(status == 0) {
     sscanf(buffer, "%lf", &plx2_value);
     }
// parallax error in in_err_plx1_col:
   status = latex_get_column_item(in_line, buffer, in_err_plx1_col, 
                                  verbose_if_error);
   if(status == 0) {
     sscanf(buffer, "%lf", &err_plx1_value);
     }
// parallax error in in_err_plx2_col:
   status = latex_get_column_item(in_line, buffer, in_err_plx2_col, 
                                  verbose_if_error);
   if(status == 0) {
     sscanf(buffer, "%lf", &err_plx2_value);
     }

// Compute parallax from plx1_value and plx2_value:
   if(plx1_value > 0.001) {
     *plx0 = plx1_value;
     }
   if(plx2_value > 0.001) {
     *plx0 = plx2_value;
     }
   if((plx1_value > 0.001) && (plx2_value > 0.001)) {
     *plx0 = (plx1_value + plx2_value)/2.;
     }
   if(err_plx1_value > 0.001) {
     *err_plx0 = err_plx1_value;
     }
   if(err_plx2_value > 0.001) {
     *err_plx0 = err_plx2_value;
     }
   if((err_plx1_value > 0.001) && (err_plx2_value > 0.001)) {
     *err_plx0 = 1. / (1. / err_plx1_value + 1. / err_plx2_value);
     }
   if((plx1_value > 0.001) && (err_plx1_value > 0.001) && (plx2_value > 0.001)
     && (err_plx2_value > 0.001)) {
     *plx0 = (plx1_value / err_plx1_value + plx2_value / err_plx2_value) 
              / (1. / err_plx1_value + 1. / err_plx2_value);
     }
#ifdef DEBUG
   printf("OK: in_line=%s \n", in_line);
#endif
   printf("OK: plx1_value=%f plx2_value=%f ", 
          plx1_value, plx2_value);
   printf("absp0=%f plx0=%f err_plx0=%f \n", 
          *absp0, *plx0, *err_plx0);

status = 0;
if(*plx0 < 0.001) status = -1;

return(status);
} 
/************************************************************************
* Scan the input table and make the modifications
*
* INPUT:
* fp_in: pointer to the input file containing the input table
* fp_out: pointer to the output Latex file
*
*************************************************************************/
static int latex_add_magnitude(FILE *fp_in, FILE *fp_out,
                               int in_mag_col, int in_absp_col,
                               int in_plx1_col, int in_plx2_col,
                               int in_err_plx1_col, int in_err_plx2_col,
                               int out_abspmag_col)
{
char in_line[MAX_LENGTH], out_line[MAX_LENGTH], buffer[MAX_LENGTH];
char plx_string[64], err_plx_string[64];
char *pc, absmag_string[64];
double absp0, plx0, err_plx0;
char target_str[128];
double Vmag, Vabsmag, M_odot, ww;
int iline, status, verbose_if_error = 0, absmag_strlen;
int i, in_line_length, nval, ncols_max;

// Maximum number of columns to be created
ncols_max = out_abspmag_col;

iline = -1;
while(!feof(fp_in)) {
  if(fgets(in_line, MAX_LENGTH, fp_in)) {
// Copy to out_line here in case it is not a digit:
  strcpy(out_line, in_line);
// Remove the end of line '\n' from input line:
  jlp_cleanup_string(out_line, MAX_LENGTH);

  if(isdigit(in_line[0]) != 0) {
    iline++;
// Remove the end of line '\n' from input line:
//    jlp_cleanup_string(in_line, MAX_LENGTH);

// Get the name in the 1st column 
   status = latex_get_column_item(in_line, buffer, 1, verbose_if_error);

// Remove dollar in input string
   latex_remove_dollar(buffer, target_str, 64);

   in_line_length = MAX_LENGTH;

   status = latex_read_parallax(in_line, in_absp_col, in_plx1_col, in_plx2_col, 
                                in_err_plx1_col, in_err_plx2_col, 
                                &plx0, &err_plx0, &absp0);

// If the apparent magnitude or the parallax could not be found, write nodata to those fields:
   sprintf(absmag_string, "\\nodata");

// Look for the parallax of the target:
   sprintf(plx_string, "\\nodata ");
   sprintf(err_plx_string, "\\nodata ");
   Vmag = 0.;
// Read the apparent magnitude from the good column
   status = latex_get_column_item(in_line, buffer, in_mag_col, 
                                  verbose_if_error);
   nval = sscanf(buffer, "%lf\n", &ww);
   if(nval == 1) {
      Vmag = ww;
      }
   if((plx0 > 0.01) && (Vmag > 0.01)) {
      sprintf(plx_string, "%.3f ", plx0);
      sprintf(err_plx_string, "%.3f ", err_plx0);
// Compute absolute magnitude Vabsmag, 
// with Vmag value and the parallax and absportion at that index:
      status = compute_abs_mag(Vmag, absp0, plx0, &Vabsmag);
      sprintf(absmag_string, "%.2f ", Vabsmag);
      }
printf("QQQQQQQQQQQ: plx0=%f Vmag=%f absmag_string=%s\n", plx0, Vmag, absmag_string);
// Copy in_line to out_line:
//   strcpy(out_line, in_line);
   latex_add_emptycols_to_ncols(in_line, out_line, MAX_LENGTH, ncols_max);

// Copy absmag to output file in out_absmag_col, i.e., the good column:
    absmag_strlen = strlen(absmag_string);
    verbose_if_error = 2;
    status = latex_set_column_item(out_line, in_line_length, absmag_string,
                                   absmag_strlen, out_abspmag_col, 
                                   verbose_if_error); 
   } // isdigit

// Save to output file:
   fprintf(fp_out, "%s\n", out_line);

  } /* EOF if fgets */
 } /* EOF while ... */
printf("latex_add_magnitude: %d lines sucessfully read and processed\n",
        iline);
return(0);
}
/**********************************************************************
* Compute absolute magnitude
*
* Absolute magnitude: apparent magnitude the the object would have
* if it was viewed from a distance of 10 parsecs without extinction
* The sun has an absolute magnitude M_V= +4.83
*
* Absolute magnitude M:
*      m - M = 5 log10(D) - 5
* D distance in parsecs : D = 1/pi
* Apparent magnitude: m = m_obs - A
* Interstellar absorption in the G band: A in magnitudes
**********************************************************************/
static int compute_abs_mag(double gmag, double absp, double plx, double *Gmag)
{
int status = -1;
 *Gmag = 0.;
 if(plx > 0.) {
// absp is positive, so I add it to gmag:
   *Gmag = gmag + absp - 5. * log10(1000. / plx) + 5.;
   status = 0;
   }
return(status);
}

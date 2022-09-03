/*************************************************************************
* Program latex_add_mass
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
static int latex_add_mass(FILE *fp_in, FILE *fp_out, int in_spect_col, 
                          int in_spec_part, int in_absmag_col, 
                          int out_mass_col);
static int decode_composite_spectral_type(char *buffer, int in_spect_part, 
                                          char *spect_string);
static int compute_abs_mag(double gmag, double absp, double plx, double *Gmag);
static int compute_mass_from_mag(double MVmag, char *spect_type, double *mass);

int main(int argc, char *argv[])
{
int i, status, nval;
int in_spect_col, in_absmag_col, out_mass_col; 
int in_spect_part;
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
  printf(" Syntax: latex_add_mass in_latex_file out_latex_file in_spect_col,in_spect_part,in_absmag_col,out_mass_col\n");
  exit(-1);
  }
else
  {
  strcpy(filein,argv[1]);
  strcpy(fileout,argv[2]);
  nval = sscanf(argv[3], "%d,%d,%d,%d", &in_spect_col, &in_spect_part,
                &in_absmag_col, &out_mass_col);
  if(nval != 4) {
    fprintf(stderr, "Error, should have 4 column numbers !\n");
    exit(-1);
    }
  }

printf(" OK: filein=%s fileout=%s \n", filein, fileout);
printf("in_spect_col=%d in_spect_part=%d \n", in_spect_col, in_spect_part);
printf("in_absmag_col=%d out_mass_col=%d \n", in_absmag_col, out_mass_col);

if((fp_in = fopen(filein,"r")) == NULL) {
  printf(" Fatal error opening input file %s \n",filein);
  exit(-1);
  }

if((fp_out = fopen(fileout,"w")) == NULL) {
  printf(" Fatal error opening output file %s \n",fileout);
  fclose(fp_in);
  exit(-1);
  }

latex_add_mass(fp_in, fp_out, in_spect_col, in_spect_part,  in_absmag_col, 
               out_mass_col);

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
   if(plx1_value > 0.) {
     *plx0 = plx1_value;
     }
   if(plx2_value > 0.) {
     *plx0 = plx2_value;
     }
   if((plx1_value > 0.) && (plx2_value > 0.)) {
     *plx0 = (plx1_value + plx2_value)/2.;
     }
   if(err_plx1_value > 0.) {
     *err_plx0 = err_plx1_value;
     }
   if(err_plx2_value > 0.) {
     *err_plx0 = err_plx2_value;
     }
   if((err_plx1_value > 0.) && (err_plx2_value > 0.)) {
     *err_plx0 = 1. / (1. / err_plx1_value + 1. / err_plx2_value);
     }
   if((plx1_value > 0.) && (err_plx1_value > 0.) && (plx2_value > 0.)
     && (err_plx2_value > 0.)) {
     *plx0 = (plx1_value / err_plx1_value + plx2_value / err_plx2_value) 
              / (1. / err_plx1_value + 1. / err_plx2_value);
     }
#ifdef DEBUG
   printf("OK: in_line=%s plx1_value=%f plx2_value=%f \n", 
          in_line, plx1_value, plx2_value);
   printf("OK: in_line=%s absp0=%f plx0=%f err_plx0=%f \n", 
          in_line, *absp0, *plx0, *err_plx0);
#endif

return(0);
} 
/************************************************************************
* Scan the input table and make the modifications
*
* INPUT:
* fp_in: pointer to the input file containing the input table
* fp_out: pointer to the output Latex file
*
*************************************************************************/
static int latex_add_mass(FILE *fp_in, FILE *fp_out, int in_spect_col, 
                          int in_spect_part, int in_absmag_col, 
                          int out_mass_col)
{
char in_line[MAX_LENGTH], out_line[MAX_LENGTH], buffer[MAX_LENGTH];
char *pc, absmag_string[64], spect_string[64], mass_string[64];
char target_str[128];
double Vabsmag, M_odot;
int iline, status, verbose_if_error = 0, absmag_strlen, mass_strlen;
int i, in_line_length, nval, ncols_max;

// Maximum number of columns to be created
ncols_max = out_mass_col;

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

// Read the spectral type from in_spect_col column
   status = latex_get_column_item(in_line, buffer, in_spect_col, 
                                  verbose_if_error);
   decode_composite_spectral_type(buffer, in_spect_part, spect_string);
 
// Read the absolute magnitude from in_absmag_col column
   status = latex_get_column_item(in_line, buffer, in_absmag_col, 
                                  verbose_if_error);
   if(status == 0) {
     nval = sscanf(buffer, "%lf\n", &Vabsmag);
     if(nval == 1) {
// Compute the mass of the star from its absolute magnitude
// compute_mass_from_mag(double MVmag, char *spect_type, double *mass);
       compute_mass_from_mag(Vabsmag, spect_string, &M_odot);
       sprintf(mass_string, "%.2f ", M_odot);
       } else {
// If the absolute magnitude could not be found, write nodata to the mass field:
       sprintf(mass_string, "\\nodata");
       }
     }
// Copy in_line to out_line:
   strcpy(out_line, in_line);
   latex_add_emptycols_to_ncols(in_line, out_line, MAX_LENGTH, ncols_max);
// Copy mass to output file in out_mass_col, i.e., the good column:
    mass_strlen = strlen(mass_string);
    status = latex_set_column_item(out_line, in_line_length, mass_string,
                                   mass_strlen, out_mass_col, 
                                   verbose_if_error); 
   } // isdigit

// Save to output file:
   fprintf(fp_out, "%s\n", out_line);

  } /* EOF if fgets */
 } /* EOF while ... */
printf("latex_add_mass: %d lines sucessfully read and processed\n",
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
/**********************************************************************
* Compute mass from absolute magnitude
*
* The sun has an absolute magnitude M_V= +4.83
*
* Absolute magnitude M:     m - M = 5 log10(D) - 5
% Bolometric magnitude: magnitude of the star based upon its total radidation
%   in all wavelengths
% Bolometric correction:  BC = M_{bol} - M_V
%   The bolometric correction is the correction made to the abolute magnitude of
%   an object in order to convert its visible magnitude to its bolometric
%   magnitude
%   It is large for stars which radiate most of their energy outside
%   the visible range
* BC of main sequence stars:
* O3 : -4.3
* G0 : -0.10
* G5 : -0.14
* K0 : -0.24
* K5 : -0.66
* M0 : -1.21
* 
* /home/text/tex/papers/scard07b_Orb/scard07b_Orb.pdf
* As proposed in Scardia et al.~(2008) have used the mass-luminosity relation:
* log_{10} Mass = -k (M_bol - M_0)
* with k = 0.1045 and M_0 = 4.60
* and where Mass is in solar mass unit
**********************************************************************/
static int compute_mass_from_mag(double MVmag, char *spect_type0, double *mass)
{
double BolomCorr, MBolom;
char sptype0, spect_type[64];
int ic, nval, sptype1, status = -1;

strcpy(spect_type, spect_type0);
jlp_compact_string(spect_type, 40);
printf("AAAA: spect_type:%s\n", spect_type);
 sptype0 = spect_type[0];
 ic = int(spect_type[1]) - int('0');
 sptype1 = 0;
 if(ic >=0 && ic <= 9) sptype1 = ic;
printf("AAAA: sptype0:%c sptype1=%d\n", sptype0, sptype1);
if(sptype0 == 'G') {
  BolomCorr = -0.10 - sptype1 * 0.04/5.;
  } else if (sptype0 == 'K') {
  BolomCorr = -0.24 - sptype1 * 0.42/5.;
  } else if (sptype0 == 'M') {
  BolomCorr = -1.21 - sptype1 * 0.42/5.;
  } else {
  BolomCorr = -0.66;
  }

MBolom = MVmag - BolomCorr;
*mass = pow(10., -0.1045 * (MBolom - 4.60));

printf("AAAA: Vmag=%f BC=%f Bo=%f mass=%f\n", MVmag, BolomCorr, MBolom, *mass);
return(status);
}
/**********************************************************************
* Decode composite spectral type into single spectral_type
*
* INPUT:
*   buffer: spectral type (example: "G5 + M3" or "K2"
*   in_spect_part: 1 if first part of input buffer or 2 if second part of buffer
* OUTPUT:
*   spect_string: spectral type of good compoonent (example "G5", M3" or "K2")
**********************************************************************/
static int decode_composite_spectral_type(char *buffer, int in_spect_part, 
                                          char *spect_string)
{
int status = -1;
char *pc;

// printf("buffer=%s\n", buffer);
strcpy(spect_string, buffer);
pc = spect_string;
while(*pc && *pc != '+') pc++;
if(in_spect_part == 1) {
    *pc = '\0';
  } else if(*pc == '+') {
  pc++;
  strcpy(spect_string, pc);
  }
// printf("in_spect_part=%d spect_string=%s\n", in_spect_part, spect_string);
return(status);
}

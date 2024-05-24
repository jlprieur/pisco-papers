/*************************************************************************
* Program cvt_bessel_julian
* 
* To process and convert Bessel to Julian epochs in Latex table
*
\hline 
00014+3937 & HLD60 & K0V+K1V & 9.09 & 9.77 & \nodata & 20.42 & 1.91 & 19.336 & 0.020 & 5.16 & & & \\
*
* JLP
* Version 10/08/2023
*************************************************************************/
#include <stdio.h>
#include <stdlib.h> // exit(-1)
#include <string.h>
#include <math.h>  // log10()
#include <ctype.h>  // isalpha(), isdigit()
#include "latex_utils.h" // latex_get_column_item(), latex_remove_column()
#include "jlp_string.h"  // (in jlplib/jlp_fits/ ) jlp_cleanup_..

#include "latex_utils.h"  // latex_get_column_item()
#include "jlp_fitsio.h" // BesselToJulian
//#include "astrom_utils1.h" 
//#include "astrom_utils2.h" 

#define MAX_LENGTH 256 

#define DEBUG

static int cvt_bessel_julian(FILE *fp_in, FILE *fp_out, int epoch_col);

int main(int argc, char *argv[])
{
int i, status, nval;
int epoch_col; 
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
  printf(" Syntax: cvt_bessel_julian in_latex_file out_latex_file epoch_col\n");
  exit(-1);
  }
else
  {
  strcpy(filein,argv[1]);
  strcpy(fileout,argv[2]);
  nval = sscanf(argv[3], "%d", &epoch_col);
  if(nval != 1) {
    fprintf(stderr, "Error, should have one column number !\n");
    exit(-1);
    }
  }

printf(" OK: filein=%s fileout=%s \n", filein, fileout);
printf("epoch_col=%d \n", epoch_col);

if((fp_in = fopen(filein,"r")) == NULL) {
  printf(" Fatal error opening input file %s \n",filein);
  exit(-1);
  }

if((fp_out = fopen(fileout,"w")) == NULL) {
  printf(" Fatal error opening output file %s \n",fileout);
  fclose(fp_in);
  exit(-1);
  }

cvt_bessel_julian(fp_in, fp_out, epoch_col);

fclose(fp_in);
fclose(fp_out);
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
static int cvt_bessel_julian(FILE *fp_in, FILE *fp_out, int epoch_col) 
{
char in_line[MAX_LENGTH], out_line[MAX_LENGTH], buffer[MAX_LENGTH];
char *pc, julian_string[64];
int iline, status, verbose_if_error = 0;
int i, in_line_length, nval, julian_strlen;
double BesselEpoch, JulianEpoch;

iline = -1;
while(!feof(fp_in)) {
  if(fgets(in_line, MAX_LENGTH, fp_in)) {
// Remove the end of line '\n' from input line:
  jlp_cleanup_string(in_line, MAX_LENGTH);
// Copy to out_line here in case it is not a digit:
  strcpy(out_line, in_line);

  if(isdigit(in_line[0]) != 0) {
    iline++;

// Get the name in the 1st column 
   status = latex_get_column_item(in_line, buffer, 1, verbose_if_error);

   in_line_length = MAX_LENGTH;

// Read the epoch from epoch_col column
   status = latex_get_column_item(in_line, buffer, epoch_col, 
                                  verbose_if_error);
   if(status == 0) {
     nval = sscanf(buffer, "%lf\n", &BesselEpoch);
     if(nval == 1) {
//       JulianEpoch = BesselEpoch;
// JLP_besselian_to_julian_epoch(double b_date, double *j_date);
       JLP_besselian_to_julian_epoch(BesselEpoch, &JulianEpoch);

       printf("epoch : bessel=%.3f julian=%.3f\n", BesselEpoch, JulianEpoch);
       sprintf(julian_string, "%.3f ", JulianEpoch);
// Copy mass to output file in out_mass_col, i.e., the good column:
       julian_strlen = strlen(julian_string);
       status = latex_set_column_item(out_line, in_line_length, julian_string,
                                      julian_strlen, epoch_col, 
                                      verbose_if_error); 
      }
    }
   } // isdigit

// Save to output file:
   fprintf(fp_out, "%s\n", out_line);

  } /* EOF if fgets */
 } /* EOF while ... */
printf("cvt_bessel_julian: %d lines sucessfully read and processed\n",
        iline);
return(0);
}

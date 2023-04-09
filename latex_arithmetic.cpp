/*************************************************************************
* Program latex_arithmetic
* 
* To process and update Latex table
*
\hline 
00014+3937 & HLD60 & K0V+K1V & 9.09 & 9.77 & \nodata & 20.42 & 1.91 & 19.336 & 0.020 & 5.16 & & & \\
*
* JLP
* Version 22/02/2023
*************************************************************************/
#include <stdio.h>
#include <stdlib.h> // exit(-1)
#include <string.h>
#include <math.h>  // log10()
#include <ctype.h>  // isalpha(), isdigit()
#include "latex_utils.h" // latex_get_column_item(), latex_remove_column()
#include "jlp_string.h"  // (in jlplib/jlp_fits/ ) jlp_cleanup_..

//#include "astrom_utils1.h" 
//#include "astrom_utils2.h" 

#define MAX_LENGTH 256 

#define DEBUG

static int latex_arithm(FILE *fp_in, FILE *fp_out, int in_col, 
                        double val_arithm);

int main(int argc, char *argv[])
{
int i, status, nval, in_col;
double val_arithm;
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

if(argc != 5)
  {
  printf(" Syntax: latex_add_mass in_latex_file out_latex_file in_col val_arithm\n");
  exit(-1);
  }
else
  {
  strcpy(filein,argv[1]);
  strcpy(fileout,argv[2]);
  sscanf(argv[3], "%d", &in_col);
  sscanf(argv[4], "%lf", &val_arithm);
  }

printf(" OK: filein=%s fileout=%s \n", filein, fileout);
printf("in_col=%d \n", in_col);
printf("val_arithm=%.6f \n", val_arithm);

if((fp_in = fopen(filein,"r")) == NULL) {
  printf(" Fatal error opening input file %s \n",filein);
  exit(-1);
  }

if((fp_out = fopen(fileout,"w")) == NULL) {
  printf(" Fatal error opening output file %s \n",fileout);
  fclose(fp_in);
  exit(-1);
  }

latex_arithm(fp_in, fp_out, in_col, val_arithm);

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
static int latex_arithm(FILE *fp_in, FILE *fp_out, int in_col, 
                        double val_arithm)
{
char in_line[MAX_LENGTH], out_line[MAX_LENGTH], buffer[MAX_LENGTH];
char *pc, buffer1[64];
char out_string[64];
int iline, status, verbose_if_error = 0;
int i, in_line_length, nval, out_col, out_strlen;
double in_value, out_value;

in_line_length = MAX_LENGTH;
out_col = in_col;
iline = -1;
while(!feof(fp_in)) {
  if(fgets(in_line, MAX_LENGTH, fp_in)) {
// Copy to out_line here in case it is not a digit:
     strcpy(out_line, in_line);
     fprintf(fp_out, "%s", in_line);
// Remove the end of line '\n' from input line:
     jlp_cleanup_string(out_line, MAX_LENGTH);

     out_value = 12345.;
//old  if(isdigit(in_line[0]) != 0) {
      iline++;

   status = latex_get_column_item(out_line, buffer, in_col, 
                                  verbose_if_error);
   if(status == 0) {
// Remove dollar in input string
     latex_remove_dollar(buffer, buffer1, 64);
     nval = sscanf(buffer1, "%lf", &in_value);
     if(nval == 1) {
       out_value = in_value * val_arithm;
       sprintf(out_string, "%.3f ", out_value);
// Copy mass to output file in out_mass_col, i.e., the good column:
       out_strlen = strlen(out_string);
       status = latex_set_column_item(out_line, in_line_length, out_string,
                                      out_strlen, out_col, verbose_if_error); 
// Save to output file:
       if(status == 0) fprintf(fp_out, "%s\n", out_line);
       }
     }
//old   } // isdigit
  } /* EOF if fgets */
 } /* EOF while ... */
printf("latex_add_mass: %d lines sucessfully read and processed\n",
        iline);
return(0);
}

/*************************************************************************
* Program latex_to_ascii
* To convert a LaTeX table to a csv file 
*
* JLP
* Version 02/09/2023
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>   /* exit() */
#include <math.h>
#include <string.h>
#include <ctype.h>  // isalpha(), isdigit()
#include "latex_utils.h" // jlp: latex_read_fvalue...

#define DEBUG
/*
*/

static int jlp_latex_to_csv(char *filein, char *fileout);

int main(int argc, char *argv[])
{
char filein[60], fileout[60];
int i;

  printf("latex_to_csv/ JLP/ Version 02/09/2023\n");
  printf("Note that this program can handle multiple LaTeX tables\n\n");

if(argc == 7 && *argv[4]) argc = 5;
if(argc == 7 && *argv[3]) argc = 4;
if(argc == 7 && *argv[2]) argc = 3;
if(argc == 7 && *argv[1]) argc = 2;
if(argc != 3)
  {
  printf("Error: argc=%d\n\n", argc);
  printf("Syntax: latex_to_ascii in_latex_table out_csv_file \n");
  exit(-1);
  }
else
  {
  strcpy(filein,argv[1]);
  strcpy(fileout,argv[2]);
  }

printf(" OK: filein=%s fileout=%s\n", filein, fileout);

/* Scan the file and make the conversion: */
jlp_latex_to_csv(filein, fileout);

return(0);
}
/*************************************************************************
* Input file: plain ascii table (possibly without header)
*
*************************************************************************/
#define NMAX 1024
static int jlp_latex_to_csv_old(char *filein, char *fileout) 
{
char in_line[NMAX], out_line[NMAX];
int status, iline, i, icol, quote_is_on;
char *pc, *pc1;
FILE *fp_in, *fp_out;

if((fp_in = fopen(filein,"r")) == NULL)
{
printf(" Fatal error opening input file %s \n",filein);
exit(-1);
}

if((fp_out = fopen(fileout,"w")) == NULL)
{
printf(" Fatal error opening output file %s \n",fileout);
fclose(fp_in);
exit(-1);
}
fprintf(fp_out,"%s\n", "WDS, Name, Epoch, Bin., $\\rho$, $\\sigma_\\rho$, $\\theta$, $\\sigma_\\theta$, $\\Delta$m, Notes, Orbit, $\\Delta \\rho$(O-C), $\\Delta \\theta$(O-C), Grade"); 
fprintf(fp_out,"%s\n", ", , , , (\\arcsec) , (\\arcsec) ,  ($^\\circ$) , ($^\\circ$) , , ,  , (\\arcsec) , ($^\\circ$)"); 

iline = -1;
quote_is_on = 0;
while(!feof(fp_in))
{
  if(fgets(in_line, NMAX, fp_in))
  {
  if(isdigit(in_line[0])) {
    iline++;
// Remove the "\\" and the "\n" at the end of the line
    pc = in_line;
    while(*pc) { 
       if((*pc == '\\') && (*(pc+1) == '\\')) break;
       if(*pc == '\n') break;
       pc++;
     }
    *pc = '\0';

#ifdef DEBUG
printf(" iline=%d Data in line: >%s<\n", iline, in_line);
#endif
// Remove the "&" and replace it by "," 
      pc = in_line;
      pc1 = out_line;
      *pc1 = '\"'; pc1++;
      quote_is_on = 1;
      icol = 1;
      while(*pc) 
         {
         if((*pc == '$') && (*(pc+1) == '-') && (*(pc+2) == '$')) {
           *pc1 = '-'; pc1++;
           pc++; pc++; pc++;
           } 
// Handle general case and avoid "\&":
         else if((*pc == '&') && (*(pc-1) != '\\')) 
           { 
// icol=1 "00321+3538"
// icol=2 "HO2Aa,Ab"
             if(icol == 1) {
                *pc1 = '\"'; pc1++;
                 quote_is_on = 0;
                *pc1 = ','; pc1++; 
                *pc1 = '\"'; pc1++;
                 quote_is_on = 1;
               } else if (icol == 2) {
                *pc1 = '\"'; pc1++;
                 quote_is_on = 0;
                *pc1 = ','; pc1++; 
               } else if (icol == 9) {
                *pc1 = ','; pc1++; 
                *pc1 = '\"'; pc1++;
                 quote_is_on = 1;
               } else if (icol == 10) {
                *pc1 = '\"'; pc1++;
                 quote_is_on = 0;
                *pc1 = ','; pc1++; 
                *pc1 = '\"'; pc1++;
                 quote_is_on = 1;
               } else if (icol == 11) {
                *pc1 = '\"'; pc1++;
                 quote_is_on = 0;
                *pc1 = ','; pc1++; 
               } else {
                *pc1 = ','; pc1++; 
               }
             icol++;
             pc++; // skip '&'
           } else {
             *pc1 = *pc;
             pc++;
             pc1++;
           }
         }
// Case of short line:
      if(quote_is_on == 1) {
         *pc1 = '\"'; pc1++;
         quote_is_on = 0;
        }
      *pc1 = '\0';
#ifdef DEBUG
printf(" iline=%d Data out line: >%s<\n", iline, out_line);
#endif
     fprintf(fp_out,"%s\n", out_line); 
   } /* EOF successful reading of buffer from the input file */
  } // if not %
} /* EOF while !feof(fp_in) */

fclose(fp_in);
fclose(fp_out);
return(0);
}
/***************************************************
* New version
****************************************************/
static int jlp_latex_to_csv(char *filein, char *fileout) 
{
char in_line[NMAX], out_line[NMAX];
int status, iline, i, icol, quote_is_on;
char *pc, *pc1;
FILE *fp_in, *fp_out;

if((fp_in = fopen(filein,"r")) == NULL)
{
printf(" Fatal error opening input file %s \n",filein);
exit(-1);
}

if((fp_out = fopen(fileout,"w")) == NULL)
{
printf(" Fatal error opening output file %s \n",fileout);
fclose(fp_in);
exit(-1);
}
fprintf(fp_out,"%s\n", "WDS, Name, Epoch, Bin., $\\rho$, $\\sigma_\\rho$, $\\theta$, $\\sigma_\\theta$, $\\Delta$m, Notes, Orbit, $\\Delta \\rho$(O-C), $\\Delta \\theta$(O-C), Grade"); 
fprintf(fp_out,"%s\n", ", , , , (\\arcsec) , (\\arcsec) ,  ($^\\circ$) , ($^\\circ$) , , ,  , (\\arcsec) , ($^\\circ$)"); 

iline = -1;
quote_is_on = 0;
while(!feof(fp_in))
{
  if(fgets(in_line, NMAX, fp_in))
  {
  if(isdigit(in_line[0])) {
    iline++;
// Remove the "\\" and the "\n" at the end of the line
    pc = in_line;
    while(*pc) { 
       if((*pc == '\\') && (*(pc+1) == '\\')) break;
       if(*pc == '\n') break;
       pc++;
     }
    *pc = '\0';

#ifdef DEBUG
printf(" iline=%d Data in line: >%s<\n", iline, in_line);
#endif
// Insert the first \" 
      pc = in_line;
      pc1 = out_line;
      *pc1 = '\"'; pc1++;
      quote_is_on = 1;
      icol = 1;
      while(*pc) 
         {
         if((*pc == '$') && (*(pc+1) == '-') && (*(pc+2) == '$')) {
           *pc1 = '-'; pc1++;
           pc++; pc++; pc++;
           } 
// Handle general case and avoid "\&":
         else if((*pc == '&') && (*(pc-1) != '\\')) 
           { 
             if(quote_is_on == 1) {
                *pc1 = '\"'; pc1++;
                 quote_is_on = 0;
                *pc1 = ','; pc1++; 
                *pc1 = '\"'; pc1++;
                 quote_is_on = 1;
               }
             icol++;
             pc++; // skip '&'
           } else {
             *pc1 = *pc;
             pc++;
             pc1++;
           }
         }
// Case of short line:
      if(quote_is_on == 1) {
         *pc1 = '\"'; pc1++;
         quote_is_on = 0;
        }
      *pc1 = '\0';
#ifdef DEBUG
printf(" iline=%d Data out line: >%s<\n", iline, out_line);
#endif
     fprintf(fp_out,"%s\n", out_line); 
   } /* EOF successful reading of buffer from the input file */
  } // if not %
} /* EOF while !feof(fp_in) */

fclose(fp_in);
fclose(fp_out);
return(0);
}

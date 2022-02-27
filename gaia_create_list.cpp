/*************************************************************************
* Program gaia_create_list 
* 
* To create a list with WDSxxx from a Latex array with measurements 
*
\hline 
00022+2705 & 000210.18+270455.6 & BU733AB & G5Vb+K5V & 5.83 & 8.9 & & & 2.00 \\
*
* JLP
* Version 22/01/2022
*************************************************************************/
#include <stdio.h>
#include <stdlib.h> // exit(-1)
#include <string.h>
#include <ctype.h>  // isalpha(), isdigit()
#include "latex_utils.h" // latex_get_column_item(), latex_remove_column()
#include "jlp_string.h"  // (in jlplib/jlp_fits/ ) jlp_cleanup_..

#include "latex_utils.h"  // latex_get_column_item()
//#include "astrom_utils1.h" 
//#include "astrom_utils2.h" 

static int create_gaia_list(FILE *fp_in, FILE *fp_out);

int main(int argc, char *argv[])
{
int i, status; 
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

if(argc != 3)
  {
  printf(" Syntax: gaia_create_list in_latex_file out_txt_file \n");
  exit(-1);
  }
else
  {
  strcpy(filein,argv[1]);
  strcpy(fileout,argv[2]);
  }

printf(" OK: filein=%s fileout=%s \n", filein, fileout);

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

create_gaia_list(fp_in, fp_out);

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
static int create_gaia_list(FILE *fp_in, FILE *fp_out)
{
char in_line[256], in_line3[256];
char buffer[256], *pc, wds_name[256];
int iposition, iline, status, ival, verbose_if_error = 0;
int i;

iline = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 256, fp_in)) {
    iline++;
// Remove the end of line '\n' from input line:
    jlp_cleanup_string(in_line, 256);

  if(isdigit(in_line[0]) != 0) {
// Get the name in the 1st column 
    status = latex_get_column_item(in_line, buffer, 1, verbose_if_error);

// Remove $ in the name:
   pc = buffer;
   i = 0;
   while(*pc ) {
     if(*pc != '$') {
       wds_name[i++] = *pc; 
       }
     pc++;
     }
   wds_name[i] = '\0';
 
// Save to output file:
    if(status == 0) fprintf(fp_out, "WDS%s\n", wds_name);
   }

  } /* EOF if fgets */
 } /* EOF while ... */
printf("create_gaia_list: %d lines sucessfully read and processed\n",
        iline);
return(0);
}

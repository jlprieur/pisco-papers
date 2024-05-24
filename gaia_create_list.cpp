/*************************************************************************
* Program gaia_create_list 
* 
* To create a list with WDSxxx from a Latex array with measurements 
*
\hline 
00022+2705 & 000210.18+270455.6 & BU733AB & G5Vb+K5V & 5.83 & 8.9 & & & 2.00 \\
*
* JLP
* Version 14/08/2023
*************************************************************************/
#include <stdio.h>
#include <stdlib.h> // exit(-1)
#include <string.h>
#include <math.h>
#include <ctype.h>  // isalpha(), isdigit()
#include "latex_utils.h" // latex_get_column_item(), latex_remove_column()
#include "WDS_catalog_utils.h"       /* Routines used to read WDS catalog */
#include "jlp_string.h"  // (in jlplib/jlp_fits/ ) jlp_cleanup_..
#define PI 3.14159

#include "latex_utils.h"  // latex_get_column_item()
//#include "astrom_utils1.h" 
//#include "astrom_utils2.h" 

static int create_gaia_list1(FILE *fp_in, FILE *fp_out);
static int create_gaia_list2(FILE *fp_in, FILE *fp_out, char *WDS_catalog);

int main(int argc, char *argv[])
{
int i, status, nval, iopt; 
char filein[128], fileout[128], WDS_catalog[128];
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
  printf(" Syntax: gaia_create_list in_latex_file WDS_catalog out_txt_file iopt\n");
  printf(" iopt=1 add 'WDS'    iopt=2 accurate coordinates\n");
  exit(-1);
  }
else
  {
  strcpy(filein, argv[1]);
  strcpy(WDS_catalog, argv[2]);
  strcpy(fileout, argv[3]);
  nval = sscanf(argv[4], "%d", &iopt);
  }

printf(" OK: filein=%s WDS_cat = %s fileout=%s iopt=%d\n", filein, WDS_catalog, fileout, iopt);

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

if(iopt == 1)
 create_gaia_list1(fp_in, fp_out);
else
 create_gaia_list2(fp_in, fp_out, WDS_catalog);

fclose(fp_in);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input table and make the modifications (add 'WDS')
*
* INPUT:
* fp_in: pointer to the input file containing the input table
* fp_out: pointer to the output Latex file
*
*************************************************************************/
static int create_gaia_list1(FILE *fp_in, FILE *fp_out)
{
char in_line[256];
char buffer[256], *pc, wds_name[256];
int iline, status, verbose_if_error = 0;
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
/************************************************************************
* Scan the input table and make the modifications (accurate coordinates)
*
* INPUT:
* fp_in: pointer to the input file containing the input table
* fp_out: pointer to the output Latex file
*
*************************************************************************/
static int create_gaia_list2(FILE *fp_in, FILE *fp_out, char *WDS_catalog)
{
char in_line[256], disc_name[64];
char buffer[256], *pc, wds_name[256], str_alpha[64], str_delta[64];
double alpha, delta, equinox, rho, theta, d_alpha, d_delta;
double alpha_scale, delta_scale, new_alpha, new_delta;
int iline, status, nval, ifound, verbose_if_error = 0;
int i, icol_wds, icol_name, icol_rho, icol_theta;

fprintf(fp_out, "str_alpha str_delta alpha(deg) delta(deg) new_alpha new_delta (equinox 2000)\n"); 

iline = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 256, fp_in)) {
    iline++;
// Remove the end of line '\n' from input line:
    jlp_cleanup_string(in_line, 256);

printf(" iline=%d\n %s\n", iline, in_line);
  if(isdigit(in_line[0]) != 0) {
// Get the wds name in the 1st column 
    icol_wds = 1;
    status = latex_get_column_item(in_line, buffer, icol_wds, verbose_if_error);

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

// Get the discover name in the 2nd column 
    icol_name = 2;
    status = latex_get_column_item(in_line, disc_name, icol_name, verbose_if_error);

// Get rho in the column 5 
    icol_rho = 5;
    status = latex_get_column_item(in_line, buffer, icol_rho, verbose_if_error);
    nval = sscanf(buffer, "%lf", &rho);
    printf("rho:%s rho=%.3f\n", buffer, rho);
 
// Get theta in the column 7 
    icol_theta = 7;
    status = latex_get_column_item(in_line, buffer, icol_theta, verbose_if_error);
    nval = sscanf(buffer, "%lf", &theta);
    printf("theta:%s theta=%.1f\n", buffer, theta);

    alpha_scale = (24. / 360.) / (60.* 60.);
    delta_scale =  1. / (60. * 60.);
    d_alpha = rho * sin(theta * PI/ 180.) * alpha_scale;
    d_delta =  - rho * cos(theta * PI/ 180.) * delta_scale;
 
/* Look for accurate coordinates in WDS catalog
* alpha, delta, equinox: coordinates of the object
*            (alpha in hours and delta in degrees)
*/
  read_coordinates_from_WDS_catalog(wds_name, WDS_catalog, str_alpha,
                                    str_delta, &alpha, &delta,
                                    &equinox, &ifound);
  new_alpha = alpha + d_alpha;
  new_delta = delta + d_delta;
// Save to output file:
    if(ifound != 0) 
        fprintf(fp_out, "%s %s %.6f %.6f %.6f %.6f %s %s (%.3f %.1f)\n", 
                str_alpha, str_delta, alpha*365/24., delta, new_alpha*365./24.,
	        new_delta, wds_name, disc_name, rho, theta);
     else 
        fprintf(fp_out, "error for: %s %s\n", wds_name, disc_name);
   }

  } /* EOF if fgets */
 } /* EOF while ... */
printf("create_gaia_list: %d lines sucessfully read and processed\n",
        iline);
return(0);
}

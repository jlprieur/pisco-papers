/************************************************************************
* "merge_calib_resid_gili.c"
*
* To merge the two Latex tables: calibrated table and residual table
*
* JLP 
* Version 25/02/2020
*************************************************************************/
#include "jlp_catalog_utils.h"
#include "jlp_string.h"
#include "latex_utils.h"

/*
#define DEBUG
#define DEBUG_1
*/

int add_residuals_to_calib(FILE *fp_calib, char *resid_fname, FILE *fp_out,
                           int gili_format); 
static int modify_LaTeX_header(char *in_line, FILE *fp_out, int gili_format);
static int latex_truncate_columns(char *in_line, char *buffer, int ncolumns);

int main(int argc, char *argv[])
{
char calib_fname[80], resid_fname[80], out_fname[80];
int gili_format;
FILE *fp_calib, *fp_out;
time_t t = time(NULL);

if(argc != 5) {
  printf("Syntax: merge_calib_resid_gili calibrated_table residual_table out_table gili_format\n");
  return(-1);
}
strcpy(calib_fname, argv[1]);
strcpy(resid_fname, argv[2]);
strcpy(out_fname, argv[3]);
sscanf(argv[4], "%d", &gili_format);

printf("OK: calib=%s resid=%s output=%s gili_format=%d\n", 
       calib_fname, resid_fname, out_fname, gili_format); 

/* Open input calibrated table: */
if((fp_calib = fopen(calib_fname, "r")) == NULL) {
   fprintf(stderr, "merge_calib_resid/Fatal error opening calib. table %s\n",
           calib_fname);
   return(-1);
  }

/* Open output table: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "merge_calib_resid/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the Latex table: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "%% Merged table from: %s and %s \n%% Created on %s", 
        calib_fname, resid_fname, ctime(&t));
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

/* Scan the input calibrated table and add the residuals 
*/
add_residuals_to_calib(fp_calib, resid_fname, fp_out, gili_format); 

/* Close opened files:
*/
fclose(fp_calib);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the calibrated table and add the residuals 
*
* INPUT:
* fp_calib: pointer to the input file containing the calibrated table 
* resid_fname: name of the file containing the Latex table with the residuals
* fp_out: pointer to the output Latex file with the full table 
*
*************************************************************************/
int add_residuals_to_calib(FILE *fp_calib, char *resid_fname, FILE *fp_out,
                           int gili_format) 
{
char in_line[256], wds_name[40]; 
char object_name[40], discov_name[40], comp_name[40]; 
char buffer[256], orbit_ref[60*50], quadrant_discrep[20]; 
char sign_Drho[20], sign_Dtheta[20];
double epoch_o, rho_o, rho_o_c[50], theta_o_c[50];
int iline, norbits_found, status, verbose_if_error = 0, length1;
int ref_slength = 60, nmax_orbits = 50, orbit_grade;
int irho, ncolumns;
register int i;

/* Read rho in column irho=5 if Gili's format: */
      if(gili_format == 1) {
        irho = 5;
      } else {
        irho = 6;
      }
iline = 0;
while(!feof(fp_calib)) {
  if(fgets(in_line, 256, fp_calib)) {
    iline++;
// Read header (assume 10 columns) and add 3 columns for output header
// In gili'format: ADS col. is replaced by Dm col., but no filter, no orbit... 
// Good lines start with a digit (WDS names...) or with \idem
// Lines starting with % are ignored
    if((in_line[0] != '%') && (isdigit(in_line[0])) ) {
// Remove the end of line '\n' from input line:
      jlp_cleanup_string(in_line, 256);

// New version (after 2011): now search for orbit in all cases 
// since the orbit flag may be wrong

/* Retrieve the residuals for this object and epoch in the residual table: */
/* Read object name and companion name from columns 2 and 3: */
      read_object_name_from_CALIB_line_gili(in_line, wds_name, discov_name, 
                                               comp_name);
/* Generate the object name: */
      strcpy(object_name, discov_name);
#ifdef DEBUG
      printf(" %s discov_name=%s \n", in_line, discov_name);
#endif

/* Read epoch in column 3: */
      status = latex_get_column_item(in_line, buffer, 3, verbose_if_error);
      if(status || (sscanf(buffer,"%lf", &epoch_o) != 1)) {
         fprintf(stderr,"Fatal error reading epoch_o: %s (line=%d)\n", 
                 in_line, iline);
         exit(-1);
         }

/* Read rho in column irho=5 if Gili's format: */
      rho_o = 0.; 
      status = latex_get_column_item(in_line, buffer, irho, verbose_if_error);
      if(status || (sscanf(buffer,"%lf", &rho_o) != 1)) {
         fprintf(stderr,"Warning: error reading rho_o: %s (line=%d) Unres ?\n", 
                 in_line, iline);
         }
#ifdef DEBUG
      printf(" object=%s comp=%s epoch=%.4f rho=%.4f \n", object_name, comp_name, epoch_o, rho_o);
#endif

/* Retrieve the residuals for this object and epoch in the residual table: */
      strcpy(quadrant_discrep,"");
      get_values_from_RESID_table(resid_fname, object_name, comp_name, epoch_o,
                                  rho_o, orbit_ref, &orbit_grade,
                                  ref_slength, rho_o_c, 
                                  theta_o_c, quadrant_discrep, &norbits_found, 
                                  nmax_orbits);
      if(!norbits_found) {
         fprintf(stderr," Residuals for %s %s not found in resid_table for epoch=%f and rho=%f in %s\n", 
                 object_name, comp_name, epoch_o, rho_o, resid_fname);
/* Simply copy the input line to the output file: */
         fprintf(fp_out, "%s \n", in_line);
         } else {
/* quadrant_discrep = "$^Q$" if Quadrant discrepancy between measure and orbit */
         for(i = 0; i < norbits_found; i++) {
           if(rho_o_c[i] < 0.) {
             rho_o_c[i] *= -1;
             strcpy(sign_Drho,"$-$");
             } else {
             strcpy(sign_Drho,"");
             }
           if(theta_o_c[i] < 0.) {
             theta_o_c[i] *= -1;
             strcpy(sign_Dtheta,"$-$");
             } else {
             strcpy(sign_Dtheta,"");
             }
// Truncate the line if too many columns:
        if(gili_format == 1)
           ncolumns = 10;
// JLP 2023: same behaviour now:
        else
           ncolumns = 10;
        status = latex_truncate_columns(in_line, buffer, ncolumns);

/*
           printf("EEZZ/ >%s<\n i=%d iline=%d length=%d -1=%c -2=%c\n", 
             in_line, i, buffer, length1, buffer[length1-1], buffer[length1-2]);
*/
           fprintf(fp_out, "%s & %s & %s%.2f & %s%.1f%s & %d \\\\\n", 
                   buffer, &orbit_ref[i*ref_slength], sign_Drho, rho_o_c[i], 
                   sign_Dtheta, theta_o_c[i], quadrant_discrep, orbit_grade);
           } // EOF for i loop
      } /* EOF for (i=0, norbits_found) */
    } else if (in_line[0] == '%') {
/* Simply copy the input line to the output file if it is a comment: */
/* Remove the end of line '\n' from input line: */
      jlp_cleanup_string(in_line, 256);
      fprintf(fp_out, "%s\n", in_line);
    } else {
      modify_LaTeX_header(in_line, fp_out, gili_format);
    }/* EOF if !isdigit ... */
  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("add_residuals_to_calib: %d lines sucessfully read and processed\n", 
        iline);
return(0);
}
/********************************************************************
* Modify the lines with Latex syntax
* to add 3 columns with the residuals (orbit reference, rho_O-C, theta_O-C)
*
*********************************************************************/
static int modify_LaTeX_header(char *in_line, FILE *fp_out, int gili_format)
{
char buffer[80];

/* Copy input line to "compacted" form in order to perform the tests safely */
strncpy(buffer, in_line,80);
buffer[79]='\0';
// Error if compact buffer !!!
//jlp_compact_string(buffer, 80);

/* begin{tabular*} */
 if(!strncmp(buffer, "\\begin{tabular*}", 16)) {
    fprintf(fp_out, "\\small \n");
    fprintf(fp_out, "\\begin{tabular*}{\\textwidth}{cllccccrrllrrc} \n");
/* Extended header */
 } else if(!strncmp(buffer, "WDS", 3)) {
if(gili_format == 1) {
    fprintf(fp_out,"WDS & Name & Epoch & Bin. & $\\rho$ \
& $\\sigma_\\rho$ & \\multicolumn{1}{c}{$\\theta$} & $\\sigma_\\theta$ \
& Dm & Notes & Orbit & {\\scriptsize $\\Delta \\rho$(O-C)} \
& {\\scriptsize $\\Delta \\theta$(O-C)} & Grade \\\\ \n");
} else {
    fprintf(fp_out,"WDS & Name & Epoch & Filt. & Eyep. & $\\rho$ \
& $\\sigma_\\rho$ & \\multicolumn{1}{c}{$\\theta$} & $\\sigma_\\theta$ \
& Notes & Orbit & {\\scriptsize $\\Delta \\rho$(O-C)} \
& {\\scriptsize $\\Delta \\theta$(O-C)} & Grade \\\\ \n");
}
/* JLPPPP/DDEBUG: To be tested */
 } else if(!strncmp(buffer, "& &     &     & (", 17)){
    fprintf(fp_out,"& & & & (\\arcsec) & (\\arcsec) &  \\multicolumn{1}{c}{($^\\circ$)} \
& ($^\\circ$) & & & (\\arcsec) & ($^\\circ$) \\\\ \n");
 } else {
/* For the other lines, simply copy the input line to the output file: */
/* Remove the end of line '\n' from input line: */
  jlp_cleanup_string(in_line, 256);
  fprintf(fp_out, "%s\n", in_line);
 }
return(0);
}
/************************************************************************
*
*
*************************************************************************/
static int latex_truncate_columns(char *in_line, char *buffer, int ncolumns)
{
int i;
char *pc;

// printf("DEBUG: input \n >%s< \n", in_line);
strcpy(buffer, in_line);
pc = buffer;
for(i = 0; i < ncolumns; i++) {
//  printf("DEBUG: column=%d \n", i);
  while(*pc && (*pc != '&') && (*pc != '\\')) pc++;
// Problem with \rlap...
  if(*pc == '\\') {
    if(*(pc+1) == '\\') {
     *pc = '\0'; 
      break;
     }
    }
  if(*pc) 
    pc++;
  else
    break;
  }
*pc = '\0';
// printf("DEBUG: output \n >%s< \n", buffer);
return(0);
}

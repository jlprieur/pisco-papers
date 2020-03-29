/************************************************************************
* "merge_calib_resid_gili.c"
*
* To merge the two Latex tables: calibrated table and residual table
*
* JLP 
* Version 25/02/2020
*************************************************************************/
#include "jlp_catalog_utils.h"

#define DEBUG
#define DEBUG_1
/*
*/

int add_residuals_to_calib(FILE *fp_calib, char *resid_fname, FILE *fp_out,
                           int gili_format); 
static int modify_LaTeX_header(char *in_line, FILE *fp_out);

int main(int argc, char *argv[])
{
char calib_fname[80], resid_fname[80], out_fname[80];
int gili_format;
FILE *fp_calib, *fp_out;
time_t t = time(NULL);


if(argc != 4) {
  printf("Syntax: merge_calib_resid_gili calibrated_table residual_table out_table\n");
  return(-1);
}
strcpy(calib_fname, argv[1]);
strcpy(resid_fname, argv[2]);
strcpy(out_fname, argv[3]);

printf("OK: calib=%s resid=%s output=%s \n", calib_fname, resid_fname, 
           out_fname); 

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
gili_format = 1;
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
char previous_object_name[40], previous_comp_name[40];
char sign_Drho[20], sign_Dtheta[20];
double epoch_o, rho_o, rho_o_c[50], theta_o_c[50];
int iline, norbits_found, status, verbose_if_error = 0, length1;
int ref_slength = 60, nmax_orbits = 50;
register int i;

previous_object_name[0] = '\0';
previous_comp_name[0] = '\0';
iline = 0;
while(!feof(fp_calib)) {
  if(fgets(in_line, 256, fp_calib)) {
    iline++;
// Read header (assume 10 columns) and add 3 columns for output header
// In gili'format: ADS col. is replaced by Dm col., but no filter, no orbit... 
    if(!strncmp(in_line,"& & & & & & & & & \\\\", 20)){
      fprintf(fp_out,   "& & & & & & & & & & & & \\\\ \n");
// Good lines start with a digit (WDS names...) or with \idem
// Lines starting with % are ignored
    } else if(in_line[0] != '%' && (isdigit(in_line[0]) 
        || !strncmp(in_line, "\\idem", 5))) {
// Remove the end of line '\n' from input line:
      cleanup_string(in_line, 256);

// New version (after 2011): now search for orbit in all cases 
// since the orbit flag may be wrong

/* Retrieve the residuals for this object and epoch in the residual table: */
/* Read object name and companion name from columns 2 and 3: */
      if(isdigit(in_line[0])) {
         read_object_name_from_CALIB_line_gili(in_line, wds_name, discov_name, 
                                               comp_name);
/* Generate the object name: */
          strcpy(object_name, discov_name);
      } else {
         strcpy(object_name, previous_object_name);
         strcpy(comp_name, previous_comp_name);
      }

/* Read epoch in column 3: */
      status = latex_get_column_item(in_line, buffer, 3, verbose_if_error);
      if(status || (sscanf(buffer,"%lf", &epoch_o) != 1)) {
         fprintf(stderr,"Fatal error reading epoch_o: %s (line=%d)\n", 
                 in_line, iline);
         exit(-1);
         }

/* Read rho in column 5: */
      rho_o = 0.; 
      status = latex_get_column_item(in_line, buffer, 5, verbose_if_error);
      if(status || (sscanf(buffer,"%lf", &rho_o) != 1)) {
         fprintf(stderr,"Warning: error reading rho_o: %s (line=%d) Unres ?\n", 
                 in_line, iline);
         }

/* Retrieve the residuals for this object and epoch in the residual table: */
      strcpy(quadrant_discrep,"");
      get_values_from_RESID_table(resid_fname, object_name, comp_name, epoch_o,
                                  rho_o, orbit_ref, ref_slength, rho_o_c, 
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
// Remove \\ from in_line and copy to buffer:
           strcpy(buffer, in_line);
           length1 = strlen(buffer);
           printf("EEZZ/ >%s<\n i=%d iline=%d length=%d -1=%c -2=%c\n", 
             in_line, i, buffer, length1, buffer[length1-1], buffer[length1-2]);
           buffer[length1-2] = '\0';
           fprintf(fp_out, "%s & %s & %s%.2f & %s%.1f%s \\\\\n", 
                   buffer, &orbit_ref[i*ref_slength], sign_Drho, rho_o_c[i], 
                   sign_Dtheta, theta_o_c[i], quadrant_discrep);
           } // EOF for i loop
      } /* EOF for (i=0, norbits_found) */
    } else if (in_line[0] == '%') {
/* Simply copy the input line to the output file if it is a comment: */
/* Remove the end of line '\n' from input line: */
      cleanup_string(in_line, 256);
      fprintf(fp_out, "%s\n", in_line);
    } else {
      modify_LaTeX_header(in_line, fp_out);
    }/* EOF if !isdigit ... */
  } /* EOF if fgets */ 
/* Load object name to handle the case of "idem" (i.e. multiple measurements
* of the same object, without repeating the object name in cols. 1 2 3)*/
 strcpy(previous_object_name, object_name);
 strcpy(previous_comp_name, comp_name);
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
static int modify_LaTeX_header(char *in_line, FILE *fp_out)
{
char buffer[80];

/* Copy input line to "compacted" form in order to perform the tests safely */
strncpy(buffer, in_line,80);
buffer[79]='\0';
compact_string(buffer, 80);

/* begin{tabular*} */
 if(!strncmp(buffer, "\\begin{tabular*}", 16)) {
    fprintf(fp_out, "\\small \n");
    fprintf(fp_out, "\\begin{tabular*}{\\textwidth}{clrcccccllllrr} \n");
/* Extended header */
 } else if(!strncmp(buffer, "WDS", 3)) {
    fprintf(fp_out,"WDS & Name & Epoch & Bin. & $\\rho$ \
& $\\sigma_\\rho$ & \\multicolumn{1}{c}{$\\theta$} & $\\sigma_\\theta$ \
& Dm & Notes & Orbit & {\\scriptsize $\\Delta \\rho$(O-C)} \
& {\\scriptsize $\\Delta \\theta$(O-C)} \\\\ \n");
/* JLPPPP/DDEBUG: To be tested */
 } else if(!strncmp(buffer, "& &     &     & (\")", 20)){
    fprintf(fp_out,"& & & & (\\arcsec) & (\\arcsec) &  \\multicolumn{1}{c}{($^\\circ$)} \
& ($^\\circ$) & & & (\\arcsec) & ($^\\circ$) \\\\ \n");
 } else {
/* For the other lines, simply copy the input line to the output file: */
/* Remove the end of line '\n' from input line: */
  cleanup_string(in_line, 256);
  fprintf(fp_out, "%s\n", in_line);
 }
return(0);
}

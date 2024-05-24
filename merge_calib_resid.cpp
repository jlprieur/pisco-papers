/************************************************************************
* "merge_calib_resid.c"
*
* To merge the two Latex tables: calibrated table and residual table
*
* JLP 
* Version 15/09/2011
*************************************************************************/
#include "jlp_catalog_utils.h"
#include "jlp_string.h"
#include "latex_utils.h"

#define DEBUG
#define DEBUG_1
/*
*/

int add_residuals_to_calib(FILE *fp_calib, char *resid_fname, FILE *fp_out); 
static int modify_LaTeX_header(char *in_line, FILE *fp_out);
static int remove_orbit_column(char *in_line, int ilen);

int main(int argc, char *argv[])
{
char calib_fname[80], resid_fname[80], out_fname[80];
FILE *fp_calib, *fp_out;
time_t t = time(NULL);


if(argc != 4) {
  printf("Syntax: merge_calib_resid calibrated_table residual_table out_table\n");
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
add_residuals_to_calib(fp_calib, resid_fname, fp_out); 

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
int add_residuals_to_calib(FILE *fp_calib, char *resid_fname, FILE *fp_out) 
{
char in_line[256], wds_name[40]; 
char object_name[40], ads_name[40], discov_name[40], comp_name[40]; 
char buffer[128], orbit_ref[60*50], quadrant_discrep[20]; 
char previous_object_name[40], previous_comp_name[40];
char sign_Drho[20], sign_Dtheta[20];
double epoch_o, rho_o, rho_o_c[50], theta_o_c[50];
int i, iline, norbits_found, status, orbit_grade, verbose_if_error = 0;
int ref_slength = 60, nmax_orbits = 50;

previous_object_name[0] = '\0';
previous_comp_name[0] = '\0';
iline = 0;
while(!feof(fp_calib)) {
  if(fgets(in_line, 256, fp_calib)) {
    iline++;
// Read header and add 3 columns for output header
    if(!strncmp(in_line,"& & & & & & & & & & & \\\\", 24)){
      fprintf(fp_out,   "& & & & & & & & & & & & & \\\\ \n");
    } 
// Good lines start with a digit (WDS names...) or with \idem
// Lines starting with % are ignored
    else if(in_line[0] != '%' && (isdigit(in_line[0]) 
        || !strncmp(in_line, "\\idem", 5))) {
// Remove the end of line '\n' from input line:
      jlp_cleanup_string(in_line, 256);

/* Old version : read orbit flag in column 11: 
      status = latex_get_column_item(in_line, buffer, 11, verbose_if_error);
      if(status || (sscanf(buffer,"%d", &orbit) != 1)) {
         fprintf(stderr,"Fatal error reading orbit flag: %s (line=%d)\n", 
                 in_line, iline);
         exit(-1);
         }
*/
// New version (after 2011): now search for orbit in all cases 
// since the orbit flag may be wrong

/* Removes the orbit column (column=11) */
      remove_orbit_column(in_line, 256);

/* Retrieve the residuals for this object and epoch in the residual table: */
/* Read object name and companion name from columns 2 and 3: */
      if(isdigit(in_line[0])) {
         read_object_name_from_CALIB_line(in_line, wds_name, discov_name, 
                                          comp_name, ads_name);
/* Generate the object name: */
        if(*ads_name) 
          strcpy(object_name, ads_name);
        else 
          strcpy(object_name, discov_name);
      } else {
         strcpy(object_name, previous_object_name);
         strcpy(comp_name, previous_comp_name);
      }

/* Read epoch in column 4: */
      status = latex_get_column_item(in_line, buffer, 4, verbose_if_error);
      if(status || (sscanf(buffer,"%lf", &epoch_o) != 1)) {
         fprintf(stderr,"Fatal error reading epoch: %s (line=%d)\n", 
                 in_line, iline);
         exit(-1);
         }

/* Read rho in column irho=5 if Gili's format: */
      rho_o = 0.;
      status = latex_get_column_item(in_line, buffer, 5, verbose_if_error);
      if(status || (sscanf(buffer,"%lf", &rho_o) != 1)) {
         fprintf(stderr,"Warning: error reading rho_o: %s (line=%d) Unres ?\n",
                 in_line, iline);
         }

/* Retrieve the residuals for this object and epoch in the residual table: */
      strcpy(quadrant_discrep,"");
      get_values_from_RESID_table(resid_fname, object_name, comp_name, epoch_o,
                                  rho_o, orbit_ref, &orbit_grade,
                                  ref_slength, rho_o_c, theta_o_c, 
                                  quadrant_discrep, &norbits_found, 
                                  nmax_orbits);
      if(!norbits_found) {
         fprintf(stderr," Residuals for %s %s not found in resid_table for epoch=%f in %s\n", 
                 object_name, comp_name, epoch_o, resid_fname);
/* Simply copy the input line to the output file: */
         fprintf(fp_out, "%s \\\\ \n", in_line);
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
           if(i == 0) {
           fprintf(fp_out, "%s & %s & %s%.2f & %s%.1f%s \\\\\n", 
                   in_line, &orbit_ref[i*ref_slength], 
                   sign_Drho, rho_o_c[i], sign_Dtheta,
                   theta_o_c[i], quadrant_discrep);
           } else { 
           fprintf(fp_out, "\\idem & \\idem & \\idem & \\idem & \
\\idem & \\idem & \\idem & \\idem & \\idem & \\idem &  & \
%s & %s%.2f & %s%.1f%s \\\\\n", 
                   &orbit_ref[i*ref_slength], 
                   sign_Drho, rho_o_c[i], sign_Dtheta,
                   theta_o_c[i], quadrant_discrep);
/*
           fprintf(fp_out, "\\idem & \\idem & \\idem & %.3f & %d \
% %.3f & %.3f & %.1f & %.1f & %s & 
%s & %s%.2f & %s%.1f%s \\\\\n", 
                   epoch_o, filter, eyepiece, rho_o, drho_o,
                   theta_o, dtheta_o, notes, &orbit_ref[i*ref_slength], 
                   sign_Drho, rho_o_c[i], sign_Dtheta,
                   theta_o_c[i], quadrant_discrep);
*/
           }
         } /* EOF for (i=0, norbits_found) */
         }
    } else if (in_line[0] == '%') {
/* Simply copy the input line to the output file if it is a comment: */
/* Remove the end of line '\n' from input line: */
      jlp_cleanup_string(in_line, 256);
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
jlp_compact_string(buffer, 80);

/* begin{tabular*} */
 if(!strncmp(buffer, "\\begin{tabular*}", 16)) {
/* JLP 2011: smaller to fit in page width */
    fprintf(fp_out, "\\small \n");
    fprintf(fp_out, "\\begin{tabular*}{\\textwidth}{clrcccccrcllrr} \n");
/* Extended header */
 } else if(!strncmp(buffer, "WDS", 3)) {
    fprintf(fp_out,"WDS & Name & ADS & Epoch & Fil. & Eyep. & $\\rho$ \
& $\\sigma_\\rho$ & \\multicolumn{1}{c}{$\\theta$} & $\\sigma_\\theta$ \
& Notes & Orbit & {\\scriptsize $\\Delta \\rho$(O-C)} \
& {\\scriptsize $\\Delta \\theta$(O-C)} \\\\ \n");
/* JLPPPP/DDEBUG: To be tested */
 } else if(!strncmp(buffer, "&&&&&(mm)", 9)){
    fprintf(fp_out,"& & & & & (mm) & (\\arcsec) & (\\arcsec) &  \\multicolumn{1}{c}{($^\\circ$)} \
& ($^\\circ$) & & & (\\arcsec) & ($^\\circ$) \\\\ \n");
 } else {
/* For the other lines, simply copy the input line to the output file: */
/* Remove the end of line '\n' from input line: */
  jlp_cleanup_string(in_line, 256);
  fprintf(fp_out, "%s\n", in_line);
 }
return(0);
}
/***********************************************************************
*
* Remove the 11th column corresponding to the orbit column
***********************************************************************/
static int remove_orbit_column(char *in_line, int ilen)
{
int verbose_if_error = 0;
char *pc, buffer[256], notes[80];

if(ilen > 256) {
  fprintf(stderr, "Fatal error: ilen=%d > 256! \n", ilen);
  exit(-1);
 }

/* Get notes from column 12: */
if(latex_get_column_item(in_line, notes, 12, verbose_if_error)){
  fprintf(stderr, "Fatal error reading the notes (in col. 12): bad syntax of line=%s \n", in_line);
  fprintf(stderr, "Fatal error reading the notes (in col. 12): notes=%s< \n", notes);
  exit(-1);
 }
/* Reduce length if full of ' ' ... */
jlp_trim_string(notes, 80);

strcpy(buffer, in_line);
/* Go to the end of Latex line (marked with "\\" or "\cr"): */
pc = buffer;
buffer[255] = '\0';
while(*pc && strncmp(pc, "\\cr", 3) && strncmp(pc, "\\\\", 2)) pc++;
if(!*pc) {
  fprintf(stderr, "Fatal error: bad syntax of line=%s \n", in_line);
  exit(-1);
 }
pc--;
while((pc != buffer) && (*pc != '&')) pc--;
pc--;
while((pc != buffer) && (*pc != '&')) pc--;
*pc = '\0';

sprintf(in_line, "%s & %s", buffer, notes);

return(0);
}

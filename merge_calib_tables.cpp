/************************************************************************
* merge_calib_tables.cpp
*
* To merge the two Latex calibrated tables
*
* JLP 
* Version 02/03/2023
*************************************************************************/
#include "jlp_catalog_utils.h"
#include "WDS_catalog_utils.h"
#include "jlp_calib_table.h"
#include "jlp_string.h"
#include "latex_utils.h"

/*
#define DEBUG
#define DEBUG_1
*/

static int merge_calib_tables(FILE *fp_calib1, char *calib2_fname, 
                              FILE *fp_out);

int main(int argc, char *argv[])
{
char calib1_fname[80], calib2_fname[80], out_fname[80];
FILE *fp_calib1, *fp_out;
time_t t = time(NULL);

if(argc != 4) {
  printf("Syntax: merge_calib_tables calib1_table calib2_table out_table \n");
  return(-1);
}
strcpy(calib1_fname, argv[1]);
strcpy(calib2_fname, argv[2]);
strcpy(out_fname, argv[3]);

printf("OK: calib=%s resid=%s output=%s \n", 
       calib1_fname, calib2_fname, out_fname); 

/* Open input calibrated table: */
if((fp_calib1 = fopen(calib1_fname, "r")) == NULL) {
   fprintf(stderr, "merge_calib_tables/Fatal error opening calib. table %s\n",
           calib1_fname);
   return(-1);
  }

/* Open output table: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "merge_calib_tables/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the Latex table: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "%% Merged table from: %s and %s \n%% Created on %s", 
        calib1_fname, calib2_fname, ctime(&t));
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

/* Scan the input calibrated table and add the residuals 
*/
merge_calib_tables(fp_calib1, calib2_fname, fp_out); 

/* Close opened files:
*/
fclose(fp_calib1);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the calibrated table and add the residuals 
*
* INPUT:
* fp_calib1: pointer to the input file containing the calibrated table 
* calib2_fname: name of the file containing the Latex table with the residuals
* fp_out: pointer to the output Latex file with the full table 
*
*************************************************************************/
static int merge_calib_tables(FILE *fp_calib1, char *calib2_fname, FILE *fp_out)
{
char in_line1[256], in_line2[256], buffer[256]; 
double WDS_alpha1, WDS_delta1, WDS_alpha2, WDS_delta2;
int iline1, iline2, status, verbose_if_error = 0, status2, nline2;
JLP_CalibTab *jlp_calib_tab2;

jlp_calib_tab2 = new JLP_CalibTab(calib2_fname);
status2 = jlp_calib_tab2->NextLineWithMeasurements(in_line2, 
                                                   &WDS_alpha2, &WDS_delta2); 
if(status2 != 0) {
  fprintf(stderr, "Fatal error reading %s\n", calib2_fname);
  exit(-1);
  }

// nline2 objects for one output page:
// 60: too many
nline2 = 50;
fprintf(fp_out, "\\jlpBeginTable0 \n");
iline1 = 0;
iline2 = 0;
while(!feof(fp_calib1)) {
  if(fgets(in_line1, 256, fp_calib1)) {
    iline1++;
// Remove the end of line '\n' from input line:
      jlp_cleanup_string(in_line1, 256);

/* Retrieve the residuals for this object and epoch in the residual table: */
/* Read object name and companion name from columns 2 and 3: */
      if(isdigit(in_line1[0])) {

/* Read WDS name in column 1: */
      status = latex_get_column_item(in_line1, buffer, 1, verbose_if_error);
//int decode_WDS_name(char *WDS_name, double *WDS_alpha, double *WDS_delta)
       decode_WDS_name(buffer, &WDS_alpha1, &WDS_delta1);
// Check if line2 should be placed before line1:
       while((WDS_alpha1 > WDS_alpha2) && (status2 == 0)) {
          fprintf(fp_out, "%s \n", in_line2);
          iline2++;
          if(iline2 % nline2 == 0) fprintf(fp_out, "\\jlpEndTable \n\\jlpBeginTable \n");
          status2 = jlp_calib_tab2->NextLineWithMeasurements(in_line2, 
                                                   &WDS_alpha2, &WDS_delta2); 
          }

/* Copy the input line to the output file: */
      fprintf(fp_out, "%s \n", in_line1);
      if(iline2 % nline2 == 0) fprintf(fp_out, "\\jlpEndTable \n \\jlpBeginTable \n");
      iline2++;
      }/* EOF if !isdigit ... */
  } /* EOF if fgets */ 
} /* EOF while ... */
printf("merge_calib_tables: %d lines sucessfully read and processed (outfile: %d lines)\n", 
        iline1, iline2);
fprintf(fp_out, "\\jlpEndTableE \n");
return(0);
}

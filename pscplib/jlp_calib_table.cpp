/********************************************************************
* jlp_calib_table.ccp
*
* JLP
* Version 03/03/2023
*********************************************************************/
#include  <stdio.h>
#include  <ctype.h>   // isdigit()
#include  <string.h>  // strcpy()
#include "jlp_calib_table.h"
#include "jlp_string.h"
#include "latex_utils.h"
#include "WDS_catalog_utils.h"

/***********************************************************************
* Constructor 
************************************************************************/
JLP_CalibTab::JLP_CalibTab(char *calib_table_fname0) {

 strcpy(calib1_fname, calib_table_fname0);
 iline1 = 0;

 /* Open input calibrated table: */
if((fp_calib1 = fopen(calib1_fname, "r")) == NULL) {
   fprintf(stderr, "JLP_CalibTab/Fatal error opening calib. table %s\n",
           calib1_fname);
   return;
  }

}
/***********************************************************************
* Look for next line with double star measurement:
************************************************************************/
int JLP_CalibTab::NextLineWithMeasurements(char *full_line, 
                             double *WDS_alpha, double *WDS_delta) {
int status = -1, verbose_if_error = 0;;
char in_line[256], buffer[64];

full_line[0] = '\0';
*WDS_alpha = 0.;
*WDS_delta = 0.;
if(fp_calib1 == NULL) return(-1);

while(!feof(fp_calib1)) {
  if(fgets(in_line, 256, fp_calib1)) {
  iline1++;

// Remove the end of line '\n' from input line:
      jlp_cleanup_string(in_line, 256);

/* Retrieve the residuals for this object and epoch in the residual table: */
/* Read object name and companion name from columns 2 and 3: */
      if(isdigit(in_line[0])) {

/* Read WDS name in column 1: */
      status = latex_get_column_item(in_line, buffer, 1, verbose_if_error);
//int decode_WDS_name(char *WDS_name, double *WDS_alpha, double *WDS_delta)
      if(status == 0) {
         decode_WDS_name(buffer, WDS_alpha, WDS_delta);
         strcpy(full_line, in_line);
         return(0);
         }
       } // if isdigit()
    } // ifgets
}

return(status);
}

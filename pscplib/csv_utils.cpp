/*************************************************************************
* csv_utils.c
* JLP
* Version 02/03/2020
*************************************************************************/
#include <ctype.h>          /* isdigit() */
#include <string.h>         /* strstr() */
#include "astrom_utils1.h"
#include "astrom_utils2.h"
#include "PISCO_catalog_utils.h" /* get_data_from_PISCO_catalog */
#include "WDS_catalog_utils.h" /* get_data_from_WDS_catalog */
#include "jlp_fitsio.h"  // get_bessel_epoch_from_fits_file()
#include "jlp_numeric.h" // JLP_QSORT, MINI, MAXI, PI, etc
#include "jlp_string.h"    // jlp_trim_string, jlp_compact_string 
#include "latex_utils.h"

#include "csv_utils.h" // prototypes defined here
/*
#define DEBUG1
*/
#define DEBUG
/*****************************************************************************
* Read the measurements from the input line from Gili's csv file 
* Example:
,,,,,,,,,,,,,,
 a830\_a , 04/01/2013 ,  ,1,"4,01","0,1","78,41","1,8"," EP=2013,0129 Q=2 dm=0,01 ",,,"0,295938","168,41","169,4",
 a920\_a , 04/01/2013 ,  ,1,"26,11","0,1","-43,07","0,4"," EP=2013,0125 Q=3? ",,1,"1,926918","226,93","227,92",

*
* INPUT:
* nobj: number of objects already entered into *obj 
*
* OUTPUT:
* obj: OBJECT structure 
* nobj: total number of objects entered into *obj
*****************************************************************************/
int csv_read_gili_measures_from_line(char *b_data, OBJECT *obj1, int *nobj1,
                                     double scale_mini)
{
int i_WDSName, i_ObjectName, i_filename, i_eyepiece; 
int i_rho, i_rho_pixels, i_theta, i_drho, i_dtheta, i_notes;
int status = -1, iverbose = 0;
double epoch1, rho1, theta1, drho1, dtheta1, dmag1, ddmag1, ww, rho1_pixels;
double scale_arcsec_pixel;
int eyepiece1, quadrant1, dquadrant1;
char WDSName1[64], ObjectName1[64], filter1[32], notes1[64];
char *pc, buffer[256], out_string[64];
MEASURE *me;


if(b_data[0] == ',') {
  return(-1);
  }

eyepiece1 = 0;
quadrant1 = 0;
dquadrant1 = 0;
epoch1 = NO_DATA;
rho1 = NO_DATA;
drho1 = NO_DATA;
theta1 = NO_DATA;
dtheta1 = NO_DATA;
dmag1 = NO_DATA;
ddmag1 = NO_DATA;

/*
* i_filename: column nber of the filename used for this measurement
* i_rho: column nber with rho values
* i_theta: column nber with theta values
* i_notes: column nber with the notes
*/

i_WDSName = 1;

// Gili scale: should multiply rho_pixels in col 5 
// by 0.0738 to obtain rho_arcsec in col 12 for autocorrelations (*_a, *_aw,..)
// by 0.455 to obtain rho_arcsec in col 12 for long exposures (*_l)
/*
 a830\_a , 04/01/2013 ,  ,1,"4,01","0,1","78,41","1,8"," EP=2013,0129 Q=2 dm=0,01 ",,,"0,295938","168,41","169,4",
 lds1074\_l , 03/12/2013 ,  ,1,"10,52","0,1","146,46","0,3"," EP=2013,9244  Dm=3,53+/-0,02",,,"4,7866","236,46","237,45",
*/

(obj1[*nobj1]).nmeas = 1;
me = &(obj1[*nobj1]).meas[0];

#ifdef DEBUG1
printf("DEBUG1/csv_read_gili_measures_from_line \n");
printf("in: b_data=%s\n", b_data);
#endif
// Decode the object name:
i_ObjectName = 1;
status = csv_read_object_from_filename(b_data, i_ObjectName, ObjectName1);
if(status == 0) {
  strcpy((obj1[*nobj1]).discov_name, ObjectName1);
  }

// Decode the rho measurement:
i_rho = 12; 
status = csv_read_dvalue(b_data, i_rho, &rho1);
if(status == 0) 
  me->rho = rho1;
else
  return(-1);

// Decode rho in pixels:
i_rho_pixels = 5; 
scale_arcsec_pixel = 0.;
status = csv_read_dvalue(b_data, i_rho_pixels, &rho1_pixels);
if(status == 0) {
  if(rho1_pixels != 0) scale_arcsec_pixel = rho1 / rho1_pixels;
  }

if(scale_arcsec_pixel > 0.4) 
  me->eyepiece = 5;
else if(scale_arcsec_pixel > 0.08) 
  me->eyepiece = 2;
else 
  me->eyepiece = 1;
#ifdef DEBUG1
printf("out: rho1=%f eyepiece=%d \n", rho1, me->eyepiece);
#endif

// Decode the drho measurement:
i_drho = 6; 
status = csv_read_dvalue(b_data, i_drho, &drho1);
// Gili scale: should multiply rho_pixels in col 5 
// by 0.0738 to obtain rho_arcsec in col 12 for autocorrelations (*_a)
// by 0.455 to obtain rho_arcsec in col 12 for long exposures (*_l)
drho1 *= scale_arcsec_pixel;
if(status == 0) 
  me->drho = drho1;

// Decode the theta measurement:
i_theta = 14; 
status = csv_read_dvalue(b_data, i_theta, &theta1);
if(status == 0) me->theta = theta1;
#ifdef DEBUG1
printf("out: theta1=%f status=%d \n", theta1, status);
#endif

// Decode the dtheta measurement:
i_dtheta = 8; 
status = csv_read_dvalue(b_data, i_dtheta, &dtheta1);
if(status == 0) me->dtheta = dtheta1;

// Decode the notes:
i_notes = 9; 
status = csv_read_string(b_data, i_notes, notes1);
if(status == 0) strcpy((obj1[*nobj1]).notes, notes1);

// Decode the epoch from the notes::
status = csv_read_epoch_from_notes(b_data, i_notes, &epoch1);
if(status == 0) me->epoch = epoch1;
#ifdef DEBUG1
printf("out: epoch1=%f status=%d \n", epoch1, status);
#endif

// Decode Dmag from the notes::
status = csv_read_dmag_from_notes(b_data, i_notes, &dmag1);
if(status == 0) me->dmag = dmag1;
#ifdef DEBUG1
printf("out: dmag1=%f status=%d\n", dmag1, status);
#endif

#ifdef DEBUG
   printf("bdata: %s", b_data);
   printf("object: >%s< , notes: >%s<\n", ObjectName1, notes1);
   printf("rho=%.3f+/-%.3f (rho=%.2f px) theta=%.1f+/-%.1f (scale=%.4f, eyep=%d) epoch=%.4f Dm=%.2f\n", 
          rho1, drho1, rho1_pixels, theta1, dtheta1, scale_arcsec_pixel, 
          me->eyepiece, epoch1, dmag1);
#endif

// EXIT for DEBUG
// if(*pc != 1234) exit(0);

/***
status = latex_read_svalue(b_data, WDSName1, i_WDSName);
if(status == 0) strcpy((obj1[*nobj1]).wds, WDSName1);

***/
// Gili scale: should multiply rho_pixels in col 5 
// by 0.0738 to obtain rho_arcsec in col 12 for autocorrelations (*_a, *_aw,..)
// by 0.455 to obtain rho_arcsec in col 12 for long exposures (*_l)
// Load only the measures with a problem:
// Typically: scale_mini=0.1
if((rho1 != NO_DATA) && (scale_arcsec_pixel > scale_mini)) (*nobj1)++;

return(0);
}
/*************************************************************
* Decode the epoch from the notes:
*************************************************************/
int csv_read_epoch_from_notes(char *b_data, int i_notes, double *epoch1)
{
char buffer[64], *pc, epoch_title[64], *epoch_str;
int status = -1, nval;

*epoch1 = 0.;

// Read the notes:
status = csv_read_string(b_data, i_notes, buffer);

// First look for "EP=":
strcpy(epoch_title, "EP=");
epoch_str = strstr(buffer, epoch_title);

status = -1;
if(epoch_str) {
// Begins at EP=
  epoch_str += 3;
  strcpy(buffer, epoch_str);
  pc = buffer;
// From ',' to '.':
  while (*pc) {if(*pc == ',') *pc = '.'; pc++;}

  while (*pc) {if(*pc == ',') *pc = '.'; pc++;}

  nval = sscanf(buffer, "%lf", epoch1);
  if(nval == 1) status = 0;
}

return(status);
}
/*************************************************************
* Decode dmag from the notes:
*************************************************************/
int csv_read_dmag_from_notes(char *b_data, int i_notes, double *dmag1)
{
char buffer[64], *pc, dmag_title[64], *dmag_str;
int status = -1, nval;

*dmag1 = 0.;

// Read the notes:
status = csv_read_string(b_data, i_notes, buffer);
// printf(" dm_buffer=%s\n", buffer);

// From 'd' to 'D':
pc = buffer;
while(*pc) {if(*pc == 'd') *pc = 'D'; pc++;}

// printf(" Dm_buffer=%s\n", buffer);
// First look for "Dm=":
strcpy(dmag_title, "Dm=");
dmag_str = strstr(buffer, dmag_title);

status = -1;
if(dmag_str) {

// Begins at Dm=
  dmag_str += 3;
  strcpy(buffer, dmag_str);
  pc = buffer;
//From ',' to '.':
  while (*pc) {if(*pc == ',') *pc = '.'; pc++;}

  while (*pc) {if(*pc == ',') *pc = '.'; pc++;}

  nval = sscanf(buffer, "%lf", dmag1);
  if(nval == 1) status = 0;
}

return(status);
}
/*************************************************************
* Decode the object name:
*************************************************************/
int csv_read_object_from_filename(char *b_data, int i_column, char *ObjectName1)
{
char buffer[64], *pc, out_string[64];
int icol, status = -1;

status = csv_read_string(b_data, i_column, out_string);

// Example: csv_read_string/out_string: > a120\_a <

 strcpy(buffer, out_string);

// Decode string: stop at "\_"
 pc = buffer;
 while(*pc && (*pc != '\\')) pc++;
 *pc = '\0';

 if(buffer[0] == ' ') strncpy(out_string, &buffer[1], 64);
   else strncpy(out_string, buffer, 64);

// Conversion to upper case:
 pc = out_string;
 while(*pc) { *pc = toupper(*pc); pc++; }
  if(out_string[0] == ' ') strcpy(ObjectName1, &out_string[1]);
     else strcpy(ObjectName1, out_string);

if(out_string[0] != '\0') status = 0;

return(status);
}
/*************************************************************
* Read a string from the input line 
*************************************************************/
int csv_read_string(char *b_data, int i_column, char *out_string)
{
char buffer[256], buf1[256], *pc;
int icol, iquote, status = -1;

/* DEBUG
printf("i_column=%d\n", i_column);
printf("b_data: >%s<\n", b_data);
*/

strcpy(buf1, b_data);
pc = buf1;
// Go to the wright column:
icol = 1;
iquote = 0;
 while(*pc) {
   if(*pc == '"') if(iquote == 0) iquote = 1; else iquote = 0; 
   if((iquote == 0) && (*pc == ',')) icol++; 
   pc++;
   if(i_column == icol) break;
   }
strcpy(buffer, pc);

// Decode string: stop at ','
pc = buffer;
iquote = 0;
 while(*pc) {
   if(*pc == '"') if(iquote == 0) iquote = 1; else iquote = 0; 
   if((iquote == 0) && (*pc == ',')) {
      *pc = '\0'; 
      break;
      }
   pc++;
  }
 if(*pc == ',') *pc = '\0';

strncpy(out_string, buffer, 64);
if(out_string[0] != '\0') status = 0;

#ifdef DEBUG1
printf("csv_read_string/out_string: >%s<\n", out_string);
#endif

return(status);
}
/*************************************************************
* Read a double value from the input line 
*************************************************************/
int csv_read_dvalue(char *b_data, int i_column, double *dvalue)
{
char buffer[64], *pc, out_string[64];
int nval, status = -1;

status = csv_read_string(b_data, i_column, out_string);

// Removes the quotes at the beginning if any: 
 if(out_string[0] == '"') strcpy(buffer, &out_string[1]);
   else strcpy(buffer, out_string);

 pc = buffer;
 while(*pc && (*pc != '"') ) 
    {
    if(*pc == ',') *pc = '.';
    pc++;
    }
 *pc = '\0';

nval = sscanf(buffer, "%lf", dvalue);
if(nval == 1) status = 0;
 else status = -1;

#ifdef DEBUG1
printf("dvalue=%f status=%d\n", *dvalue, status);
#endif

return(status);
}
/*****************************************************************************
* Read the measurements and object parameters from the input tex calib file 
*
* INPUT:
* filein1: input LaTeX file with calibrated data
* nobj1: number of objects already entered into *obj
* with_wds_data: flag set to one if wds data are in the input file
*
* OUTPUT:
* obj1: OBJECT structure
* nobj1: total number of objects entered into *obj1
*****************************************************************************/
int csv_read_gili_measures(char *filein1, OBJECT *obj1, int *nobj1, 
                           int nobj_maxi, double scale_mini)
{
char b_in[NMAX], wds_name[40], discov_name[40]; 
int status, iline;
char *pc, *pc1;
FILE *fp_in;

if((fp_in = fopen(filein1,"r")) == NULL) {
  fprintf(stderr, " Fatal error opening input file1 %s \n", filein1);
  return(-1);
  }

#ifdef DEBUG1
printf("csv_read_gili_measures: input/nobj=%d\n", *nobj1);
#endif

wds_name[0] = '\0';
discov_name[0] = '\0';
iline = 0;
*nobj1 = 0;
while(!feof(fp_in))
{
/* Maximum length for a line will be 170 characters: */
  if(fgets(b_in,170,fp_in))
  {
  b_in[169] = '\0';
/* NEW/2009: I remove ^M (Carriage Return) if present: */
  pc = b_in;
  while(*pc) {
  if(*pc == '\r') {*pc = '\0'; break;}
   pc++;
   }
  *pc = '\0';
  iline++;

/*
int tex_calib_read_measures(FILE *fp_in, int comments_wanted, OBJECT *obj,
                            int *nobj);
*/
 status = csv_read_gili_measures_from_line(b_in, obj1, nobj1, scale_mini);
if(*nobj1 > nobj_maxi-2) {
  fprintf(stderr, "Fatal exit: nobj > maxi=%d - 2\n", nobj_maxi);
  }

  } /* EOF if fgets() */
} /* EOF while loop */
printf("csv_read_gili_measures/Number of lines: nlines=%d nobj=%d\n", iline, *nobj1);

fclose(fp_in);
return(0);
}

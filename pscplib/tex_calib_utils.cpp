/*************************************************************************
* tex_calib_utils.c
* JLP
* Version 02/03/2020
*************************************************************************/
#include "ctype.h"            /* isdigit() */
#include "astrom_utils1.h"
#include "astrom_utils2.h"
#include "PISCO_catalog_utils.h" /* get_data_from_PISCO_catalog */
#include "WDS_catalog_utils.h" /* get_data_from_WDS_catalog */
#include "jlp_fitsio.h"  // get_bessel_epoch_from_fits_file()
#include "jlp_numeric.h" // JLP_QSORT, MINI, MAXI, PI, etc
#include "jlp_string.h"    // jlp_trim_string, jlp_compact_string 
#include "latex_utils.h"

#include "tex_calib_utils.h" // prototypes defined here
/*
#define DEBUG
#define DEBUG1
*/
/************************************************************************
*
************************************************************************/
int extract_companion_from_name(char *input_object_name, char *object_name0, 
                                char *comp_name0, int *object_len0)
{
int len;
char *pc;

strcpy(object_name0, input_object_name);
// Assume length of object_name0 is 64:
jlp_compact_string(object_name0, 64);

/* Extract companion name if present: */
  pc = object_name0;
  len = 0;
  while(*pc && !isdigit(*pc)) { pc++; len++; }
  while(*pc && isdigit(*pc)) { pc++; len++; }
  strcpy(comp_name0, pc);
  *pc = '\0';

#ifdef DEBUG1
printf("DEBUG1/input>%s< object >%s< companion >%s< len=%d\n", 
       input_object_name, object_name0, comp_name0, len);
#endif
*object_len0 = len;

return (0);
}

/*****************************************************************************
* Read the measurements from the input tex calib line from file 
* Example:
00014+3937  &  HLD60  & 2014.828 & 1 & 1.322 & 0.007 & 167.5\rlap{$^*$} & 0.3 &  & & Izm2019 & 0.01 & 0.3 \\
*
* INPUT:
* nobj: number of objects already entered into *obj 
* in_calib_fmt : 1 = (gili format) if no filter column
*                2 = (Calern format) if filter column
*
* OUTPUT:
* obj: OBJECT structure 
* nobj: total number of objects entered into *obj
*****************************************************************************/
int tex_calib_read_measures_from_line(char *b_data, OBJECT *obj1,
                                      int *nobj1, int *nunres1, 
                                      int in_calib_fmt)
{
int i_WDSName, i_ObjectName, i_eyepiece, i_filter; 
int i_epoch, i_rho, i_drho, i_theta, i_dtheta, i_dmag;
int status = -1, iverbose = 0;
double epoch1, rho1, drho1, theta1, dtheta1, dmag1, ddmag1, ww;
int eyepiece1, quadrant1, dquadrant1;
char WDSName1[64], ObjectName1[64], filter1[32], notes1[64];
MEASURE *me;

if(isdigit(b_data[0]) == 0) {
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
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
* i_notes: column nber with the notes
* i_dmag: column nber with dmag values
*/

i_WDSName = 1;
i_ObjectName = 2;
i_epoch = 3;
// Gili's format: no filter 
//printf("in_calib_fmt=%d (=1 if Gili, ie, without filter)\n", in_calib_fmt);
if(in_calib_fmt == 1) {
  i_filter = -1;
  i_eyepiece = 4;
// Calern format: filter 
 } else {
  i_filter = 4;
  i_eyepiece = 5;
 }
i_rho = i_eyepiece + 1; 
i_drho = i_rho + 1; 
i_theta = i_drho + 1; 
i_dtheta = i_theta + 1; 
i_dmag = i_dtheta + 1;

(obj1[*nobj1]).nmeas = 1;
me = &(obj1[*nobj1]).meas[0];
strcpy((obj1[*nobj1]).wds, "");
strcpy((obj1[*nobj1]).discov_name, "");
me->rho = NO_DATA;
me->drho = NO_DATA;
me->theta = NO_DATA;
me->dtheta = NO_DATA;
me->epoch = NO_DATA;
me->dmag = NO_DATA;
status = latex_read_svalue(b_data, WDSName1, i_WDSName);
if(status == 0) strcpy((obj1[*nobj1]).wds, WDSName1);
status = latex_read_svalue(b_data, ObjectName1, i_ObjectName);
if(status == 0) strcpy((obj1[*nobj1]).discov_name, ObjectName1);
if(i_filter > 0) {
  status = latex_read_svalue(b_data, filter1, i_filter);
  if(status == 0) strcpy(me->filter, filter1);
  }

status = latex_read_dvalue(b_data, &epoch1, i_epoch, iverbose);
if(status == 0) me->epoch = epoch1;
status = latex_read_dvalue(b_data, &rho1, i_rho, iverbose);
if((status != 0) || (rho1 == NO_DATA)) (*nunres1)++;
if(status == 0) { 
   me->rho = rho1;
#ifdef DEBUG1
   printf("bdata: %s", b_data);
   printf("rho=%.2f (NO_DATA=%.2f) i_rho=%d, status=%d\n", 
          rho1, NO_DATA, i_rho, status);
#endif
  } 
status = latex_read_dvalue(b_data, &drho1, i_drho, iverbose);
if(status == 0) me->drho = drho1;
status = latex_read_dvalue(b_data, &theta1, i_theta, iverbose);
if(status == 0) me->theta = theta1;
status = latex_read_dvalue(b_data, &dtheta1, i_dtheta, iverbose);
if(status == 0) me->dtheta = dtheta1;
status = latex_read_dvalue(b_data, &dmag1, i_dmag, iverbose);
if(status == 0) me->dmag = dmag1;
status = latex_read_dvalue(b_data, &ww, i_eyepiece, iverbose);
if(status == 0) {
   eyepiece1 = (int)(ww+0.5);
   me->eyepiece = eyepiece1;
   }
(*nobj1)++;

return(0);
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
int tex_calib_read_measures(char *filein1, OBJECT *obj1, int *nobj1, 
                            int *nunres1, int in_calib_fmt)
{
char b_in[NMAX], wds_name[40], discov_name[40]; 
int status, iline;
char *pc, *pc1;
FILE *fp_in;

if((fp_in = fopen(filein1,"r")) == NULL) {
  fprintf(stderr, " Fatal error opening input file1 %s \n", filein1);
  return(-1);
  }

#ifdef DEBUG
printf("tex_calib_read_measures: input/nobj=%d\n", *nobj1);
#endif

wds_name[0] = '\0';
discov_name[0] = '\0';
iline = 0;
*nobj1 = 0;
*nunres1 = 0;
while(!feof(fp_in))
{
/* Maximum length for a line will be 170 characters: */
  if(fgets(b_in,170,fp_in))
  {
  b_in[169] = '\0';
/* NEW/2009: I remove ^M (Carriage Return) if present: */
  pc = b_in;
  while(*pc) {
  if(*pc == '\r') *pc = ' ';
   pc++;
   }
  iline++;

/*
int tex_calib_read_measures(FILE *fp_in, int comments_wanted, OBJECT *obj,
                            int *nobj, int in_calib_fmt);
*/
 status = tex_calib_read_measures_from_line(b_in, obj1, nobj1, nunres1,
                                            in_calib_fmt);

  } /* EOF if fgets() */
} /* EOF while loop */
printf("Number of lines: nlines=%d nobj=%d\n", iline, *nobj1);

fclose(fp_in);
return(0);
}
/***************************************************************************
* tex_calib_compare
*
* INPUT:
*   rho3, err_rho3....
*
* OUTPUT:
*   compatible_meas31: 1 if compatible, 0 if not compatible
***************************************************************************/
int tex_calib_compare(double rho3, double err_rho3, double theta3,
                      double err_theta3, double rho1, double err_rho1,
                      double theta1, double err_theta1, 
                      double max_drho, double max_dtheta, 
                      int *compatible_meas31)
{
int status = -1;
double D_rho, D_theta;
 *compatible_meas31 = 0;
 D_rho = ABS(rho3 - rho1);
 if(theta3 > 180.) theta3 -= 180.;
 if(theta1 > 180.) theta1 -= 180.;
 D_theta = ABS(theta3 - theta1);
 if((D_rho < max_drho)
    && (D_theta < max_dtheta) ) {
    *compatible_meas31 = 1;
    status = 0;
    }
 if((rho1 == NO_DATA) || (theta1 == NO_DATA)
    || (rho3 == NO_DATA) || (theta3 == NO_DATA)) status = -1;

return(status);
}
/***************************************************************************
* Get measures from OBJECT
*
* INPUT:
*  WDSName0, ObjectName0
*
* OUTPUT:
*  WDSName1, ObjectName1
*  epoch1, rho1, drho1, theta1, dtheta1
***************************************************************************/
int tex_calib_get_meas_from_object(OBJECT *obj1, int nobj1, 
                                   char *WDSName0, char *ObjectName0,
                                   double epoch_obs0,
                                   char *WDSName1, char *ObjectName1,
                                   double *epoch1, char *filt1, int *eyep1, 
                                   double *rho1, double *err_rho1, 
                                   double *theta1, double *err_theta1)
{
char wds_name[64], discov_name[64], object_name0[64], comp_name0[64], *pc;
char discov_name1[64], discov_comp1[64];
double epoch_obs;
MEASURE *me;
int i, status = -1, ifound = -1, len0, len1, maxlen0;

strcpy(object_name0, ObjectName0);

/* Extract companion name if present: */
extract_companion_from_name(ObjectName0, object_name0, comp_name0, &len0);

// First go looking for WDSName:
for(i = 0; i < nobj1; i++) {
  strcpy(wds_name, obj1[i].wds);
  if(strcmp(WDSName0, wds_name) == 0){
    ifound = i;
    break;
    }
  }
// Second go looking for WDSName && ObjectName (only the first maxlen0 characters):
for(i = 0; i < nobj1; i++) {
  strcpy(wds_name, obj1[i].wds);
  strcpy(discov_name, obj1[i].discov_name);
  extract_companion_from_name(discov_name, discov_name1, discov_comp1, &len1);
  maxlen0 = MAXI(len0, len1);
  if((strcmp(WDSName0, wds_name) == 0)
     && (strncmp(object_name0, discov_name, maxlen0) == 0)){
    ifound = i;
    break;
    }
  }
// Third go looking for WDSName && ObjectName && Epoch:
for(i = 0; i < nobj1; i++) {
  strcpy(wds_name, obj1[i].wds);
  strcpy(discov_name, obj1[i].discov_name);
  extract_companion_from_name(discov_name, discov_name1, discov_comp1, &len1);
  maxlen0 = MAXI(len0, len1);
  epoch_obs = (obj1[i].meas[0]).epoch;
  if((strcmp(WDSName0, wds_name) == 0)
     && (strncmp(object_name0, discov_name, maxlen0) == 0)
     && (ABS(epoch_obs - epoch_obs0) < 0.01) ){
    ifound = i;
    break;
    }
  }
if(ifound >= 0) {
  strcpy(WDSName1, obj1[ifound].wds);
  strcpy(ObjectName1, obj1[ifound].discov_name);
  me = &(obj1[ifound]).meas[0];
  *epoch1 = me->epoch;
  strcpy(filt1, me->filter);
  *eyep1 = me->eyepiece;
  *rho1 = me->rho;
  *err_rho1 = me->drho;
  *theta1 = me->theta;
  *err_theta1 = me->dtheta;
  status = 0;
  }
return(status);
}
/***************************************************************************
* Get measures from OBJECT
* Same as tex_calib_from_object, but with fewer parameters when object was created with Gili's csv file 
*
* INPUT:
*  WDSName0, ObjectName0
*
* OUTPUT:
*  WDSName1, ObjectName1
*  epoch1, rho1, drho1, theta1, dtheta1
***************************************************************************/
int tex_calib_get_meas_from_csv_object(OBJECT *obj1, int nobj1, 
                                   char *ObjectName0, double epoch_obs0,
                                   char *ObjectName1, double *epoch1,
                                   int *eyep1, double *rho1, double *err_rho1, 
                                   double *theta1, double *err_theta1,
                                   double *dmag1)
{
char discov_name[64], object_name0[64], comp_name0[64];
char discov_name1[64], discov_comp1[64];
double epoch_obs;
MEASURE *me;
int i, status = -1, ifound = -1, len0, len1, maxlen0;

strcpy(object_name0, ObjectName0);
jlp_compact_string(object_name0, 64);

/* Extract companion name if present: */
extract_companion_from_name(ObjectName0, object_name0, comp_name0, &len0);

// First go looking for ObjectName (first maxlen0 characters):
for(i = 0; i < nobj1; i++) {
  strcpy(discov_name, obj1[i].discov_name);
  jlp_compact_string(discov_name, 64);
  extract_companion_from_name(discov_name, discov_name1, discov_comp1, &len1);
  maxlen0 = MAXI(len0, len1);
  if(strncmp(object_name0, discov_name, maxlen0) == 0){
    ifound = i;
    break;
    }
  }

// Second go looking for ObjectName && Epoch:
for(i = 0; i < nobj1; i++) {
  strcpy(discov_name, obj1[i].discov_name);
  jlp_compact_string(discov_name, 64);
  extract_companion_from_name(discov_name, discov_name1, discov_comp1, &len1);
  maxlen0 = MAXI(len0, len1);
  epoch_obs = (obj1[i].meas[0]).epoch;
  if((strncmp(object_name0, discov_name, maxlen0) == 0)
     && (ABS(epoch_obs - epoch_obs0) < 0.01) ){
    ifound = i;
    break;
    }
  }

if(ifound >= 0) {
  strcpy(ObjectName1, obj1[ifound].discov_name);
  printf("UUUU/ObjectName1=>%s< object_name=>%s< \n", ObjectName1, object_name0);
  me = &(obj1[ifound]).meas[0];
  *epoch1 = me->epoch;
  *eyep1 = me->eyepiece;
  *rho1 = me->rho;
  *err_rho1 = me->drho;
  *theta1 = me->theta;
  *err_theta1 = me->dtheta;
  *dmag1 = me->dmag;
  status = 0;
  }

return(status);
}

/*****************************************************************************
* Write the LateX array in a new format 
*
*****************************************************************************/
int astrom_write_publi_table(FILE *fp_out, int comments_wanted, OBJECT *obj, 
                             int *index_obj, int nobj, int tabular_only)
{
MEASURE *me, *me_next;
int good_q, nm, io, nlines, jj, jnext, nl_max, status;
char qflag[20], asterisk[20], exclam[20], q_error[1], wds_name[40], q_notes[20];
register int i, j;

/* Portrait: nl_max=55 (MNRAS for one page) */
/* Landscape: nl_max= 30 */
 if(comments_wanted)
   nl_max = 30;
 else
   nl_max = 53;

strcpy(asterisk,"\\rlap{$^*$}");
strcpy(exclam,"\\rlap{!}");

fprintf(fp_out,"\\def\\idem{''} \n");

nlines = -1;

for(i = 0; i < nobj; i++) {
  io = index_obj[i];
  nm = (obj[io]).nmeas;
/* New table header in publi_mode */
  if(nlines > nl_max) {
    fprintf(fp_out,"& & & & & & & & & \\\\ \n");
    fprintf(fp_out,"\\hline \n");
    fprintf(fp_out,"\\end{tabular*} \n");
/* Problem for too big tables, so I only use tabular for big tables: */
    if(!tabular_only) fprintf(fp_out,"\\end{table*} \n");
  }
  if((nlines == -1) || (nlines > nl_max)) {
/* JLP to obtain the same table number:
*/
/* Problem for too big tables, so I only use tabular for big tables: */
    if(!tabular_only) {
      if(nlines != -1) fprintf(fp_out,"\\addtocounter{table}{-1}\n");
      fprintf(fp_out,"\\begin{table*} \n");
      }
    fprintf(fp_out,"\\tabcolsep=1mm \n");
//    fprintf(fp_out,"\\begin{tabular*}{\\textwidth}{clrcccccrccl} \n");
    fprintf(fp_out,"\\begin{tabular*}{\\textwidth}{clcccccrcl} \n");
    fprintf(fp_out,"\\hline \n");
    fprintf(fp_out,"& & & & & & & & & \\\\ \n");
//    fprintf(fp_out,"WDS & Name & ADS & Epoch & Fil. & Eyep. & $\\rho$ & $\\sigma_\\rho$ & \\multicolumn{1}{c}{$\\theta$} & $\\sigma_\\theta$ & Orb. & Notes \\\\ \n");
    fprintf(fp_out,"WDS & Name & Epoch & Fil. & Eyep. & $\\rho$ & $\\sigma_\\rho$ & \\multicolumn{1}{c}{$\\theta$} & $\\sigma_\\theta$ & Notes \\\\ \n");
//    fprintf(fp_out,"& &       &       &        & (mm)& (\") & (\") & ($^\\circ$) & ($^\\circ$) & & \\\\ \n");
    fprintf(fp_out,"& &       &        & (mm)& (\") & (\") & ($^\\circ$) & ($^\\circ$) & \\\\ \n");
    fprintf(fp_out,"& & & & & & & \\\\ \n");
    fprintf(fp_out,"\\hline \n");
    fprintf(fp_out,"& & & & & & & \\\\ \n");
    nlines = 0;
  }

  jj = 0;
/************* Loop on measurements: **********************/
  for(j = 0; j < nm; j++) {
    me = &(obj[io]).meas[j];
/* BOF case not_flagged */
    if(!me->flagged_out) {
/* For DEBUG:
if(i < 4) {
  printf("i= %d io=%d nm=%d j=%d flag=%d\n", i, io, nm, j, me->flagged_out);
  printf(" WDS=%s rho=%.2f drho=%.2f theta=%.2f dtheta=%.2f eyepiece=%d\n", 
           obj[io].wds, me->rho, me->drho, me->theta, me->dtheta, me->eyepiece);
  }
*/

/* JLP2007: I check if among the next measurements there exists
* another one with the same epoch and filter
* than the current measurement: if it is the case, compute
* the mean of the two: */
    for(jnext = j+1; jnext < nm; jnext++) {
       me_next = &(obj[io]).meas[jnext];
       if(!me_next->flagged_out && (ABS(me->epoch - me_next->epoch) < 0.001) 
         && !strncmp(me->filter,me_next->filter,3)
         && (me->eyepiece == me_next->eyepiece)
         && (me_next->rho != NO_DATA) && (me_next->theta != NO_DATA) ) {
/* Too many messages in case of error:
         printf("astrom_write_publi/Warning same filter, epoch and eyepiece for object %s: I compute the mean (j=%d, jnext=%d)\n",
          obj[io].wds, j, jnext);
         printf("(filter=%s, epoch=%.4f (me_epoch=%.4f), eyepiece=%d)\n", me->filter, me_next->epoch,
                 me->epoch, me->eyepiece);
         printf("ABS=%f\n", ABS(me->epoch - me_next->epoch)); 
*/
         astrom_compute_mean_of_two_measures(obj, io, j, jnext);
         }
      }

  good_q = astrom_quadrant_is_consistent(me);
/* Case when quadrant was found and theta is OK: */
  if(good_q == 1) 
    strcpy(qflag,asterisk);
/* Case when quadrant was found and theta needs correction: */
/* For Paper III: automatic correction of theta according to quadrant: */
  else if(good_q == -1) {
     me->theta += 180.;
     if(me->theta > 360.) me->theta -= 360.;
     strcpy(qflag,asterisk);
     }
/* Case when quadrant was not found (good_q = 0)*/
/* Check if measurements in 3rd Interferometric Catalog: */ 
  else {
    strcpy(qflag," ");
    status=1;
    if((obj[io]).WY != -1 && me->theta != NO_DATA) { 
           status = astrom_correct_theta_with_WDS_CHARA(&(obj[io]),me);
           }
// JLP2021: \rlap{!} removed
//    if(status) strcpy(qflag,exclam);
    }

  *q_error = (me->dquadrant == 1) ? '?' : ' ';
#ifdef Q_IN_NOTES
    sprintf(q_notes, "q=%d%c", me->quadrant, *q_error);
#else
    *q_notes = '\0';
#endif

// Oct2020: jj is no longer used for the same object (idem removed !)
  if(jj >= 0) {
// Conversion for LaTeX:
     astrom_preformat_wds_name(obj[io].wds, wds_name);
/*** Format for first line of an object */
/* October 2008: 3 decimals for the epoch */
  if(me->rho != NO_DATA) 
//    fprintf(fp_out,"%s & %s%s & %s & %.3f & %s & %d & %.3f & %.3f & %.1f%s & %.1f &  & %s %s\\\\\n", 
    fprintf(fp_out,"%s & %s%s & %.3f & %s & %d & %.3f & %.3f & %.1f%s & %.1f & %s %s\\\\\n", 
         wds_name, obj[io].discov_name, obj[io].comp_name, 
//         obj[io].ads, me->epoch, me->filter, 
         me->epoch, me->filter, 
         me->eyepiece, me->rho, me->drho, me->theta, qflag, me->dtheta, 
         me->notes, q_notes);
  else
//    fprintf(fp_out,"%s & %s%s & %s & %.3f & %s & %d & \\nodata & \\nodata & \\nodata & \\nodata &  & %s \\\\\n", 
    fprintf(fp_out,"%s & %s%s & %.3f & %s & %d & \\nodata & \\nodata & \\nodata & \\nodata &  & %s \\\\\n", 
         wds_name, obj[io].discov_name, obj[io].comp_name, 
//         obj[io].ads, me->epoch, me->filter, 
         me->epoch, me->filter, 
         me->eyepiece, me->notes);
  } /* EOF j == 0 */
// Oct2020: jj is no longer used for the same object (idem removed !)
/*** Format for subsequent lines of an object 
  else {
  if(me->rho != NO_DATA) 
//    fprintf(fp_out,"\\idem & \\idem & \\idem & %.3f & %s & %d & %.3f & %.3f & %.1f%s & %.1f & & %s %s\\\\\n", 
    fprintf(fp_out,"\\idem & \\idem & %.3f & %s & %d & %.3f & %.3f & %.1f%s & %.1f & %s %s\\\\\n", 
         me->epoch, me->filter, 
         me->eyepiece, me->rho, me->drho, me->theta, qflag, me->dtheta, 
         me->notes, q_notes);
  else
//    fprintf(fp_out,"\\idem & \\idem & \\idem & %.3f & %s & %d & \\nodata & \\nodata & \\nodata & \\nodata & & %s \\\\\n", 
    fprintf(fp_out,"\\idem & \\idem & %.3f & %s & %d & \\nodata & \\nodata & \\nodata & \\nodata & %s \\\\\n", 
         me->epoch, me->filter, me->eyepiece, me->notes);
  } // EOF j != 0 
*/
 nlines++;
// Oct2020: jj is no longer used for the same object (idem removed !)
 jj++;
 } /* EOF case not_flagged */
 } /* EOF loop on j */
} /* EOF loop on i */

 fprintf(fp_out,"& & & & & & & \\\\ \n");
 fprintf(fp_out,"\\hline \n");
 fprintf(fp_out,"\\end{tabular*} \n \n");
 fprintf(fp_out,"Note: In column 9, the exponent $^*$ indicates that the position angle\n");
 fprintf(fp_out,"$\\theta$ could be determined without the 180$^\\circ$ ambiguity.\n");
\
/* Problem for too big tables, so I only use tabular for big tables: */
 if(!tabular_only) fprintf(fp_out,"\\end{table*} \n");

return(0);
}
/*****************************************************************************
* Write the LateX array in a gili's format 
* without the filter, the ADS name 
*****************************************************************************/
int astrom_write_publi_table_gili(FILE *fp_out, int comments_wanted, 
                                  OBJECT *obj, int *index_obj, int nobj, 
                                  int tabular_only)
{
MEASURE *me, *me_next;
int good_q, nm, io, nlines, jj, jnext, nl_max, status;
char qflag[20], asterisk[20], exclam[20], q_error[1], wds_name[40]; 
char q_notes[20], dmag_string[16], nd_notes[16];
register int i, j;

/* Portrait: nl_max=55 (MNRAS for one page) */
/* Landscape: nl_max= 30 */
 if(comments_wanted)
   nl_max = 30;
 else
   nl_max = 53;

strcpy(asterisk,"\\rlap{$^*$}");
strcpy(exclam,"\\rlap{!}");

fprintf(fp_out,"\\def\\idem{''} \n");

nlines = -1;


for(i = 0; i < nobj; i++) {
  io = index_obj[i];
  nm = (obj[io]).nmeas;
/* New table header in publi_mode */
  if(nlines > nl_max) {
    fprintf(fp_out,"& & & & & & & & & \\\\ \n");
    fprintf(fp_out,"\\hline \n");
    fprintf(fp_out,"\\end{tabular*} \n");
/* Problem for too big tables, so I only use tabular for big tables: */
    if(!tabular_only) fprintf(fp_out,"\\end{table*} \n");
  }
  if((nlines == -1) || (nlines > nl_max)) {
/* JLP to obtain the same table number:
*/
/* Problem for too big tables, so I only use tabular for big tables: */
    if(!tabular_only) {
      if(nlines != -1) fprintf(fp_out,"\\addtocounter{table}{-1}\n");
      fprintf(fp_out,"\\begin{table*} \n");
      }
    fprintf(fp_out,"\\tabcolsep=1mm \n");
//    fprintf(fp_out,"\\begin{tabular*}{\\textwidth}{clrcccccrccl} \n");
    fprintf(fp_out,"\\begin{tabular*}{\\textwidth}{clccccrccl} \n");
    fprintf(fp_out,"\\hline \n");
    fprintf(fp_out,"& & & & & & & & & \\\\ \n");
    fprintf(fp_out,"WDS & Name & Epoch & Bin. & $\\rho$ & $\\sigma_\\rho$ & \\multicolumn{1}{c}{$\\theta$} & $\\sigma_\\theta$ & Dm & Notes \\\\ \n");
    fprintf(fp_out,"& &     &     & (\") & (\") & ($^\\circ$) & ($^\\circ$) & & \\\\ \n");
    fprintf(fp_out,"& & & & & & & & & \\\\ \n");
    fprintf(fp_out,"\\hline \n");
    fprintf(fp_out,"& & & & & & & & & \\\\ \n");
    nlines = 0;
  }

  jj = 0;
/************* Loop on measurements: **********************/
  for(j = 0; j < nm; j++) {
    me = &(obj[io]).meas[j];
// Conversion for LaTeX:
     astrom_preformat_wds_name(obj[io].wds, wds_name);
/* BOF case not_flagged */
    if(!me->flagged_out) {
      if(me->rho <= 0.1) {
// If too small, invalidate it:
         me->rho = NO_DATA;
// Unresolved case:
         fprintf(fp_out,"%s & %s%s & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata &  & Unres. \\\\\n", 
             wds_name, obj[io].discov_name, obj[io].comp_name, 
             me->epoch, me->eyepiece); //  me->notes);
      } else {
/* For DEBUG:
if(i < 4) {
  printf("i= %d io=%d nm=%d j=%d flag=%d\n", i, io, nm, j, me->flagged_out);
  printf(" WDS=%s rho=%.2f drho=%.2f theta=%.2f dtheta=%.2f\n", 
           obj[io].wds, me->rho, me->drho, me->theta, me->dtheta, me->eyepiece);
  }
*/

/* JLP2007: I check if among the next measurements there exists
* another one with the same epoch and filter
* than the current measurement: if it is the case, compute
* the mean of the two: */
    for(jnext = j+1; jnext < nm; jnext++) {
       me_next = &(obj[io]).meas[jnext];
       if(!me_next->flagged_out && (ABS(me->epoch - me_next->epoch) < 0.001) 
//         && !strncmp(me->filter,me_next->filter,3)
         && (me->eyepiece == me_next->eyepiece)
         && (me_next->rho != NO_DATA) && (me_next->theta != NO_DATA) ) {
/* Too many messages in case of error:
         printf("astrom_write_publi_gili/Warning same filter, epoch and eyepiece for object %s: I compute the mean (j=%d, jnext=%d)\n",
          obj[io].wds, j, jnext);
         printf("(filter=%s, epoch=%.4f (me_epoch=%.4f), eyepiece=%d)\n", me->filter, me_next->epoch,
                 me->epoch, me->eyepiece);
         printf("ABS=%f\n", ABS(me->epoch - me_next->epoch)); 
*/
         astrom_compute_mean_of_two_measures(obj, io, j, jnext);
         }
      }

  good_q = astrom_quadrant_is_consistent(me);
/* Case when quadrant was found and theta is OK: */
  if(good_q == 1) 
    strcpy(qflag,asterisk);
/* Case when quadrant was found and theta needs correction: */
/* For Paper III: automatic correction of theta according to quadrant: */
  else if(good_q == -1) {
     me->theta += 180.;
     if(me->theta > 360.) me->theta -= 360.;
     strcpy(qflag,asterisk);
     }
/* Case when quadrant was not found (good_q = 0)*/
/* Check if measurements in 3rd Interferometric Catalog: */ 
  else {
    strcpy(qflag," ");
    status=1;
    if((obj[io]).WY != -1 && me->theta != NO_DATA) { 
           status = astrom_correct_theta_with_WDS_CHARA(&(obj[io]),me);
           }
// JLP 2021 \rlap{!} removed
//    if(status) strcpy(qflag,exclam);
    }

  *q_error = (me->dquadrant == 1) ? '?' : ' ';
#ifdef Q_IN_NOTES
    sprintf(q_notes, "q=%d%c", me->quadrant, *q_error);
#else
    *q_notes = '\0';
#endif

// ND new double:
  if(me->is_new_double == 1) strcpy(nd_notes, "ND");
   else strcpy(nd_notes, "");

// dmag:
  if(me->dmag >= 0.) sprintf(dmag_string, "%.2f", me->dmag);
   else strcpy(dmag_string, "");

  if(jj >= 0) {
// Conversion for LateX
     astrom_preformat_wds_name(obj[io].wds, wds_name);
/*** Format for first line of an object */
/* October 2008: 3 decimals for the epoch */
  if(me->rho != NO_DATA) 
/******* WWWWWWWWWWWWW
WWWWWWWWWWWWWWWW **************/
    fprintf(fp_out,"%s & %s%s & %.3f & %d & %.3f & %.3f & %.1f%s & %.1f & %s & %s\\\\\n", 
         wds_name, obj[io].discov_name, obj[io].comp_name,
         me->epoch, me->eyepiece, me->rho, me->drho, 
         me->theta, qflag, me->dtheta, dmag_string, nd_notes); //, me->notes, q_notes);
  else {
// Doesn't seem to occur:
         fprintf(fp_out,"%s & %s%s & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata &  & Unres. \\\\\n", 
             wds_name, obj[io].discov_name, obj[io].comp_name, 
             me->epoch, me->eyepiece); //  me->notes);
     }
  } /* EOF j >= 0 */
/*** Format for subsequent lines of an object */
/***
  else {
  if(me->rho != NO_DATA) 
    fprintf(fp_out,"\\idem & \\idem & %.3f & %d & %.3f & %.3f & %.1f%s & %.1f & %s & %s \\\\\n", 
         me->epoch,  me->eyepiece, me->rho, me->drho, me->theta, 
         qflag, me->dtheta, dmag_string, nd_notes); //, me->notes, q_notes);
  else
    fprintf(fp_out,"\\idem & \\idem & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata & \nodata & \\\\\n", 
         me->epoch, me->eyepiece); //, me->notes);
    } // EOF j != 0 
***/
  } // End of resolved case
 nlines++;
 jj++;
 } /* EOF case not_flagged */
 } /* EOF loop on j */
} /* EOF loop on i */

 fprintf(fp_out,"& & & & & & & & & \\\\ \n");
 fprintf(fp_out,"\\hline \n");
 fprintf(fp_out,"\\end{tabular*} \n \n");
 fprintf(fp_out,"Note: In column 7, the exponent $^*$ indicates that the position angle\n");
 fprintf(fp_out,"$\\theta$ could be determined without the 180$^\\circ$ ambiguity.\n");
\
/* Problem for too big tables, so I only use tabular for big tables: */
 if(!tabular_only) fprintf(fp_out,"\\end{table*} \n");

return(0);
}
/*************************************************************************
* tex_calib_write_miniheader
*
*************************************************************************/
int tex_calib_write_miniheader(FILE *fp_out)
{
  fprintf(fp_out,"\\documentclass[10pt]{article} \n");
  fprintf(fp_out,"\\voffset=-3.cm \n");
  fprintf(fp_out,"\\hoffset=-3cm \n");
  fprintf(fp_out,"\\textwidth=26cm \n");
  fprintf(fp_out,"\\textheight=18cm \n");
  fprintf(fp_out,"\\parindent=0pt \n");
  fprintf(fp_out,"\\def\\nodata{---} \n");
  fprintf(fp_out,"\\begin{document} \n");
  fprintf(fp_out,"\\begin{table*} \n");
  fprintf(fp_out,"\\tabcolsep=1mm \n");
  fprintf(fp_out,"\\caption{Relative astrometry of binary stars} \n");
  fprintf(fp_out,"\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|} \n");
  fprintf(fp_out,"\\hline \n");
return (0);
}

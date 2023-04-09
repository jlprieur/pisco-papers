/*************************************************************************
* calib_csv_compare
* Program to compare the measurements of full_tab_calib_tex and Gili csv files
*
* Format of input files:
WDS & Name & Epoch & Bin. & $\rho$ & $\sigma_\rho$ & \multicolumn{1}{c}{$\theta$} & $\sigma_\theta$ & Dm & Notes \\
00014+3937  &  HLD60  & 2014.828 & 1 & 1.334 & 0.009 & 167.5  & 0.3 &  & \\
%
WDS & Name & Epoch & Bin. & $\rho$ & $\sigma_\rho$ & \multicolumn{1}{c}{$\theta$} & $\sigma_\theta$ & Dm & Notes & Orbit & {\scriptsize $\Delta \rho$(O-C)} & {\scriptsize $\Delta \theta$(O-C)} \\

00014+3937  &  HLD60  & 2014.828 & 1 & 1.322 & 0.007 & 167.5\rlap{$^*$} & 0.3 &  & & Izm2019 & 0.01 & 0.3 \\

* and Gili csv file:
,,,,,,,,,,,,,,
 a200\_a , 29/10/2014 ,  ,1,"8,04","0,12","-14,57","0,3"," EP=2014,8280 Q=1 dm=>3 ",,,"0,593352","75,43","76,42",

with filename in column 1, comments in col. 9, rho in col. 12, theta in col. 14
*
* JLP
* Version 14/01/2023
*************************************************************************/
#include "tex_calib_utils.h" 
#include "csv_utils.h" 
#include "latex_utils.h" 
#include "astrom_utils1.h" 
#include "astrom_utils2.h" 
#include "jlp_string.h"  // jlp_compact_string

static int calib_csv_compare_files(char *file_tex, char *file_csv, 
                                   char *fulltable_fileout, 
                                   char *discrep_cmp_out,
                                   int calib_fmt1, int calib_fmt2,
                                   double scale_mini);
/*****************
static int calib_tex_get_meas_from_file(char *file_tex, char *obj_name1,
                       double epoch1, 
                       int i_filename, int i_date, int i_filter,
                       int i_eyepiece, int i_rho, int i_drho,
                       int i_theta, int i_dtheta, int i_notes,
                       double *rho2, double *err_rho2, double *theta2,
                       double *err_theta2);
static int calib_tex_write_publi_cmp_table_gili(FILE *fp_out, 
                                  OBJECT *obj, int *index_obj, int nobj,
                                  int tabular_only, int quadrant_correction);
static int calib_tex_change_new_theta(double *theta2, double theta1); 
static int calib_tex_write_cmpfile(char *file_tex, char *file_csv, 
                                   char *discrep_cmp_out, OBJECT *obj1, 
                                   int nobj1);
static int calib_tex_write_header1(FILE *fp_out, int nlines, int tabular_only);
****************/

#define DEBUG
/*
*/

/*********************************************************************/
int main(int argc, char *argv[])
{
char file_tex[64], file_csv[64], fulltable_fileout[64], discrep_cmp_out[64];
int calib_fmt1, calib_fmt2;
double scale_mini;

/* If command line with "runs" */
if(argc == 7){
 if(*argv[5]) argc = 6;
 else if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
 else argc = 1;
 }

if(argc != 6)
  {
  printf(" Syntax: calib_csv_compare in_calib_tex_file in_csv_file fulltable_out discrepant_cmp_out calib_fmt1,calib_fmt2 scale_mini\n");
  printf(" Example: runs calib_csv_compare full_tab_calib.tex m1314_gili.csv table_5a_5b.tex discrep_cmp_5a_5b.txt 1,2 0.1\n");
  printf("argc=%d\n", argc);
  exit(-1);
  }
else
  {
  strcpy(file_tex,argv[1]);
  strcpy(file_csv,argv[2]);
  strcpy(fulltable_fileout,argv[3]);
  strcpy(discrep_cmp_out,argv[4]);
  sscanf(argv[5], "%d,%d", &calib_fmt1, &calib_fmt2);
  sscanf(argv[6], "%lf", &scale_mini);
  }

printf(" OK: file_tex=%s file_csv=%s \n", file_tex, file_csv);
printf(" OK: fulltable_fileout=%s discrep_cmp_out=%s\n", 
         fulltable_fileout, discrep_cmp_out);
printf(" OK: calib_fmt1 = %d calib_fmt2 = %d scale_mini=%f\n", 
         calib_fmt1, calib_fmt2, scale_mini);

/* Scan the file and add epoch from FITS autocorrelation files:
* (in "astrom_utils2.c")
*/ 
  calib_csv_compare_files(file_tex, file_csv, fulltable_fileout, 
                          discrep_cmp_out, calib_fmt1, calib_fmt2,
                          scale_mini); 

return(0);
}
/*************************************************************************
*
*
*************************************************************************/
static int calib_csv_compare_files(char *file_tex, char *file_csv, 
                               char *fulltable_fileout, char *discrep_cmp_out,
                               int calib_fmt1, int calib_fmt2, 
                               double scale_mini)
{
char b_in[NMAX], wds_name[40], discov_name[40], orbit__string[64];     
int status, istat, iline, compatible_meas31, orbit_status, verbose, eyep1;
int nunres3;
char *pc, *pc1, WDSName3[64], ObjectName3[64], filt1[32];
char ObjectName1[64];
double epoch1, rho1, err_rho1, theta1, err_theta1, dmag1;
double epoch3, rho3, err_rho3, theta3, err_theta3;
double max_drho, max_dtheta;
FILE *fp_in2, *fp_discrep_out, *fp_out;
OBJECT *obj1, *obj3;
int i, nobj1 = 0, nobj3 = 0, ncompat = 0;

if((fp_discrep_out = fopen(discrep_cmp_out,"w")) == NULL)
 {
 fprintf(stderr, "calib_csv_compare_files/Fatal error opening output file %s \n",discrep_cmp_out);
 exit(-1);
 }
if((fp_out = fopen(fulltable_fileout,"w")) == NULL)
 {
 fprintf(stderr, "calib_csv_compare_files/Fatal error opening output file %s \n",fulltable_fileout);
 exit(-1);
 }
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(fp_out,"%%%% compare_files: file_tex=%s and file_csv=%s\n", file_tex, file_csv);
fprintf(fp_out,"%%%% JLP / Version of 14/01/2023 \n");
fprintf(fp_discrep_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(fp_discrep_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(fp_discrep_out,"%%%% discrepant/compare file_tex=%s and file_csv=%s\n", file_tex, file_csv);
fprintf(fp_discrep_out,"%%%% JLP / Version of 14/01/2023 \n");

if((obj1 = (OBJECT *)malloc(NOBJ_MAX * sizeof(OBJECT))) == NULL) {
  printf("calib_tex_compare_files/Fatal error allocating memory space for OBJECT: nobj_max=%d\n",
          NOBJ_MAX);
  exit(-1);
  }
// Initialize the number of measurements to zero
for(i = 0; i < NOBJ_MAX; i++) (obj1[i]).nmeas = 0;

  printf("OK1 calling csv_read_gili_measures for %s\n", file_csv);
// Read the input Gili's csv meas. and create the catalog of meas. OBJECT obj1
 csv_read_gili_measures(file_csv, obj1, &nobj1, NOBJ_MAX, scale_mini);
#ifdef DEBUG
printf("calib_tex_compare_files: in_csv=%s nobj=%d scale_mini=%.3f\n", 
         file_csv, nobj1, scale_mini);
#endif

if((obj3 = (OBJECT *)malloc(2 * sizeof(OBJECT))) == NULL) {
  printf("calib_tex_compare_files/Fatal error allocating memory space for OBJECT: nobj_max=1 \n");
  exit(-1);
  }

// Scan the LaTeX table of measurements
if((fp_in2 = fopen(file_tex,"r")) == NULL) {
  fprintf(stderr, " Fatal error opening input tex file %s \n", file_tex);
  return(-1);
  }

wds_name[0] = '\0';
discov_name[0] = '\0';
iline = 0;
while(!feof(fp_in2))
{
/* Maximum length for a line will be 170 characters: */
  if(fgets(b_in,170,fp_in2))
  {
  b_in[169] = '\0';
/* NEW/2009: I remove ^M (Carriage Return) if present: */
  pc = b_in;
  while(*pc) {
  if(*pc == '\r') *pc = ' ';
   pc++;
   }
  iline++;
  if((iline % 1000) == 0)printf("iline=%d\n", iline);
// Copy input line to output file:
  fprintf(fp_out, "%s", b_in);
  nunres3 = 0;
  nobj3 = 0;
  obj3[0].nmeas = 0;
  status = tex_calib_read_measures_from_line(b_in, obj3, &nobj3, &nunres3,
                                             calib_fmt1);
  if(status == 0 ) {
    strcpy(WDSName3, (obj3[0]).wds);
    strcpy(ObjectName3, (obj3[0]).discov_name);
    epoch3 = (obj3[0].meas[0]).epoch;
    rho3 = (obj3[0].meas[0]).rho;
    err_rho3 = (obj3[0].meas[0]).drho;
    theta3 = (obj3[0].meas[0]).theta;
    verbose = 0;
    orbit_status = latex_get_column_item(b_in, orbit__string, 11, verbose);
    err_theta3 = (obj3[0].meas[0]).dtheta;
    if((rho3 > 0 && rho3 != NO_DATA) && (theta3 != NO_DATA)) {
// Look for the measurements of the same object in the catalog obj1 
// using ObjectName3 and epoch3:
      istat = tex_calib_get_meas_from_csv_object(obj1, nobj1, ObjectName3, 
                                       epoch3, ObjectName1, &epoch1, &eyep1,
                                       &rho1, &err_rho1, &theta1, &err_theta1, 
                                       &dmag1);
      if(istat == 0) {
        max_drho = 0.3;
        max_dtheta = 100.;
        tex_calib_compare(rho3, err_rho3, theta3, err_theta3, 
                          rho1, err_rho1, theta1, err_theta1,
                          max_drho, max_dtheta,
                          &compatible_meas31); 
/*
printf("ObjectName3=%s epoch3=%f\n",
        ObjectName3, epoch3);
printf("rho3=%f err_rho3=%f theta3=%f err_theta3=%f \n",
           rho3, err_rho3, theta3, err_theta3);
printf("ObjectName1=%s epoch1=%f eyep1=%d dmag1=%f \n",
        ObjectName1, epoch1, eyep1, dmag1);
printf("rho1=%f err_rho1=%f theta1=%f err_theta1=%f \n",
           rho1, err_rho1, theta1, err_theta1);
*/

        jlp_compact_string(WDSName3, 64);
        jlp_compact_string(ObjectName1, 64);
        fprintf(fp_out, "888 %s & %s & %.3f & %d & %.3f & %.3f & %.1f & %.1f & %.2f & \\\\ \n",
               WDSName3, ObjectName1, epoch1, eyep1, rho1, err_rho1, 
               theta1, err_theta1, dmag1);
        if(compatible_meas31 == 1) {
           ncompat++;
// Do not output objects with residuals since I assume they have been checked
          } else if((rho1 != NO_DATA) && (orbit_status != 0)) {
// Save discrepant measurements:
           fprintf(fp_discrep_out, "%s", b_in);
           fprintf(fp_discrep_out, "888 %s & %s & %.3f & %d & %.3f & %.3f & %.1f & %.1f & %.2f & \\\\ \n",
               WDSName3, ObjectName1, epoch1, eyep1, rho1, err_rho1, 
               theta1, err_theta1, dmag1);
          }
       } // EOF istat == 0
      } // EOF rho3 != NO_DATA
    } // EOF if(status==0)
/*
int tex_calib_get_meas_from_object(OBJECT *obj1, int *nobj1,
          char *WDSName0, char *ObjectName0, char *WDSName1, char *ObjectName1,
          double *epoch1, char *filt1, int *eyep1, double *rho1, double *err_rho1, double *theta1, 
          double *err_theta1);
*/

  } /* EOF if fgets() */
} /* EOF while loop */
printf("Output: nlines=%d nobj1=%d, ncompat=%d\n", 
        iline, nobj1, ncompat);

#ifdef DEBUG1
if(nobj1 > 1) {
printf("Returned by calib_tex_read_measures: nobj1 = %d nmeas=%d rho=%.2f theta=%.2f\n", 
        nobj1, obj1[nobj1-1].nmeas, obj1[nobj1-1].meas[0].rho, 
        obj1[nobj1-1].meas[0].theta);
}
#endif

free(obj1);
free(obj3);
fclose(fp_in2);
fclose(fp_out);
fclose(fp_discrep_out);
return(0);
}

/********************************************************************
OLD VERSION: 

   status = calib_tex_check_measure(b_data, i_eyepiece, i_rho, i_drho,
                                       i_theta, i_dtheta);
#ifdef DEBUG1
printf("calib_tex_read_measures/Adding new measurement for object #i_obj=%d nm=%d (nobj=%d, calib_tex_check_measure: status=%d)\n",
         i_obj, (obj1[i_obj]).nmeas, *nobj1, status);
#endif
          if(status == 0) {
           comments_wanted = 0;
           calib_tex_add_new_measure(b_data, obj, i_obj, i_filename, i_date,
                           i_filter, i_eyepiece, i_rho, i_drho, i_theta,
                           i_dtheta, i_notes, comments_wanted);
// theta is always in [0, 360.]:
          if(obj[i_obj].meas[0].theta < 0.) 
                obj[i_obj].meas[0].theta += 360.;
          else if(obj[i_obj].meas[0].theta > 360.) 
                obj[i_obj].meas[0].theta -= 360.;
#ifdef DEBUG
           printf("calib_tex_read_meas/fname=%s nmeas=%d rho=%.2f+/-%.2f theta=%.2f+/-%.2f\n", 
                  obj[i_obj].name, obj[i_obj].meas[0].rho, 
                  obj[i_obj].meas[0].drho, 
                  (obj[i_obj]).nmeas, obj[i_obj].meas[0].theta,
                  obj[i_obj].meas[0].dtheta);
#endif
           strcpy(obj_name1, (obj[i_obj]).name);
           epoch1 = (obj[i_obj]).meas[0].epoch;
// Return 0 when obj_name1 was found in file_csv
           status = calib_tex_get_meas_from_file(file_csv, obj_name1,
                       epoch1, i_filename, i_date, i_filter,
                       i_eyepiece, i_rho, i_drho, i_theta, i_dtheta, i_notes,
                       &rho2, &err_rho2, &theta2, &err_theta2);
           if(status == 0) {
#ifdef DEBUG
            printf("calib_tex_get_meas/%s in %s: rho=%.2f+/-%.2f theta=%.2f+/-%.2f\n",
                  obj_name1, file_csv, rho2, err_rho2, theta2, err_theta2);
#endif
             obj[i_obj].cmp_meas[0].rho = rho2; 
             obj[i_obj].cmp_meas[0].drho = err_rho2;
// Change the theta value if needed:
             if((obj[i_obj].meas[0].theta != NO_DATA) && (theta2 != NO_DATA))
                calib_tex_change_new_theta(&theta2, obj[i_obj].meas[0].theta); 
             obj[i_obj].cmp_meas[0].theta = theta2;
             obj[i_obj].cmp_meas[0].dtheta = err_theta2;
            } else {
// Set the measurements to -NO_DATA if the file was not found
             obj[i_obj].cmp_meas[0].rho = -NO_DATA; 
             obj[i_obj].cmp_meas[0].drho = -NO_DATA;
             obj[i_obj].cmp_meas[0].theta = -NO_DATA;
             obj[i_obj].cmp_meas[0].dtheta = -NO_DATA;
#ifdef DEBUG
            printf("calib_tex_read_meas/Warning: object %s not found in %s\n", 
                   obj_name1, file_csv);
#endif
            }
            } // case status from check_measure = 0
         }
****************/

#ifdef TTT
/*********************************************************************
* Return 0 when obj_name1 was found in file_csv
*
*********************************************************************/
static int calib_tex_get_meas_from_file(char *file_csv, char *obj_name1,
                       double epoch1, 
                       int i_filename, int i_date, int i_filter,
                       int i_eyepiece, int i_rho, int i_drho,
                       int i_theta, int i_dtheta, int i_notes,
                       double *rho2, double *err_rho2, double *theta2,
                       double *err_theta2)
{
char b_in[NMAX], b_data[NMAX], filename2[64];
char *pc, filter2[64], notes2[64];
double epoch2; 
int eyepiece2, status, status0 = -1;
FILE *fp_in2;

// Remove the extension ".fits" for instance:
jlp_remove_ext_string(obj_name1, 64);
jlp_compact_string(obj_name1, 64);
#ifdef DEBUG1
printf("\n calib_tex_get_meas_from_file/looking for fname=>%s< in %s\n",
         obj_name1, file_csv);
#endif

if((fp_in2 = fopen(file_csv,"r")) == NULL) {
  fprintf(stderr, " Fatal error opening input csv file %s \n", file_csv);
  return(-1);
  }

while(!feof(fp_in2))
{
/* Maximum length for a line will be 170 characters: */
  if(fgets(b_in,170,fp_in2))
  {
  if(b_in[0] != '%') {
     b_in[169] = '\0';
/* NEW/2009: I remove ^M (Carriage Return) if present: */
    pc = b_in;
    while(*pc) {
    if(*pc == '\r') *pc = ' ';
    pc++;
    }
    strcpy(b_data, b_in);
    status = calib_tex_check_measure(b_data, i_eyepiece, i_rho, i_drho, i_theta, 
                                i_dtheta);
    if(status == 0) {
#ifdef DEBUG1
printf("calib_tex_get_meas_from_file/calib_tex_check_measure: b_data=%s status=%d)\n",
         b_data, status);
#endif
          calib_tex_read_new_measure(b_data, i_filename, i_date, i_filter, 
                                  i_eyepiece, i_rho, i_drho, i_theta,
                                  i_dtheta, i_notes, filename2, notes2, 
                                  &epoch2, filter2, &eyepiece2,
                                  rho2, err_rho2, theta2, err_theta2);
// Remove the extension ".fits" for instance:
          jlp_remove_ext_string(filename2, 64);
          jlp_compact_string(filename2, 64);
          if(!strcmp(obj_name1, filename2)) {
#ifdef DEBUG1
            printf("calib_tex_read_measures/in %s: fname=%s rho=%.2f+/-%.2f theta=%.2f+/-%.2f\n",
                  file_csv, filename2, *rho2, *err_rho2, *theta2, *err_theta2);
#endif
            status0 = 0;
            break;
            }

       }
     } // EOF if b_in[0] != '%' 
   } // EOF if fgets... 

} // EOF while
fclose(fp_in2);
return(status0);
}
/*****************************************************************************
* Write the LateX array in a gili's format
* without the filter, the ADS name and the orbit
*****************************************************************************/
static int calib_tex_write_publi_cmp_table_gili(FILE *fp_out, 
                                  OBJECT *obj, int *index_obj, int nobj,
                                  int tabular_only, int quadrant_correction)
{
MEASURE *me, *cmp_me, *me_next;
int good_q, nm, io, nlines, jj, jnext, nl_max, status;
char qflag[20], asterisk[20], exclam[20], q_error[1], wds_name[40];
char q_notes[20], dmag_string[16];
register int i, j;

// Max number of objects per page:
   nl_max = 26;

strcpy(asterisk,"\\rlap{$^*$}");
strcpy(exclam,"\\rlap{!}");

fprintf(fp_out,"\\def\\idem{''} \n");

nlines = -1;

for(i = 0; i < nobj; i++) {
  io = index_obj[i];
  nm = (obj[io]).nmeas;
/* New table header in publi_mode */
  if(nlines > nl_max) {
//    fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
    fprintf(fp_out,"& & & & & & & & & \\\\ \n");
    fprintf(fp_out,"\\hline \n");
    fprintf(fp_out,"\\end{tabular*} \n");
/* Problem for too big tables, so I only use tabular for big tables: */
    if(!tabular_only) fprintf(fp_out,"\\end{table*} \n");
  }
  if((nlines == -1) || (nlines > nl_max)) {
    calib_tex_write_header1(fp_out, nlines, tabular_only);
    nlines = 0;
  }

  jj = 0;
/************* Loop on measurements: **********************/
  for(j = 0; j < nm; j++) {
    me = &(obj[io]).meas[j];
    cmp_me = &(obj[io]).cmp_meas[j];
/* BOF case not_flagged */
    if(!me->flagged_out) {
      if(me->rho <= 0.1) {
// Set to no_data if rho is too small
         me->rho = NO_DATA;
// Unresolved case:
         fprintf(fp_out,"%s & %s & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata &  & Unres. \\\\\n",
             wds_name, obj[io].name, me->epoch, me->eyepiece); //  me->notes);
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
         printf("calib_tex_write_publi_gili/Warning same filter, epoch and eyepiece for object %s: I compute the mean (j=%d, jnext=%d)\n",
          obj[io].wds, j, jnext);
         printf("(filter=%s, epoch=%.4f (me_epoch=%.4f), eyepiece=%d)\n", me->filter, me_next->epoch,
                 me->epoch, me->eyepiece);
         printf("ABS=%f\n", ABS(me->epoch - me_next->epoch));
*/
         calib_tex_compute_mean_of_two_measures(obj, io, j, jnext);
         }
      }

// Check if theta value (me->theta) is compatible with the quadrant value (me->quadrant):
  good_q = calib_tex_quadrant_is_consistent(me);
// printf("Before quadrant correction: me->theta=%.2f cmp_me->theta=%.2f (quadrant=%d) good_q=%d\n", 
//        me->theta, cmp_me->theta, me->quadrant, good_q);

/* Case when quadrant was found and theta is OK: */
  if(good_q == 1)
    strcpy(qflag,asterisk);
/* Case when quadrant was found and theta needs correction: */
/* Automatic correction of theta according to quadrant: */
  else if(good_q == -1) {
    if(quadrant_correction == 1) {
      if(me->theta != NO_DATA) {
        me->theta += 180.;
        if(me->theta > 360.) me->theta -= 360.;
        if((cmp_me->theta != NO_DATA) && (cmp_me->theta != -NO_DATA)) {
          cmp_me->theta += 180.;
          if(cmp_me->theta > 360.) cmp_me->theta -= 360.;
          }
        }
      } // EOF quadrant_correction
       strcpy(qflag,asterisk);
/* Case when quadrant was not found (good_q = 0)*/
/* Check if measurements in 3rd Interferometric Catalog: */
  } else {
    strcpy(qflag," ");
    status=1;
      if(quadrant_correction == 1) {
        if((obj[io]).WY != -1 && me->theta != NO_DATA) {
           status = calib_tex_correct_theta_with_WDS_CHARA(&(obj[io]),me);
           if((cmp_me->theta != NO_DATA) && (cmp_me->theta != -NO_DATA))
              calib_tex_correct_theta_with_WDS_CHARA(&(obj[io]), cmp_me);
        }
        if(status) strcpy(qflag,exclam);
      } // EOF quadrant_correction
    }
// printf("After quadrant correction: me->theta=%.2f cmp->theta=%.2f (quadrant=%d)\n", 
//    me->theta, cmp_me->theta, me->quadrant);

  *q_error = (me->dquadrant == 1) ? '?' : ' ';
#ifdef Q_IN_NOTES
    sprintf(q_notes, "q=%d%c", me->quadrant, *q_error);
#else
    *q_notes = '\0';
#endif

// dmag:
  if(me->dmag >= 0.) sprintf(dmag_string, "%.2f", me->dmag);
   else strcpy(dmag_string, "");

  if(jj == 0) {
     calib_tex_preformat_wds_name(obj[io].wds, wds_name);
/*** Format for first line of an object */
/* October 2008: 3 decimals for the epoch */
  if(me->rho != NO_DATA)
    fprintf(fp_out,"%s & %s & %.3f & %d & %.3f & %.3f & %.1f%s & %.1f & %s &\\\\\n",
         wds_name, obj[io].name, me->epoch, me->eyepiece, me->rho, me->drho,
         me->theta, qflag, me->dtheta, dmag_string); //, me->notes, q_notes);
  else
    fprintf(fp_out,"%s & %s & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata & \\nodata &  \\\\\n",
         wds_name, obj[io].name, me->epoch, me->eyepiece); //, me->notes);
// File was not found when cmp_me->rho == -NO_DATA
  if(cmp_me->rho != -NO_DATA) {
  if(cmp_me->rho != NO_DATA)
    fprintf(fp_out,"cmp:%s & %s & %.3f & %d & %.3f & %.3f & %.1f%s & %.1f & %s &\\\\\n",
         wds_name, obj[io].name, me->epoch, me->eyepiece, cmp_me->rho, cmp_me->drho,
         cmp_me->theta, qflag, cmp_me->dtheta, dmag_string); //, me->notes, q_notes);
  else
    fprintf(fp_out,"cmp:%s & %s & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata & \\nodata &  \\\\\n",
         wds_name, obj[io].name, me->epoch, me->eyepiece); //, me->notes);
  }
  } /* EOF j == 0 */
/*** Format for subsequent lines of an object */
  else {
/************ BOF with idem 
  if(me->rho != NO_DATA)
    fprintf(fp_out,"\\idem & \\idem & %.3f & %d & %.3f & %.3f & %.1f%s & %.1f & %s & \\\\\n",
         me->epoch,  me->eyepiece, me->rho, me->drho, me->theta,
         qflag, me->dtheta, dmag_string); //, me->notes, q_notes);
  else
    fprintf(fp_out,"\\idem & \\idem & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata & \nodata & \\\\\n",
         me->epoch, me->eyepiece); //, me->notes);
******** EOF with idem *****************/
/************ BOF Without idem  ****/
  if(me->rho != NO_DATA)
    fprintf(fp_out,"%s & %s & %.3f & %d & %.3f & %.3f & %.1f%s & %.1f & %s &\\\\\n",
         wds_name, obj[io].name, me->epoch, me->eyepiece, me->rho, me->drho,
         me->theta, qflag, me->dtheta, dmag_string); //, me->notes, q_notes);
  else
    fprintf(fp_out,"%s & %s & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata & \\nodata &  \\\\\n",
         wds_name, obj[io].name, me->epoch, me->eyepiece); //, me->notes);
/************ EOF without idem  ****/
    } // End of resolved case
  } /* EOF j != 0 */
 nlines++;
 jj++;
 } /* EOF case not_flagged */
 } /* EOF loop on j */
} /* EOF loop on i */

// fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
 fprintf(fp_out,"& & & & & & & & & \\\\ \n");
 fprintf(fp_out,"\\hline \n");
 fprintf(fp_out,"\\end{tabular*} \n \n");
// fprintf(fp_out,"Note: In column 9, the exponent $^*$ indicates that the position angle\n");
 fprintf(fp_out,"Note: In column 7, the exponent $^*$ indicates that the position angle\n");
 fprintf(fp_out,"$\\theta$ could be determined without the 180$^\\circ$ ambiguity.\n");
\
/* Problem for too big tables, so I only use tabular for big tables: */
 if(!tabular_only) fprintf(fp_out,"\\end{table*} \n");

return(0);
}
/**********************************************************************
* Change theta2 to theta2 + 180 or to theta2 + 360 if needed
*
* INPUT:
*  theta1 (in degrees) 
*  NB: theta1 is always in [0, 360.]
*
* INPUT/OUTPUT:
*  theta2 (in degrees)
**********************************************************************/
static int calib_tex_change_new_theta(double *theta2, double theta1) 
{
double thet0 = 10., error0;

/// printf("UUU/Before: theta1=%.2f theta2=%.2f\n", theta1, *theta2);

// NB theta1 is always in [0, 360.]
// so theta2 is always in [0, 360.]
if(*theta2 < 0.) *theta2 += 360.;
else if(*theta2 > 360.) *theta2 -= 360.;

/// printf("UUU/Before666: theta1=%.2f theta2=%.2f\n", theta1, *theta2);

// First quadrant: 
if((theta1 > 0.) && (theta1 < 90.)) { 
// Correct from 3rd quadrant to first quadrant:
   if(*theta2 > 180. - thet0)  *theta2 -= 180.;
// Second quadrant: 
} else if((theta1 > 90.) && (theta1 < 180.)) {
// Correct from 4th quadrant to 2nd quadrant:
   if(*theta2 > 270. - thet0)  *theta2 -= 180.;
// Third quadrant: 
} else if((theta1 > 180.) && (theta1 < 270.)) {
// Correct from 1st quadrant to 3rd quadrant:
   if(*theta2 < 90. + thet0)  *theta2 += 180.;
// Fourth quadrant: 
} else if((theta1 > 270.) && (theta1 < 360.)) { 
// Correct from 2nd quadrant to 4th quadrant:
   if(*theta2 < 180. + thet0)  *theta2 += 180.;
}
// Look for possible other errors
   error0 = theta1 - *theta2;
   if (error0 > 180. - thet0) {
     *theta2 += 180.;
     } else if (error0 < -180. + thet0) {
     *theta2 -= 180.;
     }
   error0 = *theta2 - theta1;
   if (error0 > 180. - thet0) {
     *theta2 -= 180.;
     } else if (error0 < -180. + thet0) {
     *theta2 += 180.;
     }
/// printf("UUU/After: theta1=%.2f theta2=%.2f\n", theta1, *theta2);
return(0);
}
/****************************************************************
* Compare the two sets of measurements and output
* a file containing the LaTex table for which the compared measures differ
*
****************************************************************/
static int calib_tex_write_cmpfile(char *file_tex, char *file_csv, 
                                 char *discrep_cmp_out, OBJECT *obj1, int nobj1)
{
FILE *fp_cmp_out;
OBJECT *obj2;
MEASURE *me1, *cmp_me1, *me2, *cmp_me2;
double drho_max = 3;
double dtheta_max = 8;
int *index_obj2, tabular_only, nm1;
int i, j, k, nobj2 = 0;
int quadrant_correction;

if((obj2 = (OBJECT *)malloc(NOBJ_MAX * sizeof(OBJECT))) == NULL) {
  printf("calib_tex_compare_files/Fatal error allocating memory space for OBJECT: nobj_max=%d\n",
          NOBJ_MAX);
  exit(-1);
  }
// Initialize the number of measurements to zero
for(i = 0; i < NOBJ_MAX; i++) (obj2[i]).nmeas = 0;

if((index_obj2 = (int *)malloc((NOBJ_MAX) * sizeof(int))) == NULL) {
  printf("calib_tex_compare_files/Fatal error allocating memory space for index_obj: nobj_max=%d\n",
          NOBJ_MAX);
  exit(-1);
  }

k = 0;
for(i = 0; i < nobj1; i++) {
  nm1 = (obj1[i]).nmeas;
/************* Loop on measurements: **********************/
  for(j = 0; j < nm1; j++) {
    me1 = &(obj1[i]).meas[j];
    cmp_me1 = &(obj1[i]).cmp_meas[j];
/* BOF case not_flagged */
      if(!me1->flagged_out 
// rho=NO_DATA when not resolved...
      && (me1->rho != NO_DATA) && (me1->theta != NO_DATA) 
// rho=NO_DATA when not resolved...
      && (cmp_me1->rho != NO_DATA) && (cmp_me1->theta != NO_DATA)
// rho=-NO_DATA when file was not found...
      && (cmp_me1->rho != -NO_DATA) && (cmp_me1->theta != -NO_DATA)) {
         if((ABS(me1->rho - cmp_me1->rho) > drho_max)
            || (ABS(me1->theta - cmp_me1->theta) > dtheta_max)) {
         printf("name=%s rho=%.2f theta=%.2f Drho=%.2f Dtheta=%.2f\n",
               (obj1[i]).name, me1->rho, me1->theta,
               ABS(me1->rho - cmp_me1->rho),
               ABS(me1->theta - cmp_me1->theta));
             strcpy((obj2[k]).wds, (obj1[i]).wds);
             strcpy((obj2[k]).name, (obj1[i]).name);
             (obj2[k]).nmeas = (obj1[i]).nmeas;
             (obj2[k]).WR = (obj1[i]).WR;
             (obj2[k]).WY= (obj1[i]).WY;
             me2 = &(obj2[k]).meas[j];
             cmp_me2 = &(obj2[k]).cmp_meas[j];
             strcpy(me2->filename, me1->filename); 
             strcpy(me2->date, me1->date); 
             strcpy(me2->notes, me1->notes); 
             me2->eyepiece = me1->eyepiece; 
             me2->epoch = me1->epoch; 
             me2->quadrant = me1->quadrant; 
             me2->dquadrant = me1->dquadrant; 
             me2->dmag = me1->dmag; 
             me2->ddmag = me1->ddmag; 
             me2->rho = me1->rho; 
             me2->drho = me1->drho; 
             me2->theta = me1->theta; 
             me2->dtheta = me1->dtheta; 
             cmp_me2->rho = cmp_me1->rho; 
             cmp_me2->drho = cmp_me1->drho; 
             cmp_me2->theta = cmp_me1->theta; 
             cmp_me2->dtheta = cmp_me1->dtheta; 
             k++;
             } // > drho_max
      } // not flagged out 
    } // loop on j nm1
  } // loop on i nobj
nobj2 = k;
if(nobj2 > 1) {
printf("Returned by calib_tex_read_measures: nobj2=%d nmeas=%d rho=%.2f theta=%.2f\n",
        nobj2, obj2[nobj2-1].nmeas, obj2[nobj2-1].meas[0].rho,
        obj2[nobj2-1].meas[0].theta);
}
#ifdef DEBUG1
#endif


/* Sort the objects according to their name: */
calib_tex_name_sort_objects(obj2, index_obj2, nobj2);
if((fp_cmp_out = fopen(discrep_cmp_out,"w")) == NULL) {
  fprintf(stderr, "calib_tex_create_cmpfile/Error opening output cmp file %s \n",
        discrep_cmp_out);
  return(-1);
  }

fprintf(fp_cmp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(fp_cmp_out,"%%%% calib_tex_compare tex file=%s and csv file=%s\n", file_tex, file_csv);
fprintf(fp_cmp_out,"%%%% Measurements that differ nobj2=%d (nobj1=%d) with Drho=%.2f and Dtheta=%.2f\n", 
        nobj2, nobj1, drho_max, dtheta_max);
fprintf(fp_cmp_out,"%%%% JLP / Version of 31/03/2020 \n");
fprintf(fp_cmp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

/* For big tables, should set this parameter to 1: */
tabular_only = 0;
// NO correction for the quadrant (already done)
quadrant_correction = 0;
calib_tex_write_publi_cmp_table_gili(fp_cmp_out, obj2, index_obj2, nobj2,
                         tabular_only, quadrant_correction);

fclose(fp_cmp_out);
free(index_obj2);
free(obj2);
return(0);
}
/**************************************************************************
* Write the LateX header of the table
***************************************************************************/
static int calib_tex_write_header1(FILE *fp_out, int nlines, int tabular_only)
{
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
//    fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
    fprintf(fp_out,"& & & & & & & & & \\\\ \n");
//    fprintf(fp_out,"WDS & Name & ADS & Epoch & Fil. & Eyep. & $\\rho$ & $\\sigma_\\rho$ & \\multicolumn{1}{c}{$\\theta$} & $\\sigma_\\theta$ & Orb. & Notes \\\\ \n");
    fprintf(fp_out,"WDS & Name & Epoch & Bin. & $\\rho$ & $\\sigma_\\rho$ & \\multicolumn{1}{c}{$\\theta$} & $\\sigma_\\theta$ & Dm & Notes \\\\ \n");
//    fprintf(fp_out,"& &       &       &      &     & (mm)& (\") & (\") & ($^\\circ$) & ($^\\circ$) & \\\\ \n");
    fprintf(fp_out,"& &     &     & (\") & (\") & ($^\\circ$) & ($^\\circ$) & & \\\\ \n");
//    fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
    fprintf(fp_out,"& & & & & & & & & \\\\ \n");
    fprintf(fp_out,"\\hline \n");
//    fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
    fprintf(fp_out,"& & & & & & & & & \\\\ \n");

return(0);
}
#endif

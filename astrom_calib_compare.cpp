/*************************************************************************
* astrom_compare
* Program to compare the measurements of two astrom files
*
* Format of input files:
& 120109\_ads684ab\_Rd\_8\_a & 12/01/2009 & R & 20 & 11.71 & 0.10 & -19.37 & 0.9 & Q=? EP=2009.0852 \\
*
Possible keywords in the notes:
EP=2004.1323 (epoch)
Q=2         (Quadrant with restricted triple-correlation)
LQ=3        (Quadrant with long integration)
*
* JLP
* Version 31/03/2020
*************************************************************************/
#include "astrom_utils1.h" 
#include "astrom_utils2.h" 
#include "jlp_string.h"  // jlp_compact_string

static int astrom_compare_meas(char *filein1, char *filein2, 
                               char *fulltable_fileout, char *cmp_fileout,
                               char *unres_cmp_fileout);
static int astrom_compare_files(FILE *fp_in1, FILE *fp_out, 
                          char *filein1, char *filein2, char *cmp_fileout, 
                          char *unres_cmp_fileout,
                          int i_filename, int i_date, int i_filter,
                          int i_eyepiece, int i_rho, int i_drho,
                          int i_theta, int i_dtheta, int i_notes);
static int astrom_read_measures_cmp(FILE *fp_in, char *filein2,
                         OBJECT *obj, int *nobj, int i_filename, int i_date,
                         int i_filter, int i_eyepiece, int i_rho,
                         int i_drho, int i_theta, int i_dtheta, int i_notes,
                         int with_wds_data);
static int astrom_get_meas_from_file(char *filein2, char *obj_name1,
                       double epoch1, 
                       int i_filename, int i_date, int i_filter,
                       int i_eyepiece, int i_rho, int i_drho,
                       int i_theta, int i_dtheta, int i_notes,
                       double *rho2, double *err_rho2, double *theta2,
                       double *err_theta2);
static int astrom_write_publi_cmp_table_gili(FILE *fp_out, 
                                  OBJECT *obj, int *index_obj, int nobj,
                                  int tabular_only, int quadrant_correction);
static int astrom_change_new_theta(double *theta2, double theta1); 
static int astrom_write_cmpfile(char *filein1, char *filein2, 
                                  char *cmp_fileout, OBJECT *obj1, int nobj1);
static int astrom_write_unrescmpfile(char *filein1, char *filein2, 
                                     char *unres_cmp_fileout, OBJECT *obj1, 
                                     int nobj1);
static int astrom_write_header1(FILE *fp_out, int nlines, int tabular_only);

/*
#define DEBUG
*/

/*********************************************************************/
int main(int argc, char *argv[])
{
char filein1[64], filein2[64], fulltable_fileout[64], cmp_fileout[64];
char unres_cmp_fileout[64];

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
  printf(" Syntax: astrom_compare in_astrom_file1 in_astrom_file2 fulltable_outfile compared_outfile unres_cmp_outfile\n");
  printf(" Example: runs astrom_compare astrom05a.tex astrom05b.tex table_5a_5b.tex cmp_5a_5b.txt unres_cmp_5b.txt\n");
  exit(-1);
  }
else
  {
  strcpy(filein1,argv[1]);
  strcpy(filein2,argv[2]);
  strcpy(fulltable_fileout,argv[3]);
  strcpy(cmp_fileout,argv[4]);
  strcpy(unres_cmp_fileout,argv[5]);
  }

printf(" OK: filein1=%s filein2=%s \n", filein1, filein2);
printf(" OK: fulltable_fileout=%s cmp_fileout=%s\n", 
         fulltable_fileout, cmp_fileout);
printf(" OK: unres_cmp_fileout=%s\n", unres_cmp_fileout);

/* Scan the file and add epoch from FITS autocorrelation files:
* (in "astrom_utils2.c")
*/ 
  astrom_compare_meas(filein1, filein2, fulltable_fileout, cmp_fileout,
                      unres_cmp_fileout); 

return(0);
}
/*************************************************************************
*
*
*************************************************************************/
static int astrom_compare_meas(char *filein1, char *filein2, 
                               char *fulltable_fileout, char *cmp_fileout,
                               char *unres_cmp_fileout)
{
FILE *fp_in1, *fp_out;
int i_filename, i_date, i_filter, i_eyepiece, i_rho, i_drho, i_theta, i_dtheta;
int i_notes, publi_mode, input_with_header;


if((fp_in1 = fopen(filein1,"r")) == NULL) {
  fprintf(stderr, " Fatal error opening input file1 %s \n", filein1);
  return(-1);
  }

if((fp_out = fopen(fulltable_fileout,"w")) == NULL)
{
fprintf(stderr, "astrom_compare_meas/Fatal error opening output file %s \n",fulltable_fileout);
exit(-1);
}
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(fp_out,"%%%% astrom_compare file1=%s and file2=%s\n", filein1, filein2);
fprintf(fp_out,"%%%% JLP / Version of 31/03/2020 \n");
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

/* Scan the file and make the conversion: */
/* date in column 3
*  eyepiece in column 5
*  rho in column 6
*  drho in column 7
*  theta in column 8
*  dtheta in column 9
*/
  i_filename = 2;
  i_date = 3;
  i_filter = 4;
  i_eyepiece = 5;
  i_rho = 6;
  i_drho = 7;
  i_theta = 8;
  i_dtheta = 9;
  i_notes = 10;
/* Version "publi", to generate a laTeX array formatted for publication:
* sort out the objects, compute mean values, perform calibration, etc.
*/
  astrom_compare_files(fp_in1, fp_out, filein1, filein2, cmp_fileout, 
                       unres_cmp_fileout, i_filename, i_date, i_filter, 
                       i_eyepiece, i_rho, i_drho, i_theta, i_dtheta, i_notes); 
fclose(fp_in1);
fclose(fp_out);
return(0);
}
/*************************************************************************
*
* INPUT:
* i_filename: column nber of the filename used for this measurement
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
*
*************************************************************************/
static int astrom_compare_files(FILE *fp_in1, FILE *fp_out, 
                          char *filein1, char *filein2, char *cmp_fileout, 
                          char *unres_cmp_fileout,
                          int i_filename, int i_date, int i_filter,
                          int i_eyepiece, int i_rho, int i_drho,
                          int i_theta, int i_dtheta, int i_notes)
{
OBJECT *obj1;
int *index_obj1, tabular_only;
int i, nobj1 = 0, with_wds_data;
int quadrant_correction;

if((obj1 = (OBJECT *)malloc(NOBJ_MAX * sizeof(OBJECT))) == NULL) {
  printf("astrom_compare_files/Fatal error allocating memory space for OBJECT: nobj_max=%d\n",
          NOBJ_MAX);
  exit(-1);
  }
// Initialize the number of measurements to zero
for(i = 0; i < NOBJ_MAX; i++) (obj1[i]).nmeas = 0;

if((index_obj1 = (int *)malloc((NOBJ_MAX) * sizeof(int))) == NULL) {
  printf("astrom_compare_files/Fatal error allocating memory space for index_obj: nobj_max=%d\n",
          NOBJ_MAX);
  exit(-1);
  }
with_wds_data = 0;
astrom_read_measures_cmp(fp_in1, filein2, obj1, &nobj1, 
                         i_filename, i_date, i_filter, i_eyepiece, i_rho, 
                         i_drho, i_theta, i_dtheta, i_notes, with_wds_data);

#ifdef DEBUG1
if(nobj1 > 1) {
printf("Returned by astrom_read_measures: nobj1 = %d nmeas=%d rho=%.2f theta=%.2f\n", 
        nobj1, obj1[nobj1-1].nmeas, obj1[nobj1-1].meas[0].rho, 
        obj1[nobj1-1].meas[0].theta);
}
#endif

/* Sort the objects according to their name: */
astrom_name_sort_objects(obj1, index_obj1, nobj1);

/* For big tables, should set this parameter to 1: */
tabular_only = 0;
// Correction for the quadrant (if there is an indication of the quadrant)
quadrant_correction = 1;
astrom_write_publi_cmp_table_gili(fp_out, obj1, index_obj1, nobj1,
                                  tabular_only, quadrant_correction);

// Number of objects, quadrants, etc...
astrom_compute_statistics(fp_out, obj1, nobj1, filein1);

// Write the LaTeX table for which the compared measures differ 
astrom_write_cmpfile(filein1, filein2, cmp_fileout, obj1, nobj1);

// Write the LaTeX table with the unresolved stars in filein2 compared to filein1
astrom_write_unrescmpfile(filein1, filein2, unres_cmp_fileout, obj1, nobj1);

free(index_obj1);
free(obj1);
return(0);
}
/*****************************************************************************
* Read the measurements and object parameters from the input file
*
* INPUT:
* i_filename: column nber of the filename used for this measurement
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
* i_notes: column nber with the notes
* nobj: number of objects already entered into *obj
* with_wds_data: flag set to one if wds data are in the input file
*
* OUTPUT:
* obj: OBJECT structure
* nobj: total number of objects entered into *obj
*****************************************************************************/
static int astrom_read_measures_cmp(FILE *fp_in, char *filein2,
                         OBJECT *obj, int *nobj, int i_filename, int i_date,
                         int i_filter, int i_eyepiece, int i_rho,
                         int i_drho, int i_theta, int i_dtheta, int i_notes,
                         int with_wds_data)
{
char b_in[NMAX], b_data[NMAX];
char wds_name[40], discov_name[40];
int inside_array, line_is_opened, status, orbit, line_to_reject, i_obj;
int input_with_header, comments_wanted;
char *pc, *pc1;
char obj_name1[64];
double epoch1, rho2, err_rho2, theta2, err_theta2;

#ifdef DEBUG1
printf("astrom_read_measures: input/nobj=%d\n", *nobj);
#endif

/* JLP 2014: automatic check with the first line */
/* Read first line: */
  fgets(b_in,170,fp_in);
/* input_with_header: 1 if input Latex file has a full header
*                    0 if no header at all
*/
 if(!strncmp(b_in,"% FILE_WITH_HEADER", 18)) input_with_header = 1;
 else input_with_header = 0;

/* Assume we are inside the array if (input_with_header == 0): */
if(input_with_header == 0)
  inside_array = 1;
 else
  inside_array = 0;

line_to_reject = 0;
line_is_opened = 0;
wds_name[0] = '\0';
discov_name[0] = '\0';
orbit = 0;
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

    if(!strncmp(b_in,"\\begin{tabular}",15)){
       inside_array = 1;
       strcpy(b_in,"\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|} \n");
        }
    else if(!strncmp(b_in,"\\end{tabular}",13)){
       inside_array = 0;
       }
    else if(inside_array && (b_in[0] != '%' && b_in[1] != '%')
            && strncmp(b_in,"\\hline",6)) {
       if(!line_is_opened) {
         strcpy(b_data, b_in);
/* Fill the data array with the next line */
       } else {
/* Look for the first zero (end of string marker) in data buffer */
         b_data[119] = '\0';
         pc1 = b_data;
         while(*pc1) pc1++;
         pc1--;
/* Then copy the second line from there*/
         strcpy(pc1, b_in);
       }
/* Check if this line is ended with "\\": */
       line_is_opened = 1;
       pc = b_data;
       while(*pc) {
         if(!strncmp(pc,"\\\\",2)){
           line_is_opened = 0;
           pc += 2; *pc = '\n'; pc++; *pc = '\0';
           break;
           }
         pc++;
         }
     if(!line_is_opened) {
#ifdef DEBUG1
printf("\n astrom_read_measures/New line: >%s<\n", b_data);
printf(" astrom_read_measures/Trying to add a new object, current nobj=%d (wds_data=%d)\n",
       *nobj, with_wds_data);
#endif
/* Try to add a new object and add the parameters of this object to obj: */
       if(with_wds_data == 1)
            status = astrom_add_new_object_with_wds_data(b_data, obj, nobj,
                                                         i_notes);
       else
            status = astrom_add_new_object_without_wds_data(b_data, obj, nobj);
       if(status == 0) {
/* Index of current object in OBJECT structure array: */
          i_obj = *nobj - 1;
          (obj[i_obj]).nmeas = 0;
          }

/* Try to add a new measure: */
// JLP2020 it was != 0 : why ??
       if((status == 0) && (*nobj > 1)) {
         status = astrom_check_measure(b_data, i_eyepiece, i_rho, i_drho,
                                       i_theta, i_dtheta);
#ifdef DEBUG1
printf("astrom_read_measures/Adding new measurement for object #i_obj=%d nm=%d (nobj=%d, astrom_check_measure: status=%d)\n",
         i_obj, (obj[i_obj]).nmeas, *nobj, status);
#endif
          if(status == 0) {
           comments_wanted = 0;
           astrom_add_new_measure(b_data, obj, i_obj, i_filename, i_date,
                           i_filter, i_eyepiece, i_rho, i_drho, i_theta,
                           i_dtheta, i_notes, comments_wanted);
// theta is always in [0, 360.]:
          if(obj[i_obj].meas[0].theta < 0.) 
                obj[i_obj].meas[0].theta += 360.;
          else if(obj[i_obj].meas[0].theta > 360.) 
                obj[i_obj].meas[0].theta -= 360.;
#ifdef DEBUG
           printf("astrom_read_meas/discov_name=%s nmeas=%d rho=%.2f+/-%.2f theta=%.2f+/-%.2f\n", 
                  obj[i_obj].discov_name, 
                  obj[i_obj].nmeas, 
                  obj[i_obj].meas[0].rho, 
                  obj[i_obj].meas[0].drho, 
                  obj[i_obj].meas[0].theta,
                  obj[i_obj].meas[0].dtheta);
#endif
           strcpy(obj_name1, (obj[i_obj]).discov_name);
           epoch1 = (obj[i_obj]).meas[0].epoch;
// Return 0 when obj_name1 was found in filein2
           status = astrom_get_meas_from_file(filein2, obj_name1,
                       epoch1, i_filename, i_date, i_filter,
                       i_eyepiece, i_rho, i_drho, i_theta, i_dtheta, i_notes,
                       &rho2, &err_rho2, &theta2, &err_theta2);
           if(status == 0) {
#ifdef DEBUG
            printf("astrom_get_meas/%s in %s: rho=%.2f+/-%.2f theta=%.2f+/-%.2f\n",
                  obj_name1, filein2, rho2, err_rho2, theta2, err_theta2);
#endif
             obj[i_obj].cmp_meas[0].rho = rho2; 
             obj[i_obj].cmp_meas[0].drho = err_rho2;
// Change the theta value if needed:
             if((obj[i_obj].meas[0].theta != NO_DATA) && (theta2 != NO_DATA))
                astrom_change_new_theta(&theta2, obj[i_obj].meas[0].theta); 
             obj[i_obj].cmp_meas[0].theta = theta2;
             obj[i_obj].cmp_meas[0].dtheta = err_theta2;
            } else {
// Set the measurements to -NO_DATA if the file was not found
             obj[i_obj].cmp_meas[0].rho = -NO_DATA; 
             obj[i_obj].cmp_meas[0].drho = -NO_DATA;
             obj[i_obj].cmp_meas[0].theta = -NO_DATA;
             obj[i_obj].cmp_meas[0].dtheta = -NO_DATA;
#ifdef DEBUG
            printf("astrom_read_meas/Warning: object %s not found in %s\n", 
                   obj_name1, filein2);
#endif
            }
            } // case status from check_measure = 0
         }

       } /* EOF !line_is_opened */
    } // EOF line tabular
  } /* EOF if fgets() */
} /* EOF while loop */
return(0);
}
/*********************************************************************
* Return 0 when obj_name1 was found in filein2
*
*********************************************************************/
static int astrom_get_meas_from_file(char *filein2, char *obj_name1,
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
printf("\n astrom_get_meas_from_file/looking for fname=>%s< in %s\n",
         obj_name1, filein2);
#endif

if((fp_in2 = fopen(filein2,"r")) == NULL) {
  fprintf(stderr, " Fatal error opening input file2 %s \n", filein2);
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
    status = astrom_check_measure(b_data, i_eyepiece, i_rho, i_drho, i_theta, 
                                i_dtheta);
    if(status == 0) {
#ifdef DEBUG1
printf("astrom_get_meas_from_file/astrom_check_measure: b_data=%s status=%d)\n",
         b_data, status);
#endif
          astrom_read_new_measure(b_data, i_filename, i_date, i_filter, 
                                  i_eyepiece, i_rho, i_drho, i_theta,
                                  i_dtheta, i_notes, filename2, notes2, 
                                  &epoch2, filter2, &eyepiece2,
                                  rho2, err_rho2, theta2, err_theta2);
// Remove the extension ".fits" for instance:
          jlp_remove_ext_string(filename2, 64);
          jlp_compact_string(filename2, 64);
          if(!strcmp(obj_name1, filename2)) {
#ifdef DEBUG1
            printf("astrom_read_measures/in %s: fname=%s rho=%.2f+/-%.2f theta=%.2f+/-%.2f\n",
                  filein2, filename2, *rho2, *err_rho2, *theta2, *err_theta2);
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
static int astrom_write_publi_cmp_table_gili(FILE *fp_out, 
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
    astrom_write_header1(fp_out, nlines, tabular_only);
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
             wds_name, obj[io].discov_name, me->epoch, me->eyepiece); //  me->notes);
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

// Check if theta value (me->theta) is compatible with the quadrant value (me->quadrant):
  good_q = astrom_quadrant_is_consistent(me);
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
           status = astrom_correct_theta_with_WDS_CHARA(&(obj[io]),me);
           if((cmp_me->theta != NO_DATA) && (cmp_me->theta != -NO_DATA))
              astrom_correct_theta_with_WDS_CHARA(&(obj[io]), cmp_me);
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
     astrom_preformat_wds_name(obj[io].wds, wds_name);
/*** Format for first line of an object */
/* October 2008: 3 decimals for the epoch */
  if(me->rho != NO_DATA)
    fprintf(fp_out,"%s & %s & %.3f & %d & %.3f & %.3f & %.1f%s & %.1f & %s &\\\\\n",
         wds_name, obj[io].discov_name, me->epoch, me->eyepiece, me->rho, me->drho,
         me->theta, qflag, me->dtheta, dmag_string); //, me->notes, q_notes);
  else
    fprintf(fp_out,"%s & %s & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata & \\nodata &  \\\\\n",
         wds_name, obj[io].discov_name, me->epoch, me->eyepiece); //, me->notes);
// File was not found when cmp_me->rho == -NO_DATA
  if(cmp_me->rho != -NO_DATA) {
  if(cmp_me->rho != NO_DATA)
    fprintf(fp_out,"cmp:%s & %s & %.3f & %d & %.3f & %.3f & %.1f%s & %.1f & %s &\\\\\n",
         wds_name, obj[io].discov_name, me->epoch, me->eyepiece, cmp_me->rho, cmp_me->drho,
         cmp_me->theta, qflag, cmp_me->dtheta, dmag_string); //, me->notes, q_notes);
  else
    fprintf(fp_out,"cmp:%s & %s & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata & \\nodata &  \\\\\n",
         wds_name, obj[io].discov_name, me->epoch, me->eyepiece); //, me->notes);
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
         wds_name, obj[io].discov_name, me->epoch, me->eyepiece, me->rho, me->drho,
         me->theta, qflag, me->dtheta, dmag_string); //, me->notes, q_notes);
  else
    fprintf(fp_out,"%s & %s & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata & \\nodata &  \\\\\n",
         wds_name, obj[io].discov_name, me->epoch, me->eyepiece); //, me->notes);
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
static int astrom_change_new_theta(double *theta2, double theta1) 
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
static int astrom_write_cmpfile(char *filein1, char *filein2, 
                                 char *cmp_fileout, OBJECT *obj1, int nobj1)
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
  printf("astrom_compare_files/Fatal error allocating memory space for OBJECT: nobj_max=%d\n",
          NOBJ_MAX);
  exit(-1);
  }
// Initialize the number of measurements to zero
for(i = 0; i < NOBJ_MAX; i++) (obj2[i]).nmeas = 0;

if((index_obj2 = (int *)malloc((NOBJ_MAX) * sizeof(int))) == NULL) {
  printf("astrom_compare_files/Fatal error allocating memory space for index_obj: nobj_max=%d\n",
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
         printf("discov_name=%s rho=%.2f theta=%.2f Drho=%.2f Dtheta=%.2f\n",
               (obj1[i]).discov_name, me1->rho, me1->theta,
               ABS(me1->rho - cmp_me1->rho),
               ABS(me1->theta - cmp_me1->theta));
             strcpy((obj2[k]).wds, (obj1[i]).wds);
             strcpy((obj2[k]).discov_name, (obj1[i]).discov_name);
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
printf("Returned by astrom_read_measures: nobj2=%d nmeas=%d rho=%.2f theta=%.2f\n",
        nobj2, obj2[nobj2-1].nmeas, obj2[nobj2-1].meas[0].rho,
        obj2[nobj2-1].meas[0].theta);
}
#ifdef DEBUG1
#endif


/* Sort the objects according to their name: */
astrom_name_sort_objects(obj2, index_obj2, nobj2);
if((fp_cmp_out = fopen(cmp_fileout,"w")) == NULL) {
  fprintf(stderr, "astrom_create_cmpfile/Error opening output cmp file %s \n",
        cmp_fileout);
  return(-1);
  }

fprintf(fp_cmp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(fp_cmp_out,"%%%% astrom_compare file1=%s and file2=%s\n", filein1, filein2);
fprintf(fp_cmp_out,"%%%% Measurements that differ nobj2=%d (nobj1=%d) with Drho=%.2f and Dtheta=%.2f\n", 
        nobj2, nobj1, drho_max, dtheta_max);
fprintf(fp_cmp_out,"%%%% JLP / Version of 31/03/2020 \n");
fprintf(fp_cmp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

/* For big tables, should set this parameter to 1: */
tabular_only = 0;
// NO correction for the quadrant (already done)
quadrant_correction = 0;
astrom_write_publi_cmp_table_gili(fp_cmp_out, obj2, index_obj2, nobj2,
                         tabular_only, quadrant_correction);

fclose(fp_cmp_out);
free(index_obj2);
free(obj2);
return(0);
}
/****************************************************************
* Output of the Latex table with the unresolved stars of filein2
* that have been resolved in filein1
*
****************************************************************/
static int astrom_write_unrescmpfile(char *filein1, char *filein2, 
                                     char *unres_cmp_fileout, OBJECT *obj1, 
                                     int nobj1)
{
FILE *fp_cmp_out;
OBJECT *obj2;
MEASURE *me1, *cmp_me1, *me2, *cmp_me2;
int *index_obj2, tabular_only, nm1;
int i, j, k, nobj2 = 0;
int quadrant_correction;

if((obj2 = (OBJECT *)malloc(NOBJ_MAX * sizeof(OBJECT))) == NULL) {
  printf("astrom_compare_files/Fatal error allocating memory space for OBJECT: nobj_max=%d\n",
          NOBJ_MAX);
  exit(-1);
  }
// Initialize the number of measurements to zero
for(i = 0; i < NOBJ_MAX; i++) (obj2[i]).nmeas = 0;

if((index_obj2 = (int *)malloc((NOBJ_MAX) * sizeof(int))) == NULL) {
  printf("astrom_compare_files/Fatal error allocating memory space for index_obj: nobj_max=%d\n",
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
// rho=-NO_DATA when file was not found...
      && (cmp_me1->rho != -NO_DATA) && (cmp_me1->theta != -NO_DATA)) {
// rho=NO_DATA when not resolved...
         if((cmp_me1->rho != NO_DATA) && (cmp_me1->theta != NO_DATA)) {
             printf("discov_name=%s rho=%.2f theta=%.2f\n",
               (obj1[i]).discov_name, me1->rho, me1->theta);
             strcpy((obj2[k]).wds, (obj1[i]).wds);
             strcpy((obj2[k]).discov_name, (obj1[i]).discov_name);
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
             cmp_me2->rho = NO_DATA; 
             cmp_me2->drho = NO_DATA; 
             cmp_me2->theta = NO_DATA; 
             cmp_me2->dtheta = NO_DATA; 
             k++;
             } // > drho_max
      } // not flagged out 
    } // loop on j nm1
  } // loop on i nobj
nobj2 = k;
if(nobj2 > 1) {
printf("Returned by astrom_read_measures: nobj2=%d nmeas=%d rho=%.2f theta=%.2f\n",
        nobj2, obj2[nobj2-1].nmeas, obj2[nobj2-1].meas[0].rho,
        obj2[nobj2-1].meas[0].theta);
}
#ifdef DEBUG1
#endif


/* Sort the objects according to their name: */
astrom_name_sort_objects(obj2, index_obj2, nobj2);
if((fp_cmp_out = fopen(unres_cmp_fileout,"w")) == NULL) {
  fprintf(stderr, "astrom_write_unrescmpfile/Error opening output unrescmp file %s \n",
        unres_cmp_fileout);
  return(-1);
  }

fprintf(fp_cmp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(fp_cmp_out,"%%%% astrom_compare file1=%s and file2=%s\n", filein1, filein2);
fprintf(fp_cmp_out,"%%%% Unresolved in file2 nobj2=%d (nobj1=%d)\n", 
        nobj2, nobj1);
fprintf(fp_cmp_out,"%%%% JLP / Version of 31/03/2020 \n");
fprintf(fp_cmp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

/* For big tables, should set this parameter to 1: */
tabular_only = 0;
// NO correction for the quadrant (already done)
quadrant_correction = 0;
astrom_write_publi_cmp_table_gili(fp_cmp_out, obj2, index_obj2, nobj2,
                         tabular_only, quadrant_correction);

fclose(fp_cmp_out);
free(index_obj2);
free(obj2);
return(0);
}
/**************************************************************************
* Write the LateX header of the table
***************************************************************************/
static int astrom_write_header1(FILE *fp_out, int nlines, int tabular_only)
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

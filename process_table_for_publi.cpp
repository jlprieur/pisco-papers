/************************************************************************
* "process_table1_for_publi.cpp"
*
* Process "full_tab_calib.tex"
* - iop=0 look for the smallest and largest separations
* - iop=1 set minimum rho and theta errors for small separations (smaller than rho_diff)
* - iop=2 look for separations smaller than min_sep and extract this table
* - iop=3 look for separations smaller than min_sep and extract this table (with residuals of orbits)
* - iop=4 look for objects with same WDS name that have been observed twice or more, that are discrepant in rho and theta and extract this table
* - iop=5 look for objects with same WDS name that have been observed twice or more, with small range in rho and theta but with different discov names and extract this table
* - iop=6 look for objects with odd names and extract this table 
* - iop=7 look for new doubles (ND or nd) and extract this table 
* - iop=8 look for large residuals and extract this table 
* - iop=9 look for unsorted observations of the same objects according to the epoch 
*
* JLP 
* Version 30/10/2020
*************************************************************************/
#include <stdio.h> 
#include <stdlib.h> // exit(-1) 
#include <string.h> 
#include <ctype.h>  // isalpha(), isdigit() 
#include <time.h>   // time_t 
#include "latex_utils.h" // latex_get_column_item(), latex_remove_column()
#include "jlp_string.h"  // (in jlplib/jlp_fits/ ) jlp_cleanup_..
#include "jlp_catalog_utils.h"  // ABS 

/*
#define DEBUG
*/

typedef struct {
char wds[40]; 
char discov[40]; 
double epoch;
double eyep;
double rho_meas;
double drho_meas;
double theta_meas;
double dtheta_meas;
double dm_meas;
char notes[40];
char orbit_ref[40];
int orbit_grade;
double rho_res;
double theta_res;
} PSC_MEAS;

static int count_objects(char *in_fname, double rho_diff);
static int modif_errors_table1(char *in_fname, char *out_fname, double rho_diff,
                        double drhodiff_mini, double dthetadiff_mini);
static int extract_closest_table1(char *in_fname, char *out_fname, 
                                   double rho_c, int resid_only);
static int extract_twomeas_table1(char *in_fname, char *out_fname, int iopt);
static int extract_new_doubles_table1(char *in_fname, char *out_fname);
static int check_discrepant_rho_meas(char *obs_discov, double *obs_rho,
                                     int nobs, int same_discov, 
                                     int *out_result);
static int check_discrepant_theta_meas(char *obs_discov, double *obs_theta,
                                     int nobs, int same_discov, 
                                     int *out_result);
static int check_sorted_epoch_meas(char *obs_discov, double *obs_epoch, 
                                    int nobs, int same_discov, int *out_result);
static int check_same_discov_names(char *obs_discov, int nobs, int *out_result);
static int check_if_oddname(char *discov_name0, int nlength0, int *name_is_odd);
static int check_if_new_double(char *discov_name0, int nlength0, 
                               int *is_new_double);
static int extract_oddnames_table1(char *in_fname, char *out_fname);
static int extract_large_resid_table1(char *in_fname, char *out_fname, 
                                      double rho_res_min, double theta_res_min);
static int init_psc_meas(PSC_MEAS *psc0, char* wds0, char* dsc0, double epoch0,
                  double eyep0, double rmeas0, 
                  double drmeas0, double tmeas0, double dtmeas0, 
                  double dm0, double rres0, double tres0);
static int printf_psc_meas(PSC_MEAS psc0, char *label);
static int sort_measures_table1(char *in_fname, char *out_fname);

static int jlp_clean_dollars(char *str0, int str_len0);
/*********************************************************************
* Remove "$" in latex calib file. (used for negative values in LaTeX)
**********************************************************************/
static int jlp_clean_dollars(char *str0, int str_len0)
 {char buffer[256];
  int i, k;
  if(str_len0 > 256){ fprintf(stderr, "jlp_clean_dollars: Fatal error; str_len0=%d\n",
                               str_len0);
                      exit(-1);
                    }
  k = 0;
  for(i = 0; i < str_len0; i++) {
    if(str0[i] != '$') {
      buffer[k] = str0[i];
      k++;
      }
    }
  buffer[k] = '\0';
  strncpy(str0, buffer, str_len0);
 return(0);
 }
/***********************************************************************
* init_psc_meas
*
* typedef struct {
char wds[40]; 
char discov[40]; 
double epoch;
double eyep;
double rho_meas;
double drho_meas;
double theta_meas;
double dtheta_meas;
double dm_meas;
char notes[40];
char orbit_ref[40];
int orbit_grade;
double rho_res;
double theta_res;
* } PSC_MEAS;
*
************************************************************************/
static int init_psc_meas(PSC_MEAS *psc0, char* wds0, char* dsc0, double epoch0,
                  double eyep0, double rmeas0, 
                  double drmeas0, double tmeas0, double dtmeas0, 
                  double dm0, double rres0, double tres0)
{ 
// WDS name: 
strcpy(psc0->wds, wds0);

// Discover's name: 
strcpy(psc0->discov, dsc0);

// Epoch:
psc0->epoch = epoch0;

// Eyepiece:
psc0->eyep = eyep0;

// Rho measurements:
psc0->rho_meas = rmeas0;
psc0->drho_meas = drmeas0;

// Theta measurements:
psc0->theta_meas = tmeas0;
psc0->dtheta_meas = dtmeas0;

// Delta_mag:
psc0->dm_meas = dm0;

// Notes 
strcpy(psc0->notes, "");

// Orbit bibliographic reference: 
strcpy(psc0->orbit_ref, "");

// Orbit_grade:
psc0->orbit_grade = 0;

// Residuals:
psc0->rho_res = rres0;
psc0->theta_res = tres0;

return(0);
}
/***********************************************************************
* printf_psc_meas
*
************************************************************************/
static int printf_psc_meas(PSC_MEAS psc0, char *label)
{ 
printf("%s wds_name=%s discov_name=%s rho=%f drho=%f theta=%f dtheta=%f\n", 
       label,  psc0.wds, psc0.discov, psc0.rho_meas, psc0.drho_meas,
       psc0.theta_meas, psc0.dtheta_meas);
printf("%s wds_name=%s notes=%s orbit_ref=%f orbit_grade rho_res=%f theta_res=%f\n", 
       label,  psc0.wds, psc0.notes, psc0.orbit_ref, psc0.orbit_grade,
       psc0.rho_res, psc0.theta_res);
return(0);
}

/***********************************************************************
*
************************************************************************/
int main(int argc, char *argv[])
{
char in_fname[128], out_fname[128];
double rho_diff = 0.16, drhodiff_mini = 0., dthetadiff_mini = 0., rho_c;
double rho_res_min = 0., theta_res_min = 0.;
int iopt, resid_only;

if(argc != 4) {
  printf("Syntax:\n"); 
  printf("Option0: look for min and max separation\n");
  printf("Option1: set drho mini and dtheta mini as drhodiff_mini and dthetadiff_mini when rho less than rho_diff\n");
  printf("process_table1_for_publi old_table new_table 1,rho_diff,drhodiff_mini,dthetadiff_mini\n");
  printf("Option2: extract sub-table with separation smaller than rhodiff\n");
  printf("process_table1_for_publi old_table new_table 2,rho_diff \n");
  printf("Option3: extract sub-table with separation smaller than rhodiff (with orbits)\n");
  printf("process_table1_for_publi old_table new_table 3,rho_diff \n");
  printf("Option4: extract sub-table with WDS objects observed twice or more times with discrepant rho or theta\n");
  printf("process_table1_for_publi old_table new_table 4 \n");
  printf("Option5: extract sub-table with WDS objects observed twice or more times with discrepant discov names (and small range in rho and theta)\n");
  printf("process_table1_for_publi old_table new_table 5 \n");
  printf("Option6: extract sub-table with objects with odd names\n");
  printf("process_table1_for_publi old_table new_table 6 \n");
  printf("Option7: extract sub-table with nd or ND (new doubles)\n");
  printf("process_table1_for_publi old_table new_table 7 \n");
  printf("Option8: extract sub-table of largest residuals \n");
  printf("Option9: extract sub-table to sort the measures of the same object according to the epoch\n");
  printf("process_table1_for_publi old_table new_table 8,rho_res_min,theta_res_min \n");
  return(-1);
}
strcpy(in_fname, argv[1]);
strcpy(out_fname, argv[2]);
sscanf(argv[3], "%d", &iopt);
if(iopt == 1)
   sscanf(argv[3], "%d,%lf,%lf,%lf", &iopt, &rho_diff, &drhodiff_mini, &dthetadiff_mini);
else if((iopt == 2) || (iopt == 3))
   sscanf(argv[3], "%d,%lf", &iopt, &rho_diff);
else if(iopt == 8)
   sscanf(argv[3], "%d,%lf,%lf", &iopt, &rho_res_min, &theta_res_min);

/*
printf("%s\n", argv[3]);
printf("OK: in_table=%s output=%s\n", 
        in_fname, out_fname); 
*/
/************
if(iopt == 1) {
  printf("OK: iopt=%d rho_diff=%f drhodiff_mini=%f dthetadiff_mini=%f\n", 
         iopt, rho_diff, drhodiff_mini, dthetadiff_mini); 
} else if((iopt == 2)||(iopt == 3)) {
  printf("OK: iopt=%d rho_diff=%f\n", 
         iopt, rho_diff); 
} else if(iopt == 8) {
  printf("OK: iopt=%d rho_res_min=%f theta_res_min=%f\n", 
         iopt, rho_res_min, theta_res_min); 
} else {
  printf("OK: iopt=%d\n", 
         iopt); 
}
*********/

switch(iopt) {
  case 0:
  default:
     count_objects(in_fname, rho_diff);
     break;
  case 1:
     modif_errors_table1(in_fname, out_fname, rho_diff, drhodiff_mini, dthetadiff_mini); 
     break;
// Closest observations, with and without residuals
  case 2:
     rho_c = rho_diff;
     resid_only = 0;
     extract_closest_table1(in_fname, out_fname, rho_c, resid_only);
     break;
// Closest observations with residuals
  case 3:
     rho_c = rho_diff;
     resid_only = 1;
     extract_closest_table1(in_fname, out_fname, rho_c, resid_only);
     break;
  case 4:
  case 5:
     extract_twomeas_table1(in_fname, out_fname, iopt);
     break;
  case 6:
     extract_oddnames_table1(in_fname, out_fname);
     break;
  case 7:
     extract_new_doubles_table1(in_fname, out_fname);
     break;
  case 8:
     extract_large_resid_table1(in_fname, out_fname, 
                                rho_res_min, theta_res_min);
     break;
  case 9:
     sort_measures_table1(in_fname, out_fname);
     break;
  }
return(0);
}
/************************************************************************
* Scan the input table and make the modifications 
* - look for the smallest and largest separations
* - look for the smallest and largest sep. errors (when rho < rho_diff) 
* - count de number of objects, the number of observations
*
* INPUT:
* in_fname: input table filename
* rho_diff: telescope diffraction limit for rho, in arcseconds
*
*************************************************************************/
static int count_objects(char *in_fname, double rho_diff)
{
char in_line[256], in_line2[256], in_line3[256]; 
char discov_name[40], wds_name[40], label[64], old_discov_name[40]; 
char buffer[256], rho_meas[40], theta_meas[40], orbit_ref[40], *pc;
int iline, status, verbose_if_error = 0, nobj, nmeas, nunres, nresid, ndmag;
int nlines_meas, nquad, orbit_grade;
double rho_val, drho_val, theta_val, dtheta_val, dmag_val, dval;
char drho_meas[40], dtheta_meas[40], dmag_str[40];
double rho_min, rho_max, drhodiff_min, drhodiff_max;
double dthetadiff_min, dthetadiff_max;
int i, irho = 5, i_drho = 6, resolved;
PSC_MEAS psc_rho_min, psc_rho_max, psc_drhodiff_min, psc_drhodiff_max;
PSC_MEAS psc_dthetadiff_min, psc_dthetadiff_max;
FILE *fp_in;
time_t ttime = time(NULL);

/* Open input table: */
if((fp_in = fopen(in_fname, "r")) == NULL) {
   fprintf(stderr, "count_objects/Error opening input table %s\n",
           in_fname);
    return(-1);
  }

rho_min = 1000.;
rho_max = -1.;
drhodiff_min = 1000.;
drhodiff_max = -1.;
dthetadiff_min = 1000.;
dthetadiff_max = -1.;
nobj = 0;
nmeas = 0;
nunres = 0;
nresid = 0;
ndmag = 0;
nquad = 0;
iline = 0;
nlines_meas = 0;
resolved = -1;
strcpy(old_discov_name, "");
while(!feof(fp_in)) {
  if(fgets(in_line, 256, fp_in)) {
    iline++;
// Remove all the non-printable characters and the end of line '\n' 
// from input line:
      jlp_cleanup_string(in_line, 256);
#ifdef DEBUG
      printf("count_objects/in_line=%s\n nobj=%d nmeas=%d nunres=%d\n", in_line, nobj, nmeas, nunres);
#endif
      strcpy(in_line3, in_line);

// Good lines start with a digit (WDS names...)
// Lines starting with % are ignored
    if(isdigit(in_line[0])) { 
     nlines_meas++;
/* Get wds_name from column 1: */
     latex_get_column_item(in_line, wds_name, 1, verbose_if_error);
/* Get discov_name from column 2: */
     latex_get_column_item(in_line, discov_name, 2, verbose_if_error);
     if(strcmp(old_discov_name, discov_name) != 0) {
        nobj++;
        strcpy(old_discov_name, discov_name);
        resolved = 0;
     }
// Look fo "nodata" in line, increase the number of unresolved objects
     if(strstr(in_line, "nodata") != NULL) nunres++;
// Look for "^*" in line, increase the number of quadrants 
     if(strstr(in_line, "^*") != NULL) nquad++;

/* Get rho measure from column 5: */
/* Get drho estimate from column 6: */
      verbose_if_error = 0;
      rho_val = -1;
      drho_val = -1;

      if(latex_get_column_item(in_line, rho_meas, 5, verbose_if_error) == 0){
        if(sscanf(rho_meas, "%lf", &dval) == 1) {
          rho_val = dval;
          nmeas++;
          if(resolved == 0) resolved = 1;
          }
        if(latex_get_column_item(in_line, drho_meas, 6, verbose_if_error) == 0){
          if(sscanf(drho_meas, "%lf", &dval) == 1) drho_val = dval;
          }
      }
/* Get theta measure from column 7: */
/* Get dtheta estimate from column 8: */
      verbose_if_error = 0;
      theta_val = -1;
      dtheta_val = -1;

      if(latex_get_column_item(in_line, theta_meas, 7, verbose_if_error) == 0){
        if(sscanf(theta_meas, "%lf", &dval) == 1) theta_val = dval;
        if(latex_get_column_item(in_line, dtheta_meas, 8, verbose_if_error) == 0){
          if(sscanf(dtheta_meas, "%lf", &dval) == 1) dtheta_val = dval;
          }
      }
/* Get dmag from column 9: */
      verbose_if_error = 0;
      dmag_val = -1;
      if(latex_get_column_item(in_line, dmag_str, 9, verbose_if_error) == 0){
        if(sscanf(dmag_str, "%lf", &dval) == 1) {
          dmag_val  = dval;
#ifdef DEBUG1
        printf("discov_name=%s dmag=%s (%f)\n", discov_name, dmag_str, dmag_val);
#endif
          ndmag++;
        }
      }
/* Get orbit ref from column 11: */
      if(latex_get_column_item(in_line, orbit_ref, 11, verbose_if_error) == 0){
#ifdef DEBUG1
        printf("discov_name=%s orbit_ref=%s\n", discov_name, orbit_ref);
#endif
        nresid++;
        }
#ifdef DEBUG
       printf("wds_name=%s discov_name=%s rho=%f drho=%f\n",
               wds_name, discov_name, rho_val, drho_val);
#endif
       if((rho_val > 0) && (drho_val > 0.)) {
#ifdef DEBUG
         printf("OK rho_diff=%f, rho_min=%f rho_max=%f drhodiff_min=%f dthetadiff_min=%f\n",
                 rho_diff, rho_min, rho_max, drhodiff_min, dthetadiff_min);
#endif
       if(rho_val < rho_min) {
/*
int init_psc_meas(PSC_MEAS *psc0, char* wds0, char* dsc0, double epoch0,
                  double eyep0, double rmeas0,
                  double drmeas0, double tmeas0, double dtmeas0,
                  double dm0, double rres0, double tres0);
*/
         init_psc_meas(&psc_rho_min, wds_name, discov_name, 0., 0.,
                       rho_val, drho_val, theta_val, dtheta_val, 0., 0., 0.);
         rho_min = rho_val;
         }
       if(rho_val > rho_max) {
         init_psc_meas(&psc_rho_max, wds_name, discov_name, 0., 0.,
                       rho_val, drho_val, theta_val, dtheta_val, 0., 0., 0.);
         rho_max = rho_val;
         }
      if((rho_val < rho_diff) && (drho_val < drhodiff_min)) {
         init_psc_meas(&psc_drhodiff_min, wds_name, discov_name, 0., 0.,
                       rho_val, drho_val, theta_val, dtheta_val, 0., 0., 0.);
         drhodiff_min = drho_val;
         }
      if((rho_val < rho_diff) && (drho_val > drhodiff_max)) {
         init_psc_meas(&psc_drhodiff_max, wds_name, discov_name, 0., 0.,
                       rho_val, drho_val, theta_val, dtheta_val, 0., 0., 0.);
         drhodiff_max = drho_val;
         }
      if((rho_val < rho_diff) && (dtheta_val < dthetadiff_min)) {
         init_psc_meas(&psc_dthetadiff_min, wds_name, discov_name, 0., 0.,
                       rho_val, drho_val, theta_val, dtheta_val, 0., 0., 0.);
         dthetadiff_min = dtheta_val;
         }
      if((rho_val < rho_diff) && (dtheta_val > dthetadiff_max)) {
         init_psc_meas(&psc_dthetadiff_max, wds_name, discov_name, 0., 0.,
                       rho_val, drho_val, theta_val, dtheta_val, 0., 0., 0.);
         dthetadiff_max = dtheta_val;
         }
      }

#ifdef DEBUG
      printf("OK rho_diff=%f, rho_min=%f rho_max=%f \n", 
              rho_diff, rho_min, rho_max);
#endif

    }/* EOF if !isdigit % ... */

  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("count_objects: %d lines sucessfully read and processed, nlines_meas=%d\n", 
        iline, nlines_meas);
printf("count_objects: n_objects=%d n_meas(resolved)=%d n_unres=%d \n", 
        nobj, nmeas, nunres);
printf("count_objects: n_resid=%d n_mag=%d n_quad=%d\n", 
        nresid, ndmag, nquad);

if(rho_min < 100.) {
   strcpy(label, "smallest rho: ");
   printf_psc_meas(psc_rho_min, label);
}
if(rho_max > -1.) {
   strcpy(label, "largest rho: ");
   printf_psc_meas(psc_rho_max, label);
}
if(drhodiff_min < 100.) {
   strcpy(label, "smallest drho (when rho < rho_diff): ");
   printf_psc_meas(psc_drhodiff_min, label);
}
if(drhodiff_max > -1.) {
   strcpy(label, "largest drho (when rho < rho_diff): ");
   printf_psc_meas(psc_drhodiff_max, label);
}
if(dthetadiff_min < 100.) {
   strcpy(label, "smallest dtheta (when rho < rho_diff): ");
   printf_psc_meas(psc_dthetadiff_min, label);
}
if(dthetadiff_max > -1.) {
   strcpy(label, "largest dtheta (when rho < rho_diff): ");
   printf_psc_meas(psc_dthetadiff_max, label);
}

/* Close opened files:
*/
fclose(fp_in);

return(0);
}
/************************************************************************
* Scan the input table and make the modifications 
* - set minimum rho error (drhodiff_mini) for separations smaller than rho_diff
* - set minimum theta error (dthetadiff_mini) for separations smaller than rho_diff
*
* INPUT:
* in_fname: input table filename
* out_fname: output Latex filenme
* rho_diff: telescope diffraction limit for rho, in arcseconds
*
*************************************************************************/
static int modif_errors_table1(char *in_fname, char *out_fname, double rho_diff,
                        double drhodiff_mini, double dthetadiff_mini) 
{
char in_line[256], in_line2[256], in_line3[256]; 
char discov_name[40], wds_name[40], label[64]; 
char buffer[256], rho_meas[40], drho_meas[40], *pc, drho_item[32];
char dtheta_str[40], dtheta_item[32];
int iline, status, verbose_if_error = 0;
double rho_val, drho_val, dtheta_val, dval;
int i, irho = 5, i_drho = 6, i_dtheta = 8;
FILE *fp_in, *fp_out;
time_t ttime = time(NULL);

/* Open input table: */
if((fp_in = fopen(in_fname, "r")) == NULL) {
   fprintf(stderr, "modif_errors_table1/Fatal error opening input table %s\n",
           in_fname);
    return(-1);
  }

/* Open output table: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "modif_errors_table1/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the output Latex table: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "%% Modified table from: %s \n%% Created on %s", 
        in_fname, ctime(&ttime));
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

iline = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 256, fp_in)) {
    iline++;
// Remove all the non-printable characters and the end of line '\n' 
// from input line:
      jlp_cleanup_string(in_line, 256);
#ifdef DEBUG
      printf("in_line=%s\n", in_line);
#endif
      strcpy(in_line3, in_line);

// Good lines start with a digit (WDS names...)
// Lines starting with % are ignored
    if(isdigit(in_line[0])) { 
/* Get wds_name from column 1: */
     latex_get_column_item(in_line, wds_name, 1, verbose_if_error);
/* Get discov_name from column 2: */
     latex_get_column_item(in_line, discov_name, 2, verbose_if_error);

/* Get rho measure from column 5: */
      verbose_if_error = 0;
      rho_val = -1;
      drho_val = -1;
  if(latex_get_column_item(in_line, rho_meas, 5, verbose_if_error) == 0){
    if(sscanf(rho_meas, "%lf", &dval) == 1) rho_val = dval;
    if(latex_get_column_item(in_line, drho_meas, 6, verbose_if_error) == 0){
       if(sscanf(drho_meas, "%lf", &dval) == 1) drho_val = dval;
       }
    if((rho_val > 0) && (drho_val > 0.)) {
// Set minimum value for drho for separations smaller than rho_diff:
      if((rho_val < rho_diff) && (drho_val < drhodiff_mini)) {
/* Write new drho with 3 decimals */
        sprintf(drho_item, "%.3f", drhodiff_mini);
        latex_set_column_item(in_line3, 256, drho_item, 32, i_drho, 3);
        }
      } // rho_val > 0
/* Get dtheta measure from column i_dtheta=8: */
      verbose_if_error = 0;
      dtheta_val = -1;
     if(latex_get_column_item(in_line, dtheta_str, i_dtheta, verbose_if_error) == 0){
       if(sscanf(dtheta_str, "%lf", &dval) == 1) dtheta_val = dval;
     }
    if(dtheta_val > 0.) {
// Set minimum value for dtheta for separations smaller than rho_diff:
      if((rho_val < rho_diff) && (dtheta_val < dthetadiff_mini)) {
#ifdef DEBUG
    printf("wds_name=%s discov_name=%s rho=%f drho=%f dtheta=%f \n", 
            wds_name, discov_name, rho_val, drho_val, dtheta_val);
#endif
/* Write new dtheta with 1 decimal */
        sprintf(dtheta_item, "%.1f", dthetadiff_mini);
        latex_set_column_item(in_line3, 256, dtheta_item, 32, i_dtheta, 3);
        }
      } // dtheta_val > 0
    } // rho_meas

    }/* EOF if !isdigit % ... */

// Save to output file:
    fprintf(fp_out, "%s\n", in_line3);

  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("process_table1: %d lines sucessfully read and processed\n", 
        iline);

/* Close opened files:
*/
fclose(fp_in);
fclose(fp_out);

return(0);
}
/************************************************************************
* Scan the input table and extract the table with the smallest separations 
* Example of header:
WDS & Name & Epoch & $\rho$ & $\sigma_\rho$ & Orbit & {\scriptsize $\Delta \rho$(O-C)} & {\
scriptsize $\Delta \theta$(O-C)} \\
*
* INPUT:
* in_fname: filename of the input table 
* rho_c : separation to be used as the threhold
* resid_only : flag set to one if output objects with residuals only
*
*************************************************************************/
static int extract_closest_table1(char *in_fname, char *out_fname, 
                                  double rho_c, int resid_only)
{
char in_line[256], buffer[256]; 
char rho_meas[40], drho_meas[40], theta_meas[40], orbit_ref[40];
char discov_name[40], wds_name[40], label[64], no_data[40]; 
int iline, status, verbose_if_error = 0, nobjects, orbit_grade;
double rho_val, drho_val, theta_val, dval;
int i, irho = 5, to_output = 0;

time_t ttime = time(NULL);
FILE *fp_out, *fp_in;

strcpy(no_data, "\\nodata");
jlp_compact_string(no_data, 40);

/* Open input table: */
if((fp_in = fopen(in_fname, "r")) == NULL) {
   fprintf(stderr, "extract_closest_table1/Fatal error opening input table %s\n",
           in_fname);
    return(-1);
  }

/* Open output table: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "extract_table1/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the output Latex table: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "%% Table of rho less than %f from: %s \n%% Created on %s", 
        in_fname, rho_c, ctime(&ttime));
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

fprintf(fp_out, "\\begin{table*}\n\
\\tabcolsep=1mm\n\
\\footnotesize\n\
\\centerline{\n\
\\begin{tabular}{cllclrr}\n\
\\hline\n\
\\bigstruttup\n\
WDS & Name & Epoch & $\\rho$ & Orbit & {\\scriptsize $\\Delta \\rho$(O-C)} & {\\scriptsize $\\Delta \\theta$(O-C)} \\\\\n\
& & & (\") & & (\") & \\\\\n\
\\hline\n\
\\bigstruttup");

iline = 0;
nobjects = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 256, fp_in)) {
    iline++;
// Remove all the non-printable characters and the end of line '\n'
// from input line:
      jlp_cleanup_string(in_line, 256);
#ifdef DEBUG
      printf("in_line=%s\n", in_line);
#endif

// Good lines start with a digit (WDS names...)
// Lines starting with % are ignored
    if(isdigit(in_line[0])) {

/* Get wds_name from column 1: */
     latex_get_column_item(in_line, wds_name, 1, verbose_if_error);

/* Get discov_name from column 2: */
     latex_get_column_item(in_line, discov_name, 2, verbose_if_error);

/* Get orbit_ref from column 11: */
     if(resid_only) {
       to_output = 0;
       verbose_if_error = 0;
       status = latex_get_column_item(in_line, orbit_ref, 11, verbose_if_error);
       if((status == 0) && (orbit_ref[0] != '\0')) {
#ifdef DEBUG
         printf("wds=%s discov=%s orbit ref: >%s<\n", 
                wds_name, discov_name, orbit_ref);
#endif
         to_output = 1;
         }
     } else {
// Output all lines if resid_only is false:
       to_output = 1;
     }

/* Get rho measure from column 5: */
      verbose_if_error = 0;
      rho_val = -1;
      drho_val = -1;
  if(latex_get_column_item(in_line, rho_meas, irho, verbose_if_error) == 0){
    jlp_clean_dollars(rho_meas, 40);
    if(strcmp(rho_meas, no_data)) {
      if(sscanf(rho_meas, "%lf", &dval) == 1) rho_val = dval;
      if(latex_get_column_item(in_line, drho_meas, 6, verbose_if_error) == 0){
       if(sscanf(drho_meas, "%lf", &dval) == 1) drho_val = dval;
      }
     }
     } else {
     rho_meas[0] = '\0';
     }

/* Get theta_obs from column 7: */
    theta_val = -1.;
    if(latex_get_column_item(in_line, theta_meas, 7, verbose_if_error) == 0){
      jlp_compact_string(theta_meas, 40);
      if(strcmp(theta_meas, no_data)) {
       if(sscanf(theta_meas, "%lf", &dval) == 1) theta_val = dval;
       }
      } else {
      theta_meas[0] = '\0';
      }

#ifdef DEBUG
    printf("wds_name=%s discov_name=%s \n", 
            wds_name, discov_name);
    printf("rho_val=%f theta_val=%f\n", rho_val, theta_val);
#endif

    if((to_output == 1) && (rho_val > 0.) && (theta_val > 0.) 
       && (rho_val < rho_c)) {
// SHOULD REMOVE THE LAST COLUMNS FIRST !!!!
// Remove "Notes" column:
       latex_remove_column(in_line, 10, 256);
// Remove "Dm" column:
       latex_remove_column(in_line, 9, 256);
// Remove "dtheta" column:
       latex_remove_column(in_line, 8, 256);
// Remove "theta" column:
//       latex_remove_column(in_line, 7, 256);
// Remove "drho" column:
       latex_remove_column(in_line, 6, 256);
// Remove "bin" column:
       latex_remove_column(in_line, 4, 256);
// Save to output file:
       fprintf(fp_out, "%s\n", in_line);
       nobjects++;
       }
// ! isdigit
    }
  } // fgets()
} // while

printf("nobjects=%d\n", nobjects);
fclose(fp_in);

fprintf(fp_out, "\\hline\n\
\\end{tabular}\n\
}\n\
\\end{table*}\n");
fprintf(fp_out, "%% nobjects=%d\n", nobjects);

fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input table and extract the table with the largest residuals 
* Example of header:
WDS & Name & Epoch & $\rho$ & $\sigma_\rho$ & Orbit & {\scriptsize $\Delta \rho$(O-C)} & {\
scriptsize $\Delta \theta$(O-C)} \\
*
* INPUT:
* in_fname: filename of the input table 
* rho_res_min : minimum rho threshold for selecting large residuals
* theta_res_min : minimum theta threshold for selecting large residuals
*
*************************************************************************/
static int extract_large_resid_table1(char *in_fname, char *out_fname, 
                                      double rho_res_min, double theta_res_min)
{
char in_line[256], buffer[256]; 
char rho_res_str[40], theta_res_str[40], orbit_ref[40];
char discov_name[40], wds_name[40], label[64]; 
int iline, status, verbose_if_error = 0, nobjects, orbit_grade;
double rho_res_val, theta_res_val, dval;
int i;

time_t ttime = time(NULL);
FILE *fp_out, *fp_in;

/* Open input table: */
if((fp_in = fopen(in_fname, "r")) == NULL) {
   fprintf(stderr, "extract_large_resid_table1/Fatal error opening input table %s\n",
           in_fname);
    return(-1);
  }

/* Open output table: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "extract_large_resid_table1/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the output Latex table: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "%% Table of large residuals: Delta_rho_min=%f Delta_theta_min=%f from: %s \n%% Created on %s \n", 
        rho_res_min, theta_res_min, in_fname, ctime(&ttime));
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

fprintf(fp_out, "\\begin{table*}\n\
\\tabcolsep=1mm\n\
\\footnotesize\n\
\\centerline{\n\
\\begin{tabular}{clrcccllrr}\n\
\\hline\n\
\\bigstruttup\n\
WDS & Name & Epoch & $\\rho$ & $\\sigma_\\rho$ & \\multicolumn{1}{c}{$\\theta$} & $\\sigma_\\theta$ & Orbit & {\\scriptsize $\\Delta \\rho$(O-C)} & {\\scriptsize $\\Delta \\theta$(O-C)} \\\\\n\
& &     &     & (\") & (\") & ($^\\circ$) & ($^\\circ$) & & \\\\\n\
\\hline\n\
\\bigstruttup");

iline = 0;
nobjects = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 256, fp_in)) {
    iline++;
// Remove all the non-printable characters and the end of line '\n'
// from input line:
      jlp_cleanup_string(in_line, 256);
#ifdef DEBUG
      printf("in_line=%s\n", in_line);
#endif

// Good lines start with a digit (WDS names...)
// Lines starting with % are ignored
    if(isdigit(in_line[0])) {

/* Get wds_name from column 1: */
     latex_get_column_item(in_line, wds_name, 1, verbose_if_error);

/* Get discov_name from column 2: */
     latex_get_column_item(in_line, discov_name, 2, verbose_if_error);

/* Get orbit_ref from column 11: */
     verbose_if_error = 0;
     orbit_ref[0] = '\0';
     status = latex_get_column_item(in_line, orbit_ref, 11, verbose_if_error);
     if((status == 0) && (orbit_ref[0] != '\0')) {
#ifdef DEBUG
       printf("orbit ref: >%s<\n", orbit_ref);
#endif

/* Get rho res from column 12: */
      verbose_if_error = 0;
      rho_res_val = 0;
      theta_res_val = 0;
      if(latex_get_column_item(in_line, rho_res_str, 12, verbose_if_error) == 0){
printf("ZZZAA: rho_res=%s\n", rho_res_str);
       jlp_clean_dollars(rho_res_str, 40);
       if(sscanf(rho_res_str, "%lf", &dval) == 1) rho_res_val = dval;
printf("ZZZA: rho_res=%s : value=%f\n", rho_res_str, dval);
       if(latex_get_column_item(in_line, theta_res_str, 13, verbose_if_error) == 0){
printf("ZZZAA: theta_res=%s\n", theta_res_str);
        jlp_clean_dollars(theta_res_str, 40);
        if(sscanf(theta_res_str, "%lf", &dval) == 1) theta_res_val = dval;
printf("ZZZA: theta_res=%s : value=%f\n", theta_res_str, dval);
       }
    if((ABS(rho_res_val) > rho_res_min) 
       || (ABS(theta_res_val) > theta_res_min)) {
    printf("ZZZADEBUG/ rho_res_min=%f theta_res_min=%f\n", 
               rho_res_min, theta_res_min);
    printf("DEBUG/ wds_name=%s discov_name=%s rho_res=%f theta_res=%f\n", 
            wds_name, discov_name, rho_res_val, theta_res_val);
#ifdef DEBUG
#endif
// SHOULD REMOVE THE LAST COLUMNS FIRST !!!!
// Remove "Notes" column:
       latex_remove_column(in_line, 10, 256);
// Remove "Dm" column:
       latex_remove_column(in_line, 9, 256);
// Remove "bin" column:
       latex_remove_column(in_line, 4, 256);
// Save to output file:
       fprintf(fp_out, "%s\n", in_line);
       nobjects++;
         } // if orbit_ref...
       }
    } // if latex
// ! isdigit
    }
  } // fgets()
} // while

printf("nobjects=%d\n", nobjects);
fclose(fp_in);

fprintf(fp_out, "\\hline\n\
\\end{tabular}\n\
}\n\
\\end{table*}\n");

fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input table and extract the table with the same WDS objects
* that have been measured twice or more times 
* that are discrepant (either same WDS name and different rho/theta if iopt=4,
* or same rho/theta and different discov name if iopt=5) 
* and extract corresp. table in both cases (iopt=4/5)
*
* INPUT:
*  in_fname: filename of the input table 
*  out_fname: filename of the output table 
*  iopt : 4, extract same WDS name and discrepant rho/theta ,
*         5, or same WDS name, small range in rho/theta and different discov name
*
*************************************************************************/
#define NOBS_MAX 20 
static int extract_twomeas_table1(char *in_fname, char *out_fname, int iopt)
{
char in_line[256], buffer[256]; 
char obs_line[NOBS_MAX][256], obs_discov[40*NOBS_MAX]; 
char discov_name[40], wds_name[40], label[64]; 
char last_wds_name[40], rho_meas[40], theta_meas[40], no_data[40]; 
double obs_rho[NOBS_MAX], dval, rho_val; 
double obs_theta[NOBS_MAX], theta_val; 
int iline, status, verbose_if_error = 0;
int i, iobs, nobs, out_result, same_discov;

time_t ttime = time(NULL);
FILE *fp_out, *fp_in;

/* Open input table: */
if((fp_in = fopen(in_fname, "r")) == NULL) {
   fprintf(stderr, "extract_twomeas_table1/Fatal error opening input table %s\n",
           in_fname);
    return(-1);
  }

/* Open output table: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "extract_twomeas_table1/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the output Latex table: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "%% Table of objects with two or more measures from: %s \n%% Created on %s", 
        in_fname, ctime(&ttime));
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

fprintf(fp_out, "\\begin{table*}\n\
\\tabcolsep=1mm\n\
\\small\n\
\\begin{tabular*}{\\textwidth}{clrcccccllllrr}\n\
\\hline\n\
& & & & & & & & & & & & \\\\\n\
WDS & Name & Epoch & Bin. & $\\rho$ & $\\sigma_\\rho$ & \\multicolumn{1}{c}{$\\theta$} & $\\sigma_\\theta$ & Dm & Notes & Orbit & {\\scriptsize $\\Delta \\rho$(O-C)} & {\\scriptsize $\\Delta \\theta$(O-C)} \\\\\n\
& &     &     & (\") & (\") & ($^\\circ$) & ($^\\circ$) & & \\\\\n\
& & & & & & & & & & & & \\\\\n\
\\hline\n\
& & & & & & & & & & & & \\\\\n");

strcpy(no_data, "\\nodata");
jlp_compact_string(no_data, 40);

iline = 0;
strcpy(last_wds_name, "000");
jlp_compact_string(last_wds_name, 40);

iobs = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 256, fp_in)) {
    iline++;
// Remove all the non-printable characters and the end of line '\n'
// from input line:
      jlp_cleanup_string(in_line, 256);
#ifdef DEBUG
      printf("in_line=%s\n", in_line);
#endif

// Good lines start with a digit (WDS names...)
// Lines starting with % are ignored
    if(isdigit(in_line[0])) {
/* Get wds_name from column 1: */
     latex_get_column_item(in_line, wds_name, 1, verbose_if_error);
// compact string:
     jlp_compact_string(wds_name, 40);

/* Get discov_name from column 2: */
     latex_get_column_item(in_line, discov_name, 2, verbose_if_error);

/* Get rho_obs from column 5: */
    rho_val = -1.;
    if(latex_get_column_item(in_line, rho_meas, 5, verbose_if_error) == 0){
      jlp_compact_string(rho_meas, 40);
      if(strcmp(rho_meas, no_data)) {
       if(sscanf(rho_meas, "%lf", &dval) == 1) {
        rho_val = dval;
        }
       }
      } else {
      rho_meas[0] = '\0';
      }

/* Get theta_obs from column 7: */
    theta_val = -1.;
    if(latex_get_column_item(in_line, theta_meas, 7, verbose_if_error) == 0){
      jlp_compact_string(theta_meas, 40);
      if(strcmp(theta_meas, no_data)) {
       if(sscanf(theta_meas, "%lf", &dval) == 1) {
        theta_val = dval;
        }
       }
      } else {
      theta_meas[0] = '\0';
      }

#ifdef DEBUG
    printf("wds_name=%s discov_name=%s last_wds_name=%s\n", 
            wds_name, discov_name, last_wds_name);
    printf("rho_val=%f theta_val=%f\n", rho_val, theta_val);
#endif
// Process the cases of same WDS name with observations:
    if(!strcmp(wds_name, last_wds_name) && (rho_val != -1.)
        && (theta_val != -1.)) {
        iobs++;
#ifdef DEBUG
      printf("EQUAL WDS names: obs=%d wds_name=%s< discov_name=%s last_wds_name=%s<\n", 
            iobs, wds_name, discov_name, last_wds_name);
#endif
        if(iobs >= NOBS_MAX-1 ) {
          fprintf(stderr,"extract_twomeas_table1/fatal error: iobs=%d > NOBS_MAX=%d\n",
                   iobs, NOBS_MAX);
          exit(-1);
          }  
// Save line to buffer: 
        strcpy(obs_line[iobs], in_line);
        strcpy(&obs_discov[40*iobs], discov_name);
        obs_rho[iobs] = rho_val;
        obs_theta[iobs] = theta_val;
#ifdef DEBUG
       printf("DEBUG/iopt=%d rho=%f theta=%f \n", iopt, 
            obs_rho[iobs], obs_theta[iobs]);
#endif
// When observation of another object, save all recorded lines to file: 
       } else {
         if (iobs > 0) {
           nobs = iobs + 1;
// iopt=4, same WDS name and discrepant rho or theta ,
// iopt=5, same WDS name, small range in rho and theta and different discov name
           out_result = 0;
           if(iopt == 4) {
            same_discov = 0;
            check_discrepant_rho_meas(obs_discov, obs_rho, nobs, 
                                      same_discov, &out_result);
            if(out_result == 0) 
              check_discrepant_theta_meas(obs_discov, obs_theta, nobs, 
                                          same_discov, &out_result);
           } else if(iopt == 5) { 
// Check if all the discov names of the series are the same 
// (out_result=0 if they are the same):
              same_discov = 0;
              check_same_discov_names(obs_discov, nobs, &out_result);
              if(out_result == 0) {
                 check_discrepant_rho_meas(obs_discov, obs_rho, nobs, 
                                          same_discov, &out_result);
                 if(out_result == 0) 
                   check_discrepant_theta_meas(obs_discov, obs_theta, nobs, 
                                               same_discov, &out_result);
              }
           }
// Save current line to output file:
           if(out_result != 0) {
/*
           printf("nobs=%d out_result=%d \n", nobs, out_result);
           printf("wds_name=%s< discov_name=%s\n", wds_name, discov_name);
           printf("obs_line=%s\n", obs_line[0]);
*/
             for(i = 0; i < nobs; i++)
                 fprintf(fp_out, "%s\n", obs_line[i]);
            }
         iobs = 0;
         }
       }
    strcpy(last_wds_name, wds_name);
    if(iobs == 0) {
       strcpy(obs_line[0], in_line);
       latex_get_column_item(in_line, discov_name, 2, verbose_if_error);
       strcpy(&obs_discov[0], discov_name);
       obs_rho[0] = -1;
       obs_theta[0] = -1;
       if((latex_get_column_item(in_line, rho_meas, 5, verbose_if_error) == 0)
       && (latex_get_column_item(in_line, theta_meas, 7, verbose_if_error) 
           == 0)){
        obs_rho[0] = rho_val;
        obs_theta[0] = theta_val;
       }
      }
    jlp_compact_string(last_wds_name, 40);
// ! isdigit
    }
  } // fgets()
} // while

fclose(fp_in);

fprintf(fp_out, "\\hline\n\
\\end{tabular*}\n\
Note: In column 7, the exponent $^*$ indicates that the position angle\n\
$\\theta$ could be determined without the 180$^\\circ$ ambiguity.\\\\\n\
 In column 7,  $!$ indicates that $\\theta$ could not be determined neither with this value, nor with WDS CHARA last measurement.\\\\\n\
\\end{table*}\n");

fclose(fp_out);
return(0);
}
/************************************************************************
* Check if all the discov names of the series are the same
*
* INPUT:
*   obs_discov[40*NOBS_MAX]
*
* OUTPUT:
*   out_result: 
*               0, not discrepant (same discov name)
*               1, different discov name in the series
*************************************************************************/
static int check_same_discov_names(char *obs_discov, int nobs, int *out_result)
{
int k;
char discov0[40];
char discov1[40];

 *out_result = 0;

// Check if the discov names of the series are different:
  strcpy(discov0, &obs_discov[0]);
  jlp_compact_string(discov0, 40);
  for(k = 0; k < nobs; k++) {
     strcpy(discov1, &obs_discov[k*40]);
     jlp_compact_string(discov1, 40);
     if(strcmp(discov0, discov1)) {
        *out_result = 1;
        return(0);
       } 
     }

return(0);
}
/************************************************************************
* Look for discrepant measurements for the same WDS objects 
* that have been measured twice or more times
* Same name and different rho and theta,
*
* INPUT:
*   obs_discov[40*NOBS_MAX]
*   same_discov: flag set to one if calculation has to be done for the same discover name
*
* OUTPUT:
*   out_result: 
*               0, not discrepant (same WDS name and compatible rho measurements)
*               1, same WDS name and not compatible rho,
* 
*************************************************************************/
static int check_discrepant_rho_meas(char *obs_discov, double *obs_rho,
                                     int nobs, int same_discov, int *out_result)
{
int i, k, nobs0;
double mean_rho, error0; 
char discov0[40];

 *out_result = 0;

if(nobs <= 2) return(0);

// Loop on the names:
for(k = 0; k < nobs; k++) {
  strcpy(discov0, &obs_discov[k*40]);

// Look for discrepant rho:
// First compute the mean:
 mean_rho = 0.;
 nobs0 = 0;
   for(i = 0; i < nobs; i++) {
      if(!strcmp(discov0, &obs_discov[i*40]) || (same_discov == 0)) {
         mean_rho += obs_rho[i];
         nobs0++;
         } // strcmp
     }
 if(nobs0 > 1) {
   mean_rho /= (double)nobs0;
//   printf("mean_rho=%f, nobs0=%d\n", mean_rho, nobs0);

// Discrepant rho if relative error larger than 20%
// ABS(rho - mean_rho)/mean_rho > 0.2:
   for(i = 0; i < nobs; i++) {
      if(!strcmp(discov0, &obs_discov[i*40])) {
         error0 = ABS(obs_rho[i] - mean_rho) / (ABS(mean_rho) + 0.0001);
         if(error0 > 0.2) {
           *out_result = 1;
           return(0);
         }
       } // strcmp
    }
 } // if nobs0 > 0

} // EOF loop on k 

return(0);
}
/************************************************************************
* Look for discrepant theta and for WDS objects 
* that have been measured twice or more times
*
* INPUT:
*   obs_discov[40*NOBS_MAX]
*   same_discov: flag set to one if calculation has to be done for the same discover name
*
* OUTPUT:
*   out_result: 
*               0, not discrepant (same name and compatible theta measurements)
*               1, compatible theta and different name
* 
*************************************************************************/
static int check_discrepant_theta_meas(char *obs_discov, double *obs_theta,
                                       int nobs, int same_discov, 
                                       int *out_result)
{
int i, k;
double min_theta, max_theta, scatter0, ww; 
char discov0[40];

*out_result = 0;

if(nobs <= 2) return(0);

strcpy(discov0, &obs_discov[0]);

// Loop on the measurements:
// Compute the mean with all objects: 
min_theta = 1000.;
max_theta = -1000.;
 for(i = 0; i < nobs; i++) {
      if(!strcmp(discov0, &obs_discov[i*40]) || (same_discov == 0)) {
         if(obs_theta[i] > 190.) ww = obs_theta[i] - 180.;
         else ww = obs_theta[i];
         if(ww < min_theta) min_theta = ww;
         if(ww > max_theta) max_theta = ww;
      }
   }
scatter0 = max_theta - min_theta;
/*
printf("min_theta=%f, max_theta=%f scatter0=%f nobs=%d\n", 
        min_theta, max_theta, scatter0, nobs);
*/

// Discrepant theta if scatter larger than 10 degr. 
  if(scatter0 > 10.) {
      *out_result = 1;
      return(0);
    } 

return(0);
}
/************************************************************************
* Look for unsorted measurements (according to epoch) 
*
* INPUT:
*   obs_discov[40*nobs]
*   obs_epoch[nobs]
*   same_discov: flag set to one if calculation has to be done for the same discover name
*
* OUTPUT:
*   out_result: 
*               0, not discrepant (epoch-sorted measurements)
*               1, epoch-unsorted measurements 
* 
*************************************************************************/
static int check_sorted_epoch_meas(char *obs_discov, double *obs_epoch,
                                   int nobs, int same_discov, int *out_result)

{
int i, isorted = 1;
double ww, previous_epoch;
char discov0[40];

*out_result = 0;

if(nobs <= 2) return(0);

strcpy(discov0, &obs_discov[0]);

// Loop on the measurements:
isorted = 1;
previous_epoch = -1;
 for(i = 0; i < nobs; i++) {
/*
printf("i=%d nobs=%d same_discov=%d\n", i, nobs, same_discov);
*/
      if(!strcmp(discov0, &obs_discov[i*40]) || (same_discov == 0)) {
         if(obs_epoch[i] < previous_epoch) {
           isorted = 0;
           break;
         }
       previous_epoch = obs_epoch[i];
      }
   }

// Discrepant if unsorted 
  if(isorted == 0) {
      *out_result = 1;
    } 

return(0);
}
/************************************************************************
* Scan the input table and extract the table with the objects
* with odd names (companions starting with odd letters)
*
* INPUT:
*  in_fname: filename of the input table 
*  out_fname: filename of the output table 
*
*************************************************************************/
static int extract_oddnames_table1(char *in_fname, 
                                   char *out_fname)
{
char in_line[256], buffer[256]; 
char discov_name[40], wds_name[40], label[64]; 
int iline, status, verbose_if_error = 0;
int i, name_is_odd;

time_t ttime = time(NULL);
FILE *fp_out, *fp_in;

/* Open input table: */
if((fp_in = fopen(in_fname, "r")) == NULL) {
   fprintf(stderr, "extract_oddnames_table1/Fatal error opening input table %s\n",
           in_fname);
    return(-1);
  }

/* Open output table: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "extract_oddnames_table1/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the output Latex table: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "%% Table of objects with odd names from: %s \n%% Created on %s", 
        in_fname, ctime(&ttime));
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

iline = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 256, fp_in)) {
    iline++;
// Remove all the non-printable characters and the end of line '\n'
// from input line:
      jlp_cleanup_string(in_line, 256);
#ifdef DEBUG
      printf("in_line=%s\n", in_line);
#endif

// Good lines start with a digit (WDS names...)
// Lines starting with % are ignored
    if(isdigit(in_line[0])) {
/* Get wds_name from column 1: */
     latex_get_column_item(in_line, wds_name, 1, verbose_if_error);
/* Get discov_name from column 2: */
     latex_get_column_item(in_line, discov_name, 2, verbose_if_error);
     jlp_compact_string(discov_name, 40);
     check_if_oddname(discov_name, 40, &name_is_odd);

#ifdef DEBUG
    printf("wds_name=%s discov_name=%s \n", 
            wds_name, discov_name);
#endif
// Save line to output file: 
     if(name_is_odd) fprintf(fp_out, "%s\n", in_line);
// ! isdigit
    }
  } // fgets()
} // while

fclose(fp_in);
fclose(fp_out);
return(0);
}
/*************************************************************************
* Check if input discover's name is odd
*
**************************************************************************/
static int check_if_oddname(char *discov_name0, int nlength0, int *name_is_odd)
{
char extens[40], *pc;

*name_is_odd = 0;
jlp_compact_string(discov_name0, nlength0);

// Look for extension if present:
pc = discov_name0;
// First scan the starting letters:
while(*pc && isalpha(*pc)) pc++;

// Then scan the numbers:
while(*pc && isdigit(*pc)) pc++;

strcpy(extens, pc);
/*****
// To upper case:
 pc = extens;
 while(*pc) {
   *pc = toupper((unsigned char)*pc);
   pc++;
   }
*****/

if(extens[0] != '\0') {
//   printf("check_if_oddname/discov_name0=%s< extens=%s< \n", 
//           discov_name0, extens);
   if(strcmp(extens, "AB") && strcmp(extens, "AB,C") && strcmp(extens, "DB")
    && strcmp(extens, "A,BC") && strcmp(extens, "Aa") && strcmp(extens, "AB,CD")
    && strcmp(extens, "Aa,Ab") && strcmp(extens, "Aa,Ac") 
    && strcmp(extens, "BC") && strcmp(extens, "CD")
    && strcmp(extens, "BD") 
    && strcmp(extens, "AC") && strcmp(extens, "Ba,Bb") 
    && strcmp(extens, "Ca,Cb"))
     {
      *name_is_odd = 1;
      printf("check_if_oddname/discov_name0=%s extens=%s< name_is_odd=%d \n", 
              discov_name0, extens, *name_is_odd); 
     } 
  }

return(0);
}
/************************************************************************
* Scan the input table and extract the table with the objects
* with odd names (companions starting with odd letters)
*
* INPUT:
*  in_fname: filename of the input table 
*  out_fname: filename of the output table 
*
*************************************************************************/
static int extract_new_doubles_table1(char *in_fname, char *out_fname)
{
char in_line[256], buffer[256]; 
char discov_name[40], wds_name[40], label[64]; 
int iline, status, verbose_if_error = 0;
int i, is_new_double;

time_t ttime = time(NULL);
FILE *fp_out, *fp_in;

/* Open input table: */
if((fp_in = fopen(in_fname, "r")) == NULL) {
   fprintf(stderr, "extract_new_doubles_table1/Fatal error opening input table %s\n",
           in_fname);
    return(-1);
  }

/* Open output table: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "extract_new_doubles_table1/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the output Latex table: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "%% Table of objects with odd names from: %s \n%% Created on %s", 
        in_fname, ctime(&ttime));
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

iline = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 256, fp_in)) {
    iline++;
// Remove all the non-printable characters and the end of line '\n'
// from input line:
      jlp_cleanup_string(in_line, 256);
#ifdef DEBUG
//      printf("in_line=%s\n", in_line);
#endif

// Good lines start with a digit (WDS names...)
// Lines starting with % are ignored
    if(isdigit(in_line[0])) {
/* Get wds_name from column 1: */
     latex_get_column_item(in_line, wds_name, 1, verbose_if_error);
/* Get discov_name from column 2: */
     latex_get_column_item(in_line, discov_name, 2, verbose_if_error);
     jlp_compact_string(discov_name, 40);
//     check_if_new_double(discov_name, 40, &is_new_double);
// Search for string in string:
     if((strstr(in_line, "nd") != NULL) ||
        (strstr(in_line, "ND") != NULL)) {
#ifdef DEBUG
    printf("in_line=%s\n", in_line);
#endif
// Save line to output file: 
          fprintf(fp_out, "%s\n", in_line);
        }

// ! isdigit
    }
  } // fgets()
} // while

fclose(fp_in);
fclose(fp_out);
return(0);
}
/*************************************************************************
* Check if new double (ND or nd is found) 
*
**************************************************************************/
static int check_if_new_double(char *discov_name0, int nlength0, 
                               int *is_new_double)
{
char extens[40], *pc;

*is_new_double = 0;
jlp_compact_string(discov_name0, nlength0);

// Look for extension if present:
pc = discov_name0;
// First scan the starting letters:
while(*pc && isalpha(*pc)) pc++;

// Then scan the numbers:
while(*pc && isdigit(*pc)) pc++;

strcpy(extens, pc);
/*****
// To upper case:
 pc = extens;
 while(*pc) {
   *pc = toupper((unsigned char)*pc);
   pc++;
   }
*****/

if(extens[0] != '\0') {
//   printf("check_if_new_double/discov_name0=%s< extens=%s< \n", 
//           discov_name0, extens);
   if(strcmp(extens, "AB") && strcmp(extens, "AB,C") && strcmp(extens, "DB")
    && strcmp(extens, "A,BC") && strcmp(extens, "Aa") && strcmp(extens, "AB,CD")
    && strcmp(extens, "Aa,Ab") && strcmp(extens, "Aa,Ac") 
    && strcmp(extens, "BC") && strcmp(extens, "CD")
    && strcmp(extens, "AC") && strcmp(extens, "Ba,Bb") 
    && strcmp(extens, "Ca,Cb"))
     {
      *is_new_double = 1;
      printf("check_if_new_double/discov_name0=%s extens=%s< is_new_double=%d \n", 
              discov_name0, extens, *is_new_double); 
     } 
  }

return(0);
}
/************************************************************************
* Sort the input table accrding to the epoch for with the same discov objects
* that have been measured twice or more times 
*
* INPUT:
*  in_fname: filename of the input table 
*  out_fname: filename of the output table 
*
*************************************************************************/
static int sort_measures_table1(char *in_fname, char *out_fname)
{
char in_line[256], buffer[256]; 
char obs_line[NOBS_MAX][256], obs_discov[40*NOBS_MAX]; 
char discov_name[40], wds_name[40], label[64]; 
char last_wds_name[40], epoch_meas[40], no_data[40]; 
double obs_epoch[NOBS_MAX], dval, epoch_val; 
int iline, status, verbose_if_error = 0;
int i, iobs, nobs, out_result, same_discov;

time_t ttime = time(NULL);
FILE *fp_out, *fp_in;

/* Open input table: */
if((fp_in = fopen(in_fname, "r")) == NULL) {
   fprintf(stderr, "extract_twomeas_table1/Fatal error opening input table %s\n",
           in_fname);
    return(-1);
  }

/* Open output table: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "_table1/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the output Latex table: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "%% Table of measures sorted from: %s \n%% Created on %s", 
        in_fname, ctime(&ttime));
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

strcpy(no_data, "\\nodata");
jlp_compact_string(no_data, 40);

iline = 0;
strcpy(last_wds_name, "000");
jlp_compact_string(last_wds_name, 40);

iobs = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 256, fp_in)) {
    iline++;
// Remove all the non-printable characters and the end of line '\n'
// from input line:
      jlp_cleanup_string(in_line, 256);
#ifdef DEBUG1
      printf("sort_measures/in_line=%s\n", in_line);
#endif

// Good lines start with a digit (WDS names...)
// Lines starting with % are ignored
    if(isdigit(in_line[0])) {
/* Get wds_name from column 1: */
     latex_get_column_item(in_line, wds_name, 1, verbose_if_error);
// compact string:
     jlp_compact_string(wds_name, 40);

/* Get discov_name from column 2: */
     latex_get_column_item(in_line, discov_name, 2, verbose_if_error);

/* Get epoch from column 3: */
    epoch_val = -1.;
    if(latex_get_column_item(in_line, epoch_meas, 3, verbose_if_error) == 0){
      jlp_compact_string(epoch_meas, 40);
       if(sscanf(epoch_meas, "%lf", &dval) == 1) {
        epoch_val = dval;
        }
      } else {
      epoch_meas[0] = '\0';
      } // if latex

#ifdef DEBUG1
    printf("wds_name=%s discov_name=%s last_wds_name=%s\n", 
            wds_name, discov_name, last_wds_name);
#endif
// Process the cases of same WDS name with observations:
    if(!strcmp(wds_name, last_wds_name) && (epoch_val != -1.)) {
        iobs++;
#ifdef DEBUG
      printf("EQUAL WDS names: obs=%d wds_name=%s< discov_name=%s last_wds_name=%s<\n", 
            iobs, wds_name, discov_name, last_wds_name);
    printf("epoch_val=%f\n", epoch_val);
#endif
        if(iobs >= NOBS_MAX-1 ) {
          fprintf(stderr,"extract_twomeas_table1/fatal error: iobs=%d > NOBS_MAX=%d\n",
                   iobs, NOBS_MAX);
          exit(-1);
          }  
// Save line to buffer: 
        strcpy(obs_line[iobs], in_line);
        strcpy(&obs_discov[40*iobs], discov_name);
        obs_epoch[iobs] = epoch_val;
// When observation of another object, save all recorded lines to file: 
       } else {
         if (iobs > 0) {
           nobs = iobs + 1;
#ifdef DEBUG
          printf("Processing all the measures of the same object, nobs=%d\n", nobs);
#endif
           out_result = 0;
           same_discov = 0;
           check_sorted_epoch_meas(obs_discov, obs_epoch, nobs, 
                                   same_discov, &out_result);
// Save current line to output file:
           if(out_result != 0) {
             printf("nobs=%d out_result=%d \n", nobs, out_result);
             for(i = 0; i < nobs; i++) {
               printf("obs_line[%d]=%s obs_epoch=%f\n", 
                       i, obs_line[i], obs_epoch[i]);
               fprintf(fp_out, "%s\n", obs_line[i]);
               }
            }
         iobs = 0;
         } else {
// if not the same wds name, simply copy the recorded line
//         fprintf(fp_out, "%s\n", obs_line[0]);
         }
       } // EOF epoch_val != -1 and same wds_name
    strcpy(last_wds_name, wds_name);
    if(iobs == 0) {
       strcpy(obs_line[0], in_line);
       latex_get_column_item(in_line, discov_name, 2, verbose_if_error);
       strcpy(&obs_discov[0], discov_name);
       obs_epoch[0] = -1;
       if(latex_get_column_item(in_line, epoch_meas, 3, verbose_if_error) == 0){
        obs_epoch[0] = epoch_val;
       }
     } // if iobs==0
    jlp_compact_string(last_wds_name, 40);
    } else {     // EOF isdigit
// if not digit, simply copy this line
//      fprintf(fp_out, "%s\n", in_line);
    }
  } // fgets()
} // while

fclose(fp_in);

fclose(fp_out);
return(0);
}

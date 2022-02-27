/************************************************************************
* "process_table1_with WDS.cpp"
*
* Process "full_tab_calib.tex"
* - iop=1,delta_rho_rho_max,delta_theta_max generate table with the latest WDS (rho,theta) 
*   measurements and comparison in the last columns, with delta_rho_rho_max,dtheta thresholds
*   for output 
* - iop=2,spmin,spmax generate table with all the WDS spectroscopic 
*   classification in the last column, ans spmin,spmax for output 
*   Ex: iop=2,m1V,k7V will generate the table with the white darfs  
*
* JLP 
* Version 10/05/2021
*************************************************************************/
#include <stdio.h> 
#include <stdlib.h> // exit(-1) 
#include <string.h> 
#include <ctype.h>  // isalpha(), isdigit() 
#include <time.h>   // time_t 
#include "latex_utils.h" // latex_get_column_item(), latex_remove_column()
#include "jlp_string.h"  // (in jlplib/jlp_fits/ ) jlp_cleanup_..
#include "jlp_catalog_utils.h"  // ABS 
#include "WDS_catalog_utils.h"  // get_data_from_WDS_catalog()
#include "astrom_utils1.h" // astrom_read_object_data_for_name_only()

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
double rho_res;
double theta_res;
} PSC_MEAS;

// Contained here:
static int count_all_objects(char *in_fname);
static int add_wds_meas_to_table1(char *in_fname, char *out_fname, 
                                  char *wds_cat_fname, double delta_rho_rho_max,
                                  double delta_theta_max);
static int extract_white_dwarfs_from__table1(char *in_fname, char *out_fname,
                                             char *wds_cat_fname, char *sp_mini,
                                             char *sp_maxi);
static int check_if_new_double(char *discov_name0, int nlength0, 
                               int *is_new_double);
static int init_psc_meas(PSC_MEAS *psc0, char* wds0, char* dsc0, double epoch0,
                         double eyep0, double rmeas0, 
                         double drmeas0, double tmeas0, double dtmeas0, 
                         double dm0, double rres0, double tres0);
static int printf_psc_meas(PSC_MEAS psc0, char *label);

static int jlp_clean_dollars(char *str0, int str_len0);
/*************************************************************************
*
*************************************************************************/
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

// Residuals:
psc0->rho_res = rres0;
psc0->theta_res = tres0;

// Delta_mag:
psc0->dm_meas = dm0;

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
return(0);
}


/***********************************************************************
*
* - iop=1,delta_rho_rho_max,delta_theta_max generate table with the latest WDS (rho,theta) 
*   measurements and comparison in the last columns, with delta_rho_rho,dtheta thresholds
*   for output 
* - iop=2,spmin,spmax generate table with all the WDS spectroscopic 
*   classification in the last column, ans spmin,spmax for output 
*   Ex: iop=2,m1V,k7V will generate the table with the white darfs  
************************************************************************/
int main(int argc, char *argv[])
{
char in_fname[128], out_fname[128], wds_cat_fname[128];
char sp_mini[64], sp_maxi[64];
double delta_rho_rho_max, delta_theta_max;
int iopt, resid_only;

if(argc != 5) {
  printf("Syntax:\n"); 
  printf("Option1: latest rho,theta WDS meas. and comp. in the last cols. \n");
  printf(" 1,delta_rho_rho_max,delta_theta_max for output \n");
  printf("process_calib_table_with_WDS in_table out_table wds_catalog 1,0.2,20.\n");
  return(-1);
}
strcpy(in_fname, argv[1]);
strcpy(out_fname, argv[2]);
strcpy(wds_cat_fname, argv[3]);
sscanf(argv[4], "%d", &iopt);
if(iopt == 1)
   sscanf(argv[4], "%d,%lf,%lf", &iopt, &delta_rho_rho_max, &delta_theta_max);
else 
   sscanf(argv[4], "%d,%s,%s", &iopt, &sp_mini, &sp_maxi);

printf("OK: in_table=%s output=%s wds_cat=%s\n", 
        in_fname, out_fname, wds_cat_fname); 
if(iopt == 1) {
  printf("OK: iopt=%d delta_rho_rho_max=%f delta_theta_max=%f\n", 
         iopt, delta_rho_rho_max, delta_theta_max); 
} else { 
  printf("OK: iopt=%d sp_mini=%s sp_maxi=%s\n", 
         iopt, sp_mini, sp_maxi); 
}

count_all_objects(in_fname);
switch(iopt) {
  case 1:
  default:
     add_wds_meas_to_table1(in_fname, out_fname, wds_cat_fname, delta_rho_rho_max, 
                           delta_theta_max); 
     break;
  case 2:
     extract_white_dwarfs_from__table1(in_fname, out_fname, wds_cat_fname,
                                       sp_mini, sp_maxi);
     break;
  }
return(0);
}
/************************************************************************
* Scan the input table and make the modifications 
* - count de number of objects, the number of observations
*
* INPUT:
* in_fname: input table filename
*
*************************************************************************/
static int count_all_objects(char *in_fname)
{
char in_line[256], in_line2[256], in_line3[256]; 
char discov_name[40], wds_name[40], label[64], old_discov_name[40]; 
char buffer[256], rho_meas[40], theta_meas[40], orbit_ref[40], *pc;
int iline, status, verbose_if_error = 0, nobj, nmeas, nunres, nresid, ndmag;
int nlines_meas, nquad;
double rho_val, drho_val, theta_val, dtheta_val, dmag_val, dval;
char drho_meas[40], dtheta_meas[40], dmag_str[40];
int i, irho = 5, i_drho = 6, resolved;
FILE *fp_in;
time_t ttime = time(NULL);

/* Open input table: */
if((fp_in = fopen(in_fname, "r")) == NULL) {
   fprintf(stderr, "count_objects/Error opening input table %s\n",
           in_fname);
    return(-1);
  }

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
/*
       if((rho_val > 0) && (drho_val > 0.)) {
int init_psc_meas(PSC_MEAS *psc0, char* wds0, char* dsc0, double epoch0,
                  double eyep0, double rmeas0,
                  double drmeas0, double tmeas0, double dtmeas0,
                  double dm0, double rres0, double tres0);
         init_psc_meas(&psc_rho_min, wds_name, discov_name, 0., 0.,
                       rho_val, drho_val, theta_val, dtheta_val, 0., 0., 0.);
      }
*/

    }/* EOF if !isdigit % ... */

  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("count_objects: %d lines sucessfully read and processed, nlines_meas=%d\n", 
        iline, nlines_meas);
printf("count_objects: n_objects=%d n_meas(resolved)=%d n_unres=%d \n", 
        nobj, nmeas, nunres);
printf("count_objects: n_resid=%d n_mag=%d n_quad=%d\n", 
        nresid, ndmag, nquad);

/* Close opened files:
*/
fclose(fp_in);

return(0);
}
/************************************************************************
* Scan the calib file and add various parameters
* add WDS rho, theta measurements and extract objects for which the
* difference between maesurements and WDS measurements 
* are smaller than (delta_rho_rho_max, delta_theta_max)
*
* From WDS catalog, retrieve:
* WDS number,
* WY, WT, WR (year, theta and rho of the last observation)
*
* INPUT:
* in_fname: input table filename
* out_fname: output Latex filename
* wds_cat_fname: WDS catalog filename
*  e.g., wdsweb_summ.txt (used to obtain the last WDS observations)
*
*************************************************************************/
static int add_wds_meas_to_table1(char *in_fname, char *out_fname, 
                                  char *wds_cat_fname, double delta_rho_rho_max,
                                  double delta_theta_max) 
{
char in_line[256], in_line2[256], in_line3[256]; 
char full_discov_name[64], discov_name[64], wds_name[64], label[64]; 
char buffer[256], rho_meas[64], theta_meas[64], *pc;
char comp_name[64], spectral_type[64], wds_name_from_WDS_catalog[64]; 
double WdsLastYear, WdsLastTheta, WdsLastRho, mag_A, mag_B;
int iline, status, wds_meas_found, verbose_if_error = 0;
double rho_val, theta_val, dval, delta_rho_rho, delta_theta;
int i;
FILE *fp_in, *fp_out;
time_t ttime = time(NULL);

/* Open input table: */
if((fp_in = fopen(in_fname, "r")) == NULL) {
 fprintf(stderr, "add_wds_meas_to_table1/Fatal error opening input table %s\n",
         in_fname);
 return(-1);
 }

/* Open output table: */
if((fp_out = fopen(out_fname, "w")) == NULL) {
 fprintf(stderr, "add_wds_meas_to_table1/Fatal error opening output file: %s\n",
         out_fname);
 return(-1);
 }

/* Header of the output Latex table: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "Routine add_wds_meas_to_table1 using %s\n", wds_cat_fname);
fprintf(fp_out, "%% Table of objects from %s with delta_rho_rho,dtheta smaller than %f,%f \n%% Created on %s", 
        in_fname, delta_rho_rho_max, delta_theta_max, ctime(&ttime));
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
     latex_get_column_item(in_line, full_discov_name, 2, verbose_if_error);

// Split discov_name and comp_name:
      jlp_split_discov_comp(full_discov_name, 64, discov_name, comp_name);

/* Remove the extra-blanks: */
     jlp_trim_string(discov_name, 64);
     jlp_trim_string(comp_name, 64);

#ifdef DEBUG
printf("From calib table: discov_name=%s< comp_name=%s< wds_name=>%s< \n",
       discov_name, comp_name, wds_name);
#endif
// Change some companion names:
  if(!strcmp(comp_name, "ab")) strcpy(comp_name, "AB");
  if(!strcmp(comp_name, "bc")) strcpy(comp_name, "BC");
  if(!strcmp(comp_name, "cd")) strcpy(comp_name, "CD");
  if(!strcmp(comp_name, "ac")) strcpy(comp_name, "AC");
  if(!strcmp(comp_name, "ab-c")) strcpy(comp_name, "AB-C");

/*
* Get miscellaneous data from the WDS catalog from discov_name and comp_name
int get_data_from_WDS_catalog(char *WDS_catalog, char *discov_name,
                              char *comp_name, char *wds_name_from_WDS_catalog,
                              double *WdsLastYear,
                              double *WdsLastRho, double *WdsLastTheta,
                              double *WdsMagA, double *WdsMagB,
                              char *WdsSpectralType,
                              int *wds_meas_found);
*/
     get_data_from_WDS_catalog(wds_cat_fname, discov_name, comp_name, 
                               wds_name_from_WDS_catalog,
                               &WdsLastYear, &WdsLastRho, &WdsLastTheta,
                               &mag_A, &mag_B, spectral_type,
                               &wds_meas_found);
#ifdef DEBUG
printf("get_data_from_WDS catalog: wds_name_from_WDS_catalog=>%s< \n",
        wds_name_from_WDS_catalog);
printf("%s = %s%s WY=%.2f WR=%.2f WT=%.1f \n",
       wds_name, discov_name, comp_name, WdsLastYear, WdsLastRho, WdsLastTheta);
#endif


/* Get rho measure from column irho=5: */
    delta_rho_rho = 12345.;
    if((wds_meas_found == 1) && (WdsLastRho > 0.01)) {
      verbose_if_error = 0;
      rho_val = -1;

      if(latex_get_column_item(in_line, rho_meas, 5, verbose_if_error) == 0){
        if(sscanf(rho_meas, "%lf", &dval) == 1) {
          rho_val = dval;
          delta_rho_rho = (rho_val - WdsLastRho) / WdsLastRho;
          }
      }
    } // EOF if wds_meas_found == 1

/* Get theta measure from column itheta=7: */
    delta_theta = 12345.;
    if(wds_meas_found == 1) {
      verbose_if_error = 0;
      theta_val = -1;

      if(latex_get_column_item(in_line, theta_meas, 7, verbose_if_error) == 0){
        if(sscanf(theta_meas, "%lf", &dval) == 1) {
          theta_val = dval;
          delta_theta = theta_val - WdsLastTheta;
          }
      }
    } // EOF if wds_meas_found == 1

      if(((delta_rho_rho != 12345.) && (ABS(delta_rho_rho) > delta_rho_rho_max)) 
       || ((delta_theta != 12345.) && (ABS(delta_theta) > delta_theta_max))){ 
          printf("\nFrom calib table: %s %s %s rho=%.3f theta=%.1f\n",
          discov_name, comp_name, wds_name, rho_val, theta_val);
          printf("From WDS catalog: wds=>%s< WY=%.2f WR=%.2f WT=%.1f\n",
                  wds_name_from_WDS_catalog, WdsLastYear, WdsLastRho,
                  WdsLastTheta);
          if((delta_rho_rho != 12345.) && (delta_theta != 12345.)) {
             printf("delta_rho_rho=%f delta_theta=%f \n", delta_rho_rho, delta_theta); 
// Save to output file:
            fprintf(fp_out, "%s WY=%.2f WR=%.2f WT=%.1f drho2=%f dth=%f\n", 
                    in_line3, WdsLastYear, WdsLastRho, WdsLastTheta, 
                    delta_rho_rho, delta_theta);
            } else if(delta_rho_rho != 12345.) { 
             printf("delta_rho_rho=%f \n", delta_rho_rho); 
            } else if(delta_theta != 12345.) { 
             printf("delta_theta=%f \n", delta_theta); 
            }
        } // EOF delta_rho_rho != 12345

    }/* EOF if !isdigit % ... */


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
* Scan the input table and extract the table with the objects
* whose spectra are in the range [sp_mini, sp_maxi]
*
* INPUT:
*  in_fname: filename of the input table 
*  out_fname: filename of the output table 
*  wds_cat_fname: WDS catalog filename
*  sp_mini, sp_maxi: spectral types from WDS classification
*
*************************************************************************/
static int extract_white_dwarfs_from__table1(char *in_fname, char *out_fname,
                                             char *wds_cat_fname, char *sp_mini,
                                             char *sp_maxi)
{
char in_line[256], buffer[256]; 
char discov_name[40], wds_name[40], label[64]; 
int iline, status, verbose_if_error = 0;
int i, is_new_double;

time_t ttime = time(NULL);
FILE *fp_out, *fp_in;

/* Open input table: */
if((fp_in = fopen(in_fname, "r")) == NULL) {
 fprintf(stderr, "extract_white_dwarfs_from__table1/Fatal error opening input table %s\n",
          in_fname);
 return(-1);
 }

/* Open output table: */
if((fp_out = fopen(out_fname, "w")) == NULL) {
 fprintf(stderr, "extract_white_dwarfs_from__table1/Fatal error opening output file: %s\n",
           out_fname);
 return(-1);
 }

/* Header of the output Latex table: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "Routine add_wds_meas_to_table1 using %s\n", wds_cat_fname);
fprintf(fp_out, "%% Table of objects from %s with spectral types from: %s to %s \n%% Created on %s", 
        in_fname, sp_mini, sp_maxi, ctime(&ttime));
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

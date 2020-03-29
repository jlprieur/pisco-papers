/*************************************************************************
* astrom_utils2.c
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

#define DEBUG
/*
*/
int astrom_write_publi_table(FILE *fp_out, int comments_wanted, OBJECT *obj, 
                             int *index_obj, int nobj, int tabular_only);
int astrom_write_publi_table_gili(FILE *fp_out, int comments_wanted, 
                                  OBJECT *obj, int *index_obj, int nobj, 
                                  int tabular_only);

/*************************************************************************
* Scan the astrom file and add the epoch read from the header 
* of FITS autocorrelation files 
*
* INPUT:
* fits_directory: directory containing the FITS files 
*                (e.g., /home/data/pisco_merate/2004-2008/ )
*
**************************************************************************/
int astrom_add_epoch_from_fits_file(FILE *fp_in, FILE *fp_out, 
                                    char *fits_directory)
{
int iline, is_measurement, epoch_was_found;
int eyepiece0, eyepiece1;
double epoch0;
char b_in[NMAX], fits_filename[100], *pc, date0[20], date1[20];
char full_directory[100];

iline = 0;
while(!feof(fp_in))
{
/* Maximum length for a line will be 170 characters: */
  if(fgets(b_in,170,fp_in))
  {
  iline++;
  b_in[169] = '\0';
/* Remove ^M (Carriage Return) if present: */
  pc = b_in;
  while(*pc) {
  if(*pc == '\r') *pc = ' ';
  pc++;
  }
 
    if(b_in[0] == '&') { 
/* Check if input LateX line contains a measurement: 
*/
        astrom_check_if_measurement(b_in, fits_filename, date1, &eyepiece1,
                                    &is_measurement);
        if(is_measurement) {
#ifdef DEBUG
           printf("DDEBUG/measure: %s", b_in);
           printf("DDEBUG: fits_file >%s< (date=%s eyepiece=%d) \n", 
                  fits_filename, date1, eyepiece1);
#endif
           merate_get_full_directory(fits_directory, date1, full_directory);

/* In "FITS_utils.c" */
           get_bessel_epoch_from_fits_file(fits_filename, full_directory, 
                                           &epoch0, date0, &epoch_was_found);
 printf("OK: epoch_was_found=%d (bessel epoch0=%f)\n", epoch_was_found, epoch0);
           if(epoch_was_found) { 
/* Remove "\\" (EOF line for Latex tables) if present: */
            pc = b_in;
            while(*pc) {
            if(!strncmp(pc,"\\\\",2)) break;
            pc++;
            }
            *pc = '\0';
/* Add epoch to truncated line: */
            sprintf(b_in, "%s EP=%.4f \\\\\n", b_in, epoch0);
 printf("OK: >%s<\n", b_in);
            }
         } /* EOF is_measurement */
    }  /* EOF case of not commented line*/

// Copy current line to file:
        fputs(b_in, fp_out);
  } /* EOF if fgets() */
} /* EOF while loop */

return(0);
}
/*************************************************************************
* Scan the file and add various parameters 
* WDS number, discoverer's name and "orbit" from PISCO catalog
* WY, WT, WR (year, theta and rho of the last observation) from WDS catalog
*
* INPUT:
* PISCO_cat: zeiss_doppie_new.cat (used to obtain the WDS names and orbit info.)
* WDS_cat: wdsweb_summ.txt (used to obtain the last observations)
*
**************************************************************************/
int astrom_add_WDS(FILE *fp_in, FILE *fp_out, char *PISCO_cat, char *WDS_cat)
{
double magV, B_V, paral, err_paral, magV_A, magV_B; 
double year, WdsLastYear, WdsLastTheta, WdsLastRho, mag_A, mag_B;
int has_an_orbit, iline, contains_object_name, contains_WDS_name, status;
int found, wds_has_orbit;
char b_in[NMAX], ads_name[40], discov_name[40], *pc;
char spectral_type[40], discov_name2[40], object_name[40], comp_name[40]; 
char ads_name2[40], WDS_name2[40], orbit[10], WDS_name3[40];
char NameInPiscoCatalog[40];

iline = 0;
while(!feof(fp_in))
{
/* Maximum length for a line will be 170 characters: */
  if(fgets(b_in,170,fp_in))
  {
  iline++;
  b_in[169] = '\0';
/* Remove ^M (Carriage Return) if present: */
  pc = b_in;
  while(*pc) {
  if(*pc == '\r') *pc = ' ';
  pc++;
  }
 
/* Check if input LateX line contains the name of the object:
*/
  astrom_check_if_object_name(b_in, &contains_object_name,
                              &contains_WDS_name);
/* Add the WDS name and discoverer's name if WDS name is not present:
*/
     if(contains_object_name && !contains_WDS_name) {
#ifdef DEBUG
        printf("astrom_add_WDS/DEBUG: line=%s", b_in);
#endif
        astrom_get_name_from_2nd_col(b_in, ads_name, discov_name, comp_name, 
                                     &year);
/* Look for object_name in file PISCO_catalog_name ("zeiss_doppie.cat"), 
* and determine values of: alpha, delat, coord_equinox,
*                          has_an_orbit, discov_name0, comp_name0.
*/
        if(*ads_name) sprintf(object_name, "%s%s", ads_name, comp_name); 
          else sprintf(object_name, "%s%s", discov_name, comp_name); 
#ifdef DEBUG
        printf("DDEBUG: ads_name=%s< discov=%s< object_name=%s< comp_name=%s< year=%f", 
               ads_name, discov_name, object_name, comp_name, year);
#endif
        if(*object_name) {

/* Removes AB if alone after ADS name: ADS 3456AB */ 
          pc = ads_name;
          while(*pc && !isdigit(*pc)) pc++;
          while(*pc &&  isdigit(*pc)) pc++;
          if(!strncmp(pc,"AB ",3) || !strcmp(pc,"AB")) *pc = '\0';

/* Will look for object name in PISCO catalog: */
          strcpy(NameInPiscoCatalog, object_name);
/* Removes AB if alone after ADS name: ADS 3456AB */ 
          if(*comp_name != '\0') {
            pc = NameInPiscoCatalog; 
            while(*pc && !isdigit(*pc)) pc++;
            while(*pc &&  isdigit(*pc)) pc++;
            if(!strncmp(pc,"AB ",3) || !strcmp(pc,"AB")) *pc = '\0';
          }
#ifdef DEBUG
printf("uuuu: name_in_pisco=%s<, object_name=%s\n",NameInPiscoCatalog,
       object_name);
#endif
           
          status = get_data_from_PISCO_catalog(PISCO_cat,
                             NameInPiscoCatalog, &magV, &B_V, &paral, 
                             &err_paral, &magV_A, &magV_B, spectral_type, 
                             discov_name2, ads_name2, WDS_name2, &has_an_orbit);
          if(status) {
             fprintf(stderr,"get_data_from_PISCO_catalog/Error \
object=%s not found in PISCO catalog \n (i.e., \"%s\") \n", 
                     object_name, PISCO_cat);
/*
             fprintf(stderr,"astrom_add_WDS/Fatal error: stop here\n"); 
             exit(-1);
*/
          } else {

/* JLP2010: possibility of retrieving ADS name from PISCO catalog: */
               if(ads_name[0] == '\0') {
                  strcpy(ads_name, ads_name2);
                  }

               if(has_an_orbit) strcpy(orbit, "orb");
                  else orbit[0] = '\0';
#ifdef DEBUG
           printf("(ads_name2=%s discov_name2=%s< WDS_name2=>%s< orbit=%d< stat=%d)\n", 
                   ads_name2, discov_name2, WDS_name2, has_an_orbit, status);
#endif
               get_data_from_WDS_catalog(WDS_cat, discov_name2, WDS_name3,
                                         &WdsLastYear, &WdsLastRho, 
                                         &WdsLastTheta, &mag_A, &mag_B, 
                                         spectral_type, &wds_has_orbit, &found);
               if(found) {
                  sprintf(b_in, "%s = %s%s & %s & %d & & & & & & & %s WY=%d WT=%d WR=%3.1f \\\\\n",
                         WDS_name2, discov_name2, comp_name, 
                         ads_name, (int)year, orbit, (int)WdsLastYear, 
                         (int)WdsLastTheta, WdsLastRho);
               } else {
                  sprintf(b_in, "%s = %s%s & %s & %d & & & & & & & %s \\\\\n",
                         WDS_name2, discov_name2, comp_name,
                         ads_name, (int)year, orbit);
               } 
            } /* EOF status == 0 */
          } /* EOF object_name != 0 */
      } /* EOF contains_object_name */

// Copy current line to file:
        fputs(b_in, fp_out);
  } /* EOF if fgets() */
} /* EOF while loop */

return(0);
}
/*************************************************************************
* astrom_calib_write_miniheader
*
*************************************************************************/
int astron_calib_write_miniheader(FILE *fp_out)
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
/*************************************************************************
* astrom_calib_copy
* To generate a calibrated file, without any other modifications
*
* INPUT: 
* calib_scale1 : scale values, scale_10, scale_20, scale_32
* calib_eyepiece1 : focal lengths, 10, 20, 32
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
* comments_wanted: 1 if comments are wanted, 0 otherwise
* input_with_header: 1 if input Latex file has a full header
*                    0 if no header at all
* 
*************************************************************************/
int astrom_calib_copy(FILE *fp_in, FILE *fp_out, double *calib_scale1, 
                      int *calib_eyepiece1, int n_eyepieces1, double theta0, 
                      double sign, int i_date, int i_eyepiece, int i_rho, 
                      int i_drho, int i_theta, int i_dtheta, 
                      int comments_wanted, int input_with_header)
{
char b_in[NMAX], b_data[NMAX];
int inside_array, line_is_opened, status, iline;
int line_with_object_name, nlines, nl_max; 
int contains_object_name, contains_WDS_name;
char *pc, *pc1;

/* Landscape: nl_max= 50 */
nl_max = 50;
// nlines: number of decoded lines
nlines = 0;
// iline: line number in file
iline = 0;

/* Assume we are inside the array if (input_with_header == 0): */
if(input_with_header == 0) {
  astron_calib_write_miniheader(fp_out);
  inside_array = 1;
} else {
  inside_array = 0;
}

line_is_opened = 0;

while(!feof(fp_in))
{
/* Maximum length for a line will be 170 characters: */
  if(fgets(b_in,170,fp_in))
  {
  iline++;
  line_with_object_name = 0;
  b_in[169] = '\0';
/* NEW/2009: I remove ^M (Carriage Return) if present: */
  pc = b_in;
  while(*pc) {
  if(*pc == '\r') *pc = ' ';
  pc++;
  }
 
/* Analyse this line: */
    if(!strncmp(b_in,"\\begin{tabular}",15)){
       inside_array = 1;
#ifdef DEBUG
printf(" OK: >%s< inside_array=%d\n", b_in, inside_array);
#endif
        }
    else if(!strncmp(b_in,"\\end{tabular}",13)){
       inside_array = 0;
#ifdef DEBUG
printf(" OK: >%s< inside_array=%d\n", b_in, inside_array);
#endif
       }
// JLP2009: to correct the possible case of " %%" instead of "%%":
    else if(inside_array && (b_in[0] != '%' && b_in[1] != '%')
            && strncmp(b_in,"\\hline",6)) {
       if(!line_is_opened) {
         strcpy(b_data, b_in);
         }
/* Fill the data array with the next line */
       else {
/* Look for the first zero (end of string marker) or CR ('\r') in data buffer */
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
#ifdef DEBUG
printf(" Data line: >%s<\n", b_data);
printf(" line_is_opened=%d\n", line_is_opened);
#endif
     if(!line_is_opened) { 
       status = astrom_calib_data_copy(b_data, b_in, calib_scale1, 
                                       calib_eyepiece1, n_eyepieces1,
                                       theta0, sign, i_date, 
                                       i_eyepiece, i_rho, i_drho, i_theta, 
                                       i_dtheta);
       if(status == -1) {
          printf("Fatal error in line %d status=%d\n", iline, status);
          exit(-1);
          }
       }
    } 

/* Write to Latex output file: */
    if(!line_is_opened) {
// JLP2009: to correct the possible case of " %%" instead of "%%":
      if(comments_wanted  || (b_in[0] != '%' && b_in[1] != '%')) {
        if(inside_array && (b_in[0] != '%' && b_in[1] != '%')) nlines++;
/* End table if page is full and current line contains an object name*/
        astrom_check_if_object_name(b_in, &contains_object_name,
                                    &contains_WDS_name);
        if(nlines > nl_max && !contains_object_name) {
          fprintf(fp_out,"\\hline \n");
          fprintf(fp_out,"\\end{tabular} \n");
          fprintf(fp_out,"\\end{table*} \n");
/* Before: \\newpage, but problems with "Too many unprocessed floats"
*/
          fprintf(fp_out,"\n\\clearpage \n");
          fprintf(fp_out,"\\begin{table*} \n");
          fprintf(fp_out,"\\tabcolsep=1mm \n");
          fprintf(fp_out,"\\caption{Relative astrometry of binary stars (cont.)} \n");
          fprintf(fp_out,"\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|} \n");
/* Contained in file:
          fprintf(fp_out,"\\hline \n");
*/
/*
          fprintf(fp_out,"& & & & & & & & & \\\\ \n");
*/
          nlines = 0;
        }
// Write current line content to file:
        fputs(b_in,fp_out);
      } /* EOF comments_wanted, b_in != % */ 
    }  /* EOF !line_is_opened */
  } /* EOF if fgets() */
} /* EOF while loop */

/* Epilog: */
     fprintf(fp_out,"\\hline \n");
     fprintf(fp_out,"\\end{tabular} \n");
     fprintf(fp_out,"\\end{table*} \n");
     fprintf(fp_out,"\\end{document} \n");

printf("End at line %d\n", iline);

return(0);
}
/*************************************************************************
* Publication mode
*
* INPUT:
* calib_scale1 : scale values, scale_10, scale_20, scale_32
* calib_eyepiece1 : focal lengths, 10, 20, 32
* i_filename: column nber of the filename used for this measurement
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
* comments_wanted: 1 if comments are wanted, 0 otherwise
* 
*************************************************************************/
int astrom_calib_publi(FILE *fp_in, FILE *fp_out,  double *calib_scale1,
                       int *calib_eyepiece1, int n_eyepieces1, double theta0, 
                       double sign, int i_filename, int i_date, int i_filter,
                       int i_eyepiece, int i_rho, int i_drho,
                       int i_theta, int i_dtheta, int i_notes,
                       int comments_wanted, char *filein, int gili_format)
{
OBJECT *obj;
int *index_obj, tabular_only;
int i, nobj = 0;

if((obj = (OBJECT *)malloc(NOBJ_MAX * sizeof(OBJECT))) == NULL) {
  printf("astrom_calib_publi/Fatal error allocating memory space for OBJECT: nobj_max=%d\n", 
          NOBJ_MAX);
  exit(-1);
  }
// Initialize the number of measurements to zero
for(i = 0; i < NOBJ_MAX; i++) (obj[i]).nmeas = 0;

if((index_obj = (int *)malloc((NOBJ_MAX) * sizeof(int))) == NULL) {
  printf("astrom_calib_publi/Fatal error allocating memory space for index_obj: nobj_max=%d\n", 
          NOBJ_MAX);
  exit(-1);
  }

astrom_read_measures(fp_in, comments_wanted, obj, &nobj, i_filename, i_date, 
                     i_filter, i_eyepiece, i_rho, i_drho, i_theta, i_dtheta, 
                     i_notes);

#ifdef DEBUG
printf("Returned by astrom_read_measures: nobj = %d\n", nobj);
#endif

astrom_calibrate_measures(obj, nobj, calib_scale1, calib_eyepiece1, 
                          n_eyepieces1, theta0, 
                          sign);

/* Sort the objects according to their Right Ascension: */
astrom_ra_sort_objects(obj, index_obj, nobj);

/* Compute mean values for the objects with rho < 0.3" 
*  and discard all measurements on recorded data if rho >= 0.3" (Merate-Paper II) */
astrom_mean_for_paper2(obj, nobj);

/* For big tables, should set this parameter to 1: */
tabular_only = 0;
if(gili_format == 1)
   astrom_write_publi_table_gili(fp_out, comments_wanted, obj, index_obj, nobj,
                         tabular_only);
else
   astrom_write_publi_table(fp_out, comments_wanted, obj, index_obj, nobj,
                         tabular_only);

astrom_compute_statistics(fp_out, obj, nobj, filein);

free(index_obj);
free(obj);
return(0);
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
    fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
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
    fprintf(fp_out,"\\begin{tabular*}{\\textwidth}{clrcccccrccl} \n");
    fprintf(fp_out,"\\hline \n");
    fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
    fprintf(fp_out,"WDS & Name & ADS & Epoch & Fil. & Eyep. & $\\rho$ & $\\sigma_\\rho$ & \\multicolumn{1}{c}{$\\theta$} & $\\sigma_\\theta$ & Orb. & Notes \\\\ \n");
    fprintf(fp_out,"& &       &       &        & (mm)& (\") & (\") & ($^\\circ$) & ($^\\circ$) & & \\\\ \n");
    fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
    fprintf(fp_out,"\\hline \n");
    fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
    nlines = 0;
  }

  jj = 0;
/************* Loop on measurements: **********************/
  for(j = 0; j < nm; j++) {
    me = &(obj[io]).measure[j];
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
       me_next = &(obj[io]).measure[jnext];
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
    if(status) strcpy(qflag,exclam);
    }

  *q_error = (me->dquadrant == 1) ? '?' : ' ';
#ifdef Q_IN_NOTES
    sprintf(q_notes, "q=%d%c", me->quadrant, *q_error);
#else
    *q_notes = '\0';
#endif

  if(jj == 0) {
     astrom_preformat_wds_name(obj[io].wds, wds_name);
/*** Format for first line of an object */
/* October 2008: 3 decimals for the epoch */
  if(me->rho != NO_DATA) 
    fprintf(fp_out,"%s & %s & %s & %.3f & %s & %d & %.3f & %.3f & %.1f%s & %.1f & %d & %s %s\\\\\n", 
         wds_name, obj[io].name, obj[io].ads, me->epoch, me->filter, 
         me->eyepiece, me->rho, me->drho, me->theta, qflag, me->dtheta, 
         obj[io].orbit, me->notes, q_notes);
  else
    fprintf(fp_out,"%s & %s & %s & %.3f & %s & %d & \\nodata & \\nodata & \\nodata & \\nodata & %d & %s \\\\\n", 
         wds_name, obj[io].name, obj[io].ads, me->epoch, me->filter, 
         me->eyepiece, obj[io].orbit, me->notes);
  } /* EOF j == 0 */
/*** Format for subsequent lines of an object */
  else {
  if(me->rho != NO_DATA) 
    fprintf(fp_out,"\\idem & \\idem & \\idem & %.3f & %s & %d & %.3f & %.3f & %.1f%s & %.1f & %d & %s %s\\\\\n", 
         me->epoch, me->filter, 
         me->eyepiece, me->rho, me->drho, me->theta, qflag, me->dtheta, 
         obj[io].orbit, me->notes, q_notes);
  else
    fprintf(fp_out,"\\idem & \\idem & \\idem & %.3f & %s & %d & \\nodata & \\nodata & \\nodata & \\nodata & %d & %s \\\\\n", 
         me->epoch, me->filter, me->eyepiece, obj[io].orbit, me->notes);
  } /* EOF j != 0 */
 nlines++;
 jj++;
 } /* EOF case not_flagged */
 } /* EOF loop on j */
} /* EOF loop on i */

 fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
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
* without the filter, the ADS name and the orbit
*****************************************************************************/
int astrom_write_publi_table_gili(FILE *fp_out, int comments_wanted, 
                                  OBJECT *obj, int *index_obj, int nobj, 
                                  int tabular_only)
{
MEASURE *me, *me_next;
int good_q, nm, io, nlines, jj, jnext, nl_max, status;
char qflag[20], asterisk[20], exclam[20], q_error[1], wds_name[40]; 
char q_notes[20], dmag_string[16];
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
//    fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
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
    nlines = 0;
  }

  jj = 0;
/************* Loop on measurements: **********************/
  for(j = 0; j < nm; j++) {
    me = &(obj[io]).measure[j];
/* BOF case not_flagged */
    if(!me->flagged_out) {
      if(me->rho <= 0.1) {
// Unresolved case:
         fprintf(fp_out,"%s & %s & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata &  & Unres. \\\\\n", 
             wds_name, obj[io].name, me->epoch, me->eyepiece); //  me->notes);
      } else {
/* For DEBUG:
if(i < 4) {
  printf("i= %d io=%d nm=%d j=%d flag=%d\n", i, io, nm, j, me->flagged_out);
//  printf(" WDS=%s rho=%.2f drho=%.2f theta=%.2f dtheta=%.2f eyepiece=%d\n", 
//           obj[io].wds, me->rho, me->drho, me->theta, me->dtheta, me->eyepiece);
  printf(" WDS=%s rho=%.2f drho=%.2f theta=%.2f dtheta=%.2f\n", 
           obj[io].wds, me->rho, me->drho, me->theta, me->dtheta, me->eyepiece);
  }
*/

/* JLP2007: I check if among the next measurements there exists
* another one with the same epoch and filter
* than the current measurement: if it is the case, compute
* the mean of the two: */
    for(jnext = j+1; jnext < nm; jnext++) {
       me_next = &(obj[io]).measure[jnext];
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
    if(status) strcpy(qflag,exclam);
    }

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
/*    fprintf(fp_out,"%s & %s & %s & %.3f & %s & %d & %.3f & %.3f & %.1f%s & %.1f & %d & %s %s\\\\\n", 
         wds_name, obj[io].name, obj[io].ads, me->epoch, me->filter, 
         me->eyepiece, me->rho, me->drho, me->theta, qflag, me->dtheta, 
         obj[io].orbit, me->notes, q_notes);
*/
    fprintf(fp_out,"%s & %s & %.3f & %d & %.3f & %.3f & %.1f%s & %.1f & %s &\\\\\n", 
         wds_name, obj[io].name, me->epoch, me->eyepiece, me->rho, me->drho, 
         me->theta, qflag, me->dtheta, dmag_string); //, me->notes, q_notes);
  else
/*
    fprintf(fp_out,"%s & %s & %s & %.3f & %s & %s & %d & \\nodata & \\nodata & \\nodata & \\nodata & %d & \\\\\n", 
         wds_name, obj[io].name, obj[io].ads, me->epoch, me->filter, 
         me->eyepiece, obj[io].orbit); //, me->notes);
*/
    fprintf(fp_out,"%s & %s & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata & \\nodata &  \\\\\n", 
         wds_name, obj[io].name, me->epoch, me->eyepiece); //, me->notes);
  } /* EOF j == 0 */
/*** Format for subsequent lines of an object */
  else {
  if(me->rho != NO_DATA) 
/*
    fprintf(fp_out,"\\idem & \\idem & \\idem & %.3f & %s & %d & %.3f & %.3f & %.1f%s & %.1f & %d & \\\\\n", 
         me->epoch, me->filter, 
         me->eyepiece, me->rho, me->drho, me->theta, qflag, me->dtheta, 
         obj[io].orbit); //, me->notes, q_notes);
*/
    fprintf(fp_out,"\\idem & \\idem & %.3f & %d & %.3f & %.3f & %.1f%s & %.1f & %s & \\\\\n", 
         me->epoch,  me->eyepiece, me->rho, me->drho, me->theta, 
         qflag, me->dtheta, dmag_string); //, me->notes, q_notes);
  else
/*
    fprintf(fp_out,"\\idem & \\idem & \\idem & %.3f & %s & %d & \\nodata & \\nodata & \\nodata & \\nodata & %d & %s \\\\\n", 
         me->epoch, me->filter, me->eyepiece, obj[io].orbit, me->notes);
*/
    fprintf(fp_out,"\\idem & \\idem & %.3f & %d & \\nodata & \\nodata & \\nodata & \\nodata & \nodata & \\\\\n", 
         me->epoch, me->eyepiece); //, me->notes);
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
*
* OUTPUT:
* obj: OBJECT structure 
* nobj: total number of objects entered into *obj
*****************************************************************************/
int astrom_read_measures(FILE *fp_in, int comments_wanted, OBJECT *obj, 
                         int *nobj, int i_filename, int i_date, 
                         int i_filter, int i_eyepiece, int i_rho, 
                         int i_drho, int i_theta, int i_dtheta, int i_notes)
{
char b_in[NMAX], b_data[NMAX];
char wds_name[40], discov_name[40], ads_name[40]; 
int inside_array, line_is_opened, status, orbit, line_to_reject, i_obj;
int input_with_header;
char *pc, *pc1;

printf("ZZZ read_measures: nobj=%d\n", *nobj);

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
ads_name[0] = '\0';
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
#ifdef DEBUG
printf(" Data line: >%s<\n", b_data);
#endif
#ifdef DEBUG
printf(" Checking if it is a new object: #nobj=%d\n", *nobj);
#endif
/* Try to add a new object: */
       status = astrom_add_new_object(b_data, obj, nobj, i_notes);
       if(status == 0) {
/* Index of current object in OBJECT structure array: */
          i_obj = *nobj - 1;
          (obj[i_obj]).nmeas = 0;
          }
printf(" OK: status_obj=%d index of object #i_obj=%d nm=%d (nobj=%d)\n", 
         status, i_obj, (obj[i_obj]).nmeas, nobj);

/* Try to add a new measure: */
       if((status != 0) && (*nobj > 1)) {
         status = astrom_check_measure(b_data, i_eyepiece, i_rho, i_drho, 
                                       i_theta, i_dtheta);
#ifdef DEBUG
printf(" Adding new measurement for object #i_obj=%d nm=%d (nobj=%d)\n", 
         i_obj, (obj[i_obj]).nmeas, nobj);
#endif
          if(status == 0) {
           astrom_add_new_measure(b_data, obj, i_obj, i_filename, i_date, 
                           i_filter, i_eyepiece, i_rho, i_drho, i_theta, 
                           i_dtheta, i_notes, comments_wanted);
           } 
         }

       } /* EOF !line_is_opened */
    } // EOF line tabular 
  } /* EOF if fgets() */
} /* EOF while loop */
return(0);
}
/**************************************************************************
* Check if current line is that of a new object.
* If so, load corresponding data into new OBJECT structure
*
* Check if current line is compatible with the syntax of a new object
* Example:
* 16564+6502 = STF 2118 AB & ADS 10279 & 2004. & & & & & & & orb \\
* Or:
*  & ADS 10279 & 2004. & & & & & & & \\
* Or:
  COU 2118 &  & 2004. & & & & & & & orb \\
*
*
* i_notes: column nber with the notes
**************************************************************************/
int astrom_add_new_object(char *b_data, OBJECT *obj, int *nobj, int i_notes)
{ 
OBJECT *ob;
double WR, WT, WY;
char wds_name[40], discov_name[40], ads_name[40], *pc, buffer[40];
int status, orbit, ih, im;

if(*nobj > NOBJ_MAX - 1) {
  fprintf(stderr, 
          "astrom_add_new_object/Fatal error: maximum number of objects=%d is overtaken! \n", 
           NOBJ_MAX);
  exit(-1);
  }

ob = &obj[*nobj];

status = astrom_read_object_data(b_data, wds_name, discov_name, ads_name, 
                                 &orbit, &WR, &WT, &WY, i_notes);

printf("ZZZ astrom_add_new_object/status=%d nobj = %d nm=%d\n", 
       status, *nobj, ob->nmeas);

if(status == 0) {
  strcpy(ob->wds, wds_name);
  strcpy(ob->ads, ads_name);
  strcpy(ob->name, discov_name);
  ob->nmeas = 0;
  ob->orbit = orbit;


/* Decode Right Ascension from WDS name: */
  ob->ra = 0.;
  strcpy(buffer, wds_name);
  pc = buffer;
  while(*pc && (*pc != '+') && (*pc != '-')) pc++;
  if((*pc == '+') || (*pc == '-')) {
    *pc = '\0';
    sscanf(buffer,"%2d%3d", &ih, &im);
/* Divide "im" by 600, since "im" is measured in tenths of minutes: */
    ob->ra = (double)ih + ((double)im)/600.;
    }
/* Decode Declination from WDS name: */
  ob->dec = 0.;
  pc = wds_name;
  while(*pc && (*pc != '+') && (*pc != '-')) pc++;
  if((*pc == '+') || (*pc == '-')) {
    strcpy(buffer,pc);
    sscanf(buffer,"%d", &ih);
    ob->dec = ih;
    }
  
/* WDS data, if available: */
  ob->WR = WR;
  ob->WT = WT;
  ob->WY = WY;

// Increase the number of objects:
  (*nobj)++;

#ifdef DEBUG
  printf(" New object added: nobj = %d\n", *nobj);
  printf(" name= %s %s %s orbit=%d", ob->wds, ob->name, ob->ads, ob->orbit);
  printf(" ra= %f dec=%d WR=%f WT=%f WY=%f nm=%d\n", ob->ra, ob->dec, 
           ob->WR, ob->WT, ob->WY, ob->nmeas);
#endif
  } else {
#ifdef DEBUG
  fprintf(stderr," Failure adding new object (nobj = %d) status =%d\n", 
          *nobj, status);
#endif
  }

return(status);
}
/*************************************************************************
* Compute mean values for the objects with rho < 0.3" 
*  and discard all measurements on recorded data if rho >= 0.3" (Merate-Paper II)
**************************************************************************/
int astrom_mean_for_paper2(OBJECT *obj, int nobj)
{
MEASURE *me, *me_prev;
int rec_dir, rec_dir_prev, nm, no_data_prev, no_data;
register int i, j;


for(i = 0; i < nobj; i++) {
  nm = (obj[i]).nmeas;
  me_prev = NULL;
  for(j = 0; j < nm; j++) {
    me = &(obj[i]).measure[j];
/* rec_dir = 1  if recorded file
*          = -1 if direct file
*          = 0 otherwise
*/
    if(astrom_is_record_file(me->filename) == 1)
      rec_dir = 1;
    else if(astrom_is_direct_file(me->filename) == 1)
      rec_dir = -1;
    else
      rec_dir = 0;
    no_data = ((me->rho == NO_DATA) || (me->theta == NO_DATA)) ? 1 : 0; 
/* Store current parameters (they become "previous" parameters in next 
* iteration if same object) */
    if(me_prev == NULL) {
      rec_dir_prev = rec_dir;
      no_data_prev = no_data;
      me_prev = me;
      } else {
/* Process measurements when two measurements refer to same epoch,
* with direct/recorded files: */ 
      if((me->epoch == me_prev->epoch) 
         && (me->filter[0] == me_prev->filter[0])
         && (me->eyepiece == me_prev->eyepiece)
         && (rec_dir * rec_dir_prev) == -1) {
/* If two measurements have been made: */
        if(!no_data  && !no_data_prev) {
/* Handle quadrant: 
* quadrant: 0 if Q was not present
*           -1 if Q=?
*           1, 2, 3 or 4 if Q=1, Q=2, Q=3 or Q=4
*/
          if(me->quadrant == -1 && me_prev->quadrant == 0) me_prev->quadrant = -1; 
          if(me->quadrant == 0 && me_prev->quadrant == -1) me->quadrant = -1; 
          if(me->quadrant > 0 && me_prev->quadrant <= 0) {
              me_prev->quadrant = me->quadrant; 
              me_prev->dquadrant = me->dquadrant; 
              }
          if(me_prev->quadrant > 0 && me->quadrant <= 0) {
              me->quadrant = me_prev->quadrant; 
              me->dquadrant = me_prev->dquadrant; 
              }
/* Handle rho and theta: */
          if(me->rho > 0.3) {
/* if rec and dir_prev */
            if(rec_dir == 1 && rec_dir_prev == -1) me->flagged_out = 1;
/* if dir and rec_prev */
            if(rec_dir == -1 && rec_dir_prev == 1) me_prev->flagged_out = 1;
            }
          else astrom_compute_mean_of_two_measures(obj, i, j-1, j);
        } /* EOF !no_data && !no_data_prev */
/* Special handling when NO_DATA were present for both measurements: */
/* Removes the second measurement (NB: the 2nd is the same as the first): */
        else if(no_data  && no_data_prev) me->flagged_out = 1;
/* Case when only one mesurement over two files: */
        else {
         if(no_data) me->flagged_out = 1;
         else if(no_data_prev) me_prev->flagged_out = 1;
         }
/* Second file was linked to first, I neutralize it for further inspection */
      me_prev = NULL;
       } /* EOF epoch == epoch_prev */
/* Second file was not linked to first, it may be linked to next: */
  else {
       rec_dir_prev = rec_dir;
       no_data_prev = no_data;
       me_prev = me;
       }
     } /* EOF me_prev != NULL */

/* Go to next measurement */
  } /* EOF loop on j (measurements) */
/* Go to next object */
} /* EOF loop on i  (objects) */

return(0);
}
/*************************************************************************
* Compute mean values for the objects with rho < 0.3" 
*  and discard all measurements on recorded data if rho >= 0.3" (Merate-Paper II)
**************************************************************************/
int astrom_mean_for_full_table(OBJECT *obj, int nobj)
{
MEASURE *me, *me_now;
int nj0, j0[256], j_used[256], quadrant_is_found, quadrant_value; 
int nm, ndone;
register int i, j;


for(i = 0; i < nobj; i++) {
  nm = (obj[i]).nmeas;
/* Flags set to one if the corresponding measurement was used 
* for computing the mean in one of the filters */
  for(j = 0; j < nm; j++) j_used[j] = 0;
  ndone = 0;

/* First root is the first point: */
  for(j = 0; j < nm; j++) j0[j] = 0;
  me_now = &(obj[i]).measure[0];
  j0[0] = 0;
  nj0 = 1;
  j_used[0] = 1;

  quadrant_value = me_now->quadrant;
  if(me_now->quadrant > 0) {
       quadrant_is_found = 1;
   } else {
       quadrant_is_found = 0;
   }

/* Scan all the possible roots */
  while(ndone < nm) {
/* Scan the available measurements looking for the observations made at the same epoch
* with the same filter
*/
  for(j = 1; j < nm; j++) {
    if(!j_used[j]) {
       me = &(obj[i]).measure[j];
/* Record all the measurements refering to same epoch,
*/
      if((me->epoch == me_now->epoch) 
         && (me->filter[0] == me_now->filter[0])
         && (me->eyepiece == me_now->eyepiece)
         && !me->flagged_out ) {
         j_used[j] = 1;
         j0[nj0] = j;
         nj0++; 
/* Handle quadrant: 
* quadrant: 0 if Q was not present
*           -1 if Q=?
*           1, 2, 3 or 4 if Q=1, Q=2, Q=3 or Q=4
*/
/* Store quadrant if present: */
          if(me->quadrant != 0) {
              if(me->quadrant != -1) {
                 me_now->quadrant = me->quadrant; 
                 quadrant_is_found = 1;
                 quadrant_value = me_now->quadrant;
              } else if(me_now->quadrant == 0) me_now->quadrant = -1; 
            }
      } /* EOF case of same epoch, same filter, ... */
     } /* EOF !j_done[j] */
/* Go to next measurement */
  } /* EOF loop on j (measurements) */
  if(nj0 > 1) compute_mean_for_filter(obj, i, j0, nj0, quadrant_value);
  ndone++;

/* Next root is the next free point: */
  for(j = 0; j < nm; j++) j0[j] = 0;
  for(j = 0; j < nm; j++) {
    me = &(obj[i]).measure[j];
    if(!j_used[j] && !me->flagged_out) break;
    }
/* Exit from loop if all points have been used: */
  if(j >= nm) {
   break;
/* Otherwise select the new root: */
  } else {
  me_now = &(obj[i]).measure[j];
  j_used[j] = 1;
  j0[0] = j;
  nj0 = 1;
  }

  } /* EOF while(1) */
/* Go to next object */
} /* EOF loop on i  (objects) */

return(0);
}
/********************************************************************** 
* Compute the mean of all measurements concerning this filter: 
*
* INPUT:
*  obj: list of all measurements
*  i: index of current object
*  j0[]: array of indices of measurements made with the current filter
*  nj0: number of measurements made with the current filter
*  quadrant value: value determined in this filter
**********************************************************************/
int compute_mean_for_filter(OBJECT *obj, int i, int *j0, int nj0, 
                            int quadrant_value)
{
MEASURE *me;
double rho, theta, drho, dtheta;
int j, no_data, n_measures;

n_measures = 0;
rho = 0.; theta = 0.;
drho = 0.; dtheta = 0.;
/* Compute the sum of all measurements concerning this filter: */
  for(j = 0; j < nj0; j++) {
    me = &(obj[i]).measure[j0[j]];
    no_data = ((me->rho == NO_DATA) || (me->theta == NO_DATA)) ? 1 : 0; 
    if(!no_data) {
       rho += me->rho;
       drho += me->drho;
       theta += me->theta;
       dtheta += me->dtheta;
       n_measures++;
       }
   } /* EOF loop on j */
/* Compute the mean if possible: */
if(n_measures > 0) {
  rho /= (double)n_measures;
  drho /= (double)n_measures;
  theta /= (double)n_measures;
  dtheta /= (double)n_measures;
  } else {
  rho = NO_DATA;
  drho = NO_DATA;
  theta = NO_DATA;
  dtheta = NO_DATA;
  }

/* Store result on the first element and disable the others: */
  me = &(obj[i]).measure[j0[0]];
  me->rho = rho;
  me->drho = drho;
  me->theta = theta;
  me->dtheta = dtheta;
  me->quadrant = quadrant_value;
  for(j = 1; j < nj0; j++) {
    me = &(obj[i]).measure[j0[j]];
    me->flagged_out = 1;
    }

return(0);
}
/************************************************************************
* Compute the mean of two measures
* and load the result in the first measure
*
* io: object index in obj
* jm1, jm2: measurement indices of the two measures in obj[io].measure
************************************************************************/
int astrom_compute_mean_of_two_measures(OBJECT *obj, int io, int jm1, int jm2)
{
MEASURE *me1, *me2;
double w1, w2, sum;
int no_data1, no_data2;

me1 = &(obj[io]).measure[jm1];
me2 = &(obj[io]).measure[jm2];

no_data1 = ((me1->rho == NO_DATA) || (me1->theta == NO_DATA)) ? 1 : 0; 
no_data2 = ((me2->rho == NO_DATA) || (me2->theta == NO_DATA)) ? 1 : 0; 

/* Neutralizes the measurement that lead to a non-detection
* and only keeps the good one (or the first one if two non-detections): */
if(no_data2) {
me2->flagged_out = 1;
} else if(no_data1) {
me1->flagged_out = 1;
} else {
w1 = (me2->drho / (me1->drho + me2->drho))
     + (me2->dtheta / (me1->dtheta + me2->dtheta));
w2 = (me1->drho / (me1->drho + me2->drho))
     + (me1->dtheta / (me1->dtheta + me2->dtheta));
sum = w1 + w2;
  if(sum == 0) {
     printf("astrom_compute_mean_of_two_measures/Fatal: rms error is null for io=%d jm1=%d jm2=%d\n", io, jm1, jm2);
     exit(-1);
     }
w1 /= sum; 
w2 /= sum; 

#ifdef DEBUG
printf(" astrom_compute_mean_of_two_measures/WDS=%s ADS=%s NAME=%s io=%d jm1=%d jm2=%d dt1=%f dt2=%f dr1=%f dr2=%f w1=%f w2 =%f \n", 
       (obj[io]).wds, (obj[io]).ads, (obj[io]).name, io, jm1, jm2, me1->dtheta, 
       me2->dtheta, me1->drho, me2->drho, w1, w2);
#endif

/* Load the mean onto the first measure: */
me1->rho = me1->rho * w1 + me2->rho * w2;
me1->theta = me1->theta * w1 + me2->theta * w2;
me1->drho = sqrt((SQUARE(me1->drho) * w1 + SQUARE(me2->drho) * w2));
me1->dtheta = sqrt((SQUARE(me1->dtheta) * w1 + SQUARE(me2->dtheta) * w2));

/* Flag out the second measure to prevent further use: */
me2->flagged_out = 1;

} /* EOF !no_data1 && !no_data2 */

return(0);
}
/*****************************************************************************
* Example:
* & ADS 1630 bc & 2010 & & & & & & & \\
*
* OUTPUT:
*
* for this example:
* ads_name: "ADS 1630"
* discov_name: "" 
* comp_name: "bc"
* year: "2010"
*
* Other example:
* & MLR 377 Aa-B & 2001 & & & & & & & \\
* ads_name : "" 
* discov_name: "MLR 377"
* comp_name: "Aa-B"
* year: "2001"
*****************************************************************************/
int astrom_get_name_from_2nd_col(char *b_in, char *ads_name, char *discov_name,
                                 char *comp_name, double *year)
{
char buff[80], *pc;
int icol, istat, iverbose = 0;

ads_name[0] = '\0';
discov_name[0] = '\0';
*year = 0.;

/* Reads year in 3rd column: */
icol = 3;
istat = latex_read_fvalue(b_in, year, icol, iverbose);
if(istat) {
  fprintf(stderr, "astrom_get_name_from_2nd_col/Fatal error reading year in line: %s\n", b_in);
  exit(-1);
  }

/* Reads ADS name in 2nd column: */
icol = 2;
istat = latex_read_svalue(b_in, buff, icol);
if(!istat){
 pc = buff;
 while(*pc && strncmp(pc,"ADS",3)) pc++;
 if(!strncmp(pc,"ADS",3)){
    pc = buff;
/* Skip the leading blanks: */
    while(*pc && *pc == ' ') pc++;
    strcpy(ads_name, pc);
/* Extract companion name if present: */
    pc = ads_name;
    while(*pc && !isdigit(*pc)) pc++;
    while(*pc && isdigit(*pc)) pc++;
    strcpy(comp_name, pc);
    *pc = '\0';
/* Case when ADS is not present: */
    } else {
/* Example: MLR 377 Aa-B
*/
    pc = buff;
/* Skip the leading blanks: */
    while(*pc && *pc == ' ') pc++;
    strcpy(discov_name, pc);
/* Extract companion name if present: */
    pc = discov_name;
    while(*pc && !isdigit(*pc)) pc++;
    while(*pc && isdigit(*pc)) pc++;
    strcpy(comp_name, pc);
    *pc = '\0';
    }
 }
/* Removes leading and trailing blanks:
*/
  jlp_trim_string(comp_name, 40);
  if(!strcmp(comp_name,"ab")) strcpy(comp_name,"AB");
  if(!strcmp(comp_name,"ac")) strcpy(comp_name,"AC");
  if(!strcmp(comp_name,"bc")) strcpy(comp_name,"BC");

/* Removes blanks in discover name: */
  jlp_compact_string(discov_name, 40);

return(0);
}
/*************************************************************************
*
* Example of input:
& 270109\_ads1933\_sfd\_8\_a & 27/01/2009 & sf & 20 & 25.82 & 0.62 & -54.12 & 1.2 & Q=3? \\
*
*************************************************************************/
int astrom_check_if_measurement(char *b_in, char *fits_filename, char *date1,
                                int *eyepiece1, int *is_measurement)
{
char buff[80], *pc;
double ww;
int k, dd, mm, yy, ival, istat, icol, iverbose = 0;

*is_measurement = 0;
*fits_filename = '\0';
*date1 = '\0';
*eyepiece1 = 0;

icol = 2;
istat = latex_read_svalue(b_in, buff, icol);
if(istat) return(-1);
/* Cleanup fits file
* From "270909\_ads12334\_sfd\_8\_a" to "270909_ads12334_sfd_8_a"
*/
jlp_compact_string(buff, 80);
pc = buff;
k = 0;
while(*pc) {
 if(*pc != '\\') fits_filename[k++] = *pc;
 pc++;
 }
fits_filename[k] = '\0';

/* Do not exit here in case of error, since not always available ? */
/* Input format: 27/09/2009 */
icol = 3;
istat = latex_read_svalue(b_in, buff, icol);
jlp_compact_string(buff, 80);
buff[19] = '\0';
ival = sscanf(buff, "%02d/%02d/%4d", &dd, &mm, &yy);
/* Check if format is good: */
if(ival == 3) strcpy(date1, buff); 

/* Do not exit here in case of error, since not always available ? */
icol = 5;
istat = latex_read_ivalue(b_in, eyepiece1, icol);

/* Check if rho measurement is present (or "\nodata") */ 
icol = 6;
istat = latex_read_svalue(b_in, buff, icol);
jlp_compact_string(buff, 80);
if(!strncmp(buff, "\\nodata", 7)) {
  *is_measurement = 1;
 } else {
  istat = latex_read_fvalue(b_in, &ww, icol, iverbose);
  if(!istat) {
     *is_measurement = 1;
     }
  }

return(istat);
}
/*************************************************************************
* merate_get_full_directory
*
* INPUT:
* fits_directory (example: "/home/data/pisco_merate/2004-2008/" )
* date (example: "20/01/2009")
*
* OUTPUT:
* full_directory (example: "/home/data/pisco_merate/2004-2008/gen2009/" )
*************************************************************************/
int merate_get_full_directory(char *fits_directory, char *date, 
                              char *full_directory)
{
int nval, dd, mm, yy;
char month[12];

strcpy(full_directory, fits_directory);

/* Decode the date: */
nval = sscanf(date, "%2d/%2d/%4d", &dd, &mm, &yy);
if(nval != 3) return(-1);

switch(mm)
 {
 case 1:
   strcpy(month, "gen");
   break;
 case 2:
   strcpy(month, "feb");
   break;
 case 3:
   strcpy(month, "mar");
   break;
 case 4:
   strcpy(month, "apr");
   break;
 case 5:
   strcpy(month, "mag");
   break;
 case 6:
   strcpy(month, "giu");
   break;
 case 7:
   strcpy(month, "lug");
   break;
 case 8:
   strcpy(month, "ago");
   break;
 case 9:
   strcpy(month, "set");
   break;
 case 10:
   strcpy(month, "ott");
   break;
 case 11:
   strcpy(month, "nov");
   break;
 case 12:
   strcpy(month, "dic");
   break;
 default: 
   fprintf(stderr,"merate_get_full_directory/Fatal error reading month ind date: >%s< \n", date);
   exit(-1); 
   break;
 }

sprintf(full_directory, "%s%s%4d/", fits_directory, month, yy);

printf("DDEBUG/merate_get_full_directory:  full directory= >%s< \n", full_directory);

return(0);
} 

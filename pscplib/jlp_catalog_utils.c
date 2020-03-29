/************************************************************************
* "jlp_catalog_utils.c"
* Set of routines used 
* for CALIBrated table (from latex_calib.c) and Latex tables
*
* JLP 
* Version 05/05/2010
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>    // exit()
#include <string.h>    // strcpy()
#include "jlp_numeric.h" // JLP_QSORT, MINI, MAXI, PI, etc
#include "jlp_string.h" // jlp_trim_string (in jlplib/jlp_fits 

/* The prototypes of routines included here
* are defined in "jlp_catalog_utils.h":
*/
#include "jlp_catalog_utils.h"

#ifndef MAXI
#define MAXI(a,b) ((a) < (b)) ? (b) : (a)
#endif

/*
#define DEBUG 
#define DEBUG_1 
*/
static int read_series_of_measures(char *in_line, FILE *fp_calib_table, 
                                   char *object_name, 
                                   char *comp_name, char *name1,
                                   int comp_is_AB, double *epoch_o, 
                                   double *rho_o, double *theta_o,
                                   double *err_rho_o, double *err_theta_o, 
                                   int *nmeas);
static int read_series_of_measures_gili(char *in_line, FILE *fp_calib_table, 
                                   char *object_name, 
                                   char *comp_name, char *name1,
                                   int comp_is_AB, double *epoch_o, 
                                   double *rho_o, double *theta_o,
                                   double *err_rho_o, double *err_theta_o, 
                                   int *nmeas);

/*********************************************************************
* Check if str1 is contained in in_line1
*********************************************************************/
int is_in_line(char *in_line1, char *str1)
{
int found;
char *pc;
pc = in_line1;
found = 0;
while(*pc && !found) {
  if(!strncmp(str1, pc, strlen(str1))) found = 1;
  pc++;
}
return(found);
}
/***********************************************************************
* Retrieve all the measurements obtained at different epochs 
* for a given object in the calibrated Latex table
*
* INPUT:
* calib_fname: file name of the calibrated LateX table
* object_name: name of object (either ADS 123 or COU 432)
* comp_name: name of the companion (i.e. AB, Aa, CD, etc) 
*
* OUTPUT:
* epoch_o, rho_o (arcsec), theta_o (deg): observed measurements 
* err_rho_o (arcsec), err_theta_o (deg): errors of the observed measurements
* nmeas: number of measurements
***********************************************************************/
int get_measures_from_CALIB_table(char *calib_fname, char *object_name,
                                  char *comp_name, double *epoch_o, 
                                  double *rho_o, double *theta_o,
                                  double *err_rho_o, double *err_theta_o, 
                                  int *nmeas)
{
char in_line[300], name1[40], name2[40], comp_name1[40], comp_name2[40]; 
int icol, iline, object_is_ADS, object2_is_ADS, status;
int comp1_is_AB, comp2_is_AB;
char comp_really_compacted[40], comp2_really_compacted[40];
FILE *fp_calib_table;

    strcpy(comp_name1, comp_name);
    jlp_trim_string(comp_name1, 40);

/* Test the structure of the object name (either ADS 234 or COU 345, f.i.)
* Look for ADS name in the 3rd column*/
if(!strncmp(object_name, "ADS", 3)) {
    strcpy(name1, &object_name[3]);
    jlp_trim_string(name1, 40);
    object_is_ADS = 1;
    icol = 3;
/* Else look for the discoverer name in the 2nd column*/
    } else {
    strcpy(name1, object_name);
    jlp_compact_string(name1, 40);
    icol = 2;
    object_is_ADS = 0;
    }
/* Check if AB companion (default) or something else */
   comp1_is_AB = 0;
   if((comp_name1[0] == '\0') || !strcmp(comp_name1,"AB")
       || !strncmp(comp_name1,"Aa-B",4)) comp1_is_AB = 1;

/* Open LaTeX table: */
if((fp_calib_table = fopen(calib_fname, "r")) == NULL) {
  fprintf(stderr, "get_measures_from_CALIB_table/Error opening %s\n", 
          calib_fname);
  return(-1);
 }

/* Scan all the LaTeX file looking for the object name */ 
/* Example:
11000$-$0328  &  STF1500  &  8007  & 2010.387 & R  & 20 & 1.370 & 0.015 & 299.9  & 0.3 & 0 &   \\
*/
iline = 0;
*nmeas = 0;
while(!feof(fp_calib_table)){
  if(fgets(in_line,300,fp_calib_table)) {
  iline++;
/* Process only the meaningful lines (skipping the header...)*/
  if(isdigit(in_line[0])) {
    read_full_name_from_CALIB_line(in_line, name2, comp_name2, 
                                   &object2_is_ADS, &comp2_is_AB);
    jlp_trim_string(name2, 40);
    if(!strcmp(name1, name2)) {
     jlp_really_compact_companion(comp_name2, comp2_really_compacted, 40);
     jlp_really_compact_companion(comp_name, comp_really_compacted, 40);

/* If component are the same, add the corresponding measurements */
     if(!strcmp(comp_name, comp_name2) 
            || (comp1_is_AB && comp2_is_AB)
            || (*comp_really_compacted != '\0' && 
               !strcmp(comp_really_compacted, comp2_really_compacted))){
/* Read a series of measures including all the lines starting with \idem
* that follows the object
*/
             read_series_of_measures(in_line, fp_calib_table, 
                                     object_name, comp_name, 
                                     name1, comp1_is_AB, epoch_o, rho_o, 
                                     theta_o, err_rho_o, err_theta_o, nmeas);
            }
    } /* EOF case name1 == name2 */
  } /* EOF isdigit(in_line[0]) */
  } /* EOF fgets() */
} /* EOF while */

if(*nmeas == 0) {
fprintf(stderr, "get_measures_from_CALIB_table/Object >%s< comp=>%s< not found in calibrated Latex table!\n", object_name, comp_name);
fprintf(stderr, "(name1=%s icol=%d object_is_ADS=%d)\n", name1, icol, 
       object_is_ADS);
status = -1;
} else {
status = 0;
}

fclose(fp_calib_table);
return(status);
}
/***********************************************************************
* Retrieve all the measurements obtained at different epochs 
* for a given object in the calibrated Latex table
*
* INPUT:
* calib_fname: file name of the calibrated LateX table
* object_name: name of object (e.g., COU 432)
* comp_name: name of the companion (i.e. AB, Aa, CD, etc) 
*
* OUTPUT:
* epoch_o, rho_o (arcsec), theta_o (deg): observed measurements 
* err_rho_o (arcsec), err_theta_o (deg): errors of the observed measurements
* nmeas: number of measurements
***********************************************************************/
int get_measures_from_CALIB_table_gili(char *calib_fname, char *object_name,
                                  char *comp_name, double *epoch_o, 
                                  double *rho_o, double *theta_o,
                                  double *err_rho_o, double *err_theta_o, 
                                  int *nmeas)
{
char in_line[300], name1[40], name2[40], comp_name1[40], comp_name2[40]; 
int icol, iline, status;
int comp1_is_AB, comp2_is_AB;
char comp_really_compacted[40], comp2_really_compacted[40];
FILE *fp_calib_table;

    strcpy(comp_name1, comp_name);
    jlp_trim_string(comp_name1, 40);

/* Test the structure of the object name (COU 345, f.i.)
/* Look for the discoverer name in the 2nd column*/
    strcpy(name1, object_name);
    jlp_compact_string(name1, 40);
    icol = 2;
/* Check if AB companion (default) or something else */
   comp1_is_AB = 0;
   if((comp_name1[0] == '\0') || !strcmp(comp_name1,"AB")
       || !strncmp(comp_name1,"Aa-B",4)) comp1_is_AB = 1;
// JLP2020: reduce the name if comp_is_AB is true:
   if(comp1_is_AB == 1) remove_AB_from_object_name(name1);

/* Open LaTeX table: */
if((fp_calib_table = fopen(calib_fname, "r")) == NULL) {
  fprintf(stderr, "get_measures_from_CALIB_table_gili/Error opening %s\n", 
          calib_fname);
  return(-1);
 }

/* Scan all the LaTeX file looking for the object name */ 
/* Example:
11000$-$0328  &  STF1500  & 2010.387 & 2 & 1.370 & 0.015 & 299.9  & 0.3 & 1.20 &   \\
*/
iline = 0;
*nmeas = 0;
while(!feof(fp_calib_table)){
  if(fgets(in_line,300,fp_calib_table)) {
  iline++;
/* Process only the meaningful lines (skipping the header...)*/
  if(isdigit(in_line[0])) {
    read_full_name_from_CALIB_line_gili(in_line, name2, comp_name2, 
                                        &comp2_is_AB);
// JLP2020: reduce the name if comp_is_AB is true:
    if(comp2_is_AB == 1) remove_AB_from_object_name(name2);

    jlp_trim_string(name2, 40);
    if(!strcmp(name1, name2)) {
     jlp_really_compact_companion(comp_name2, comp2_really_compacted, 40);
     jlp_really_compact_companion(comp_name, comp_really_compacted, 40);

/* If component are the same, add the corresponding measurements */
     if(!strcmp(comp_name, comp_name2) 
            || (comp1_is_AB && comp2_is_AB)
            || (*comp_really_compacted != '\0' && 
               !strcmp(comp_really_compacted, comp2_really_compacted))){
/* Read a series of measures including all the lines starting with \idem
* that follows the object
*/
             read_series_of_measures_gili(in_line, fp_calib_table, 
                                     object_name, comp_name, 
                                     name1, comp1_is_AB, epoch_o, rho_o, 
                                     theta_o, err_rho_o, err_theta_o, nmeas);
            }
    } /* EOF case name1 == name2 */
  } /* EOF isdigit(in_line[0]) */
  } /* EOF fgets() */
} /* EOF while */

if(*nmeas == 0) {
fprintf(stderr, "get_measures_from_CALIB_table/Object >%s< comp=>%s< not found in calibrated Latex table!\n", object_name, comp_name);
fprintf(stderr, "(name1=%s icol=%d)\n", name1, icol);
status = -1;
} else {
status = 0;
}

fclose(fp_calib_table);
return(status);
}
/******************************************************************************
* Read the full name in a LaTeX table
00550+2338  &  STF73  &  755  & 2010.050 & R  & 20 & 1.007 & 0.008 & 322.7\rlap{$^*$} & 0.5 & 1 &   \\
******************************************************************************/
int read_full_name_from_CALIB_line(char *in_line, char *name2, 
                                   char *comp_name2, int *object2_is_ADS,
                                   int *comp2_is_AB)
{
int status, verbose_if_error = 0;
char *pc1, name22[40];

/* Get the name in the 3rd column if ADS name, or second column otherwise: */
    status = latex_get_column_item(in_line, name2, 3, verbose_if_error);
    jlp_trim_string(name2, 40);
    if(status == 0 && (strncmp(name2,"\\nodata", 7) || (*name2 == '\0'))) {
      *object2_is_ADS = 1;
    } else {
      *object2_is_ADS = 0; 
      status = latex_get_column_item(in_line, name2, 2, verbose_if_error);
      if(status) {
        fprintf(stderr, "read_full_name_from_CALIB_line/Fatal error: in_line=%s\n", 
                in_line);
        exit(-1);
        }
    }

      jlp_compact_string(name2, 40);
/* Get discover name */
      status = latex_get_column_item(in_line, name22, 2, verbose_if_error);
      jlp_compact_string(name22, 40);
/* Compagnon name (AB, Aa, BC, etc) from discover name 
* (after the numbers) */
      comp_name2[0] = '\0';
      if(status) {
          fprintf(stderr,"get_full_name_from_CALIB_table/latex_get_column_item/Fatal error reading in_line (name22=%s)\n", name22);
          exit(-1);
          }
       if(*name22) {
          pc1 = name22;
          while(*pc1  && (isalpha(*pc1) || *pc1 == ' ')) pc1++;
          while(*pc1  && (isdigit(*pc1) || *pc1 == ' ')) pc1++;
          strcpy(comp_name2, pc1);
          *pc1 = '\0';
        }
        jlp_compact_string(comp_name2, 40);
        *comp2_is_AB = 0;
        if((comp_name2[0] == '\0') || !strcmp(comp_name2,"AB")
           || !strncmp(comp_name2,"Aa-B",4)) *comp2_is_AB = 1;
return(0);
}
/******************************************************************************
* Read the full name in a LaTeX table (Gili's format)
00550+2338  &  STF73  &  2010.050 & 2 & 1.007 & 0.008 & 322.7\rlap{$^*$} & 0.5 & 1.23 &   \\
******************************************************************************/
int read_full_name_from_CALIB_line_gili(char *in_line, char *name2, 
                                        char *comp_name2, int *comp2_is_AB)
{
int status, verbose_if_error = 0;
char *pc1, name22[40];

/* Get the name in the second column: */
  status = latex_get_column_item(in_line, name2, 2, verbose_if_error);
  if(status) {
    fprintf(stderr, "read_full_name_from_CALIB_line/Fatal error: in_line=%s\n", 
            in_line);
    exit(-1);
    }
  jlp_compact_string(name2, 40);
/* Get discover name */
  status = latex_get_column_item(in_line, name22, 2, verbose_if_error);
  jlp_compact_string(name22, 40);
/* Compagnon name (AB, Aa, BC, etc) from discover name 
* (after the numbers) */
  comp_name2[0] = '\0';
  if(status) {
      fprintf(stderr,"read_full_name_from_CALIB_table_gili/latex_get_column_item/Fatal error reading in_line (name22=%s)\n", name22);
      exit(-1);
      }
  if(*name22) {
      pc1 = name22;
      while(*pc1  && (isalpha(*pc1) || *pc1 == ' ')) pc1++;
      while(*pc1  && (isdigit(*pc1) || *pc1 == ' ')) pc1++;
      strcpy(comp_name2, pc1);
      *pc1 = '\0';
    }
  jlp_compact_string(comp_name2, 40);
  *comp2_is_AB = 0;
  if((comp_name2[0] == '\0') || !strcmp(comp_name2,"AB")
         || !strncmp(comp_name2,"Aa-B",4)) *comp2_is_AB = 1;
   

return(0);
}
/****************************************************************************
* Read a series of measures with all the lines starting with \idem
* that follows the object in a LaTeX calib table
* 
****************************************************************************/
static int read_series_of_measures(char *in_line, FILE *fp_calib_table, 
                                   char *object_name, 
                                   char *comp_name, char *name,
                                   int comp_is_AB, double *epoch_o, 
                                   double *rho_o, double *theta_o,
                                   double *err_rho_o, double *err_theta_o, 
                                   int *nmeas)
{
register int kk;
int comp2_is_AB, object2_is_ADS;
int line_is_OK, verbose_if_error = 0, status;
char buffer[80], name1[40], comp_name1[40], name2[40], comp_name2[40];
char comp_really_compacted[40], comp2_really_compacted[40];

strcpy(name1, name);
strcpy(comp_name1, comp_name);
jlp_compact_string(name1, 40);
jlp_compact_string(comp_name1, 40);

/* Do while line_is_OK: */
kk = *nmeas;
line_is_OK = 1;
do {
status = read_measures_from_CALIB_line(in_line, &epoch_o[kk],
                                       &rho_o[kk], &err_rho_o[kk],
                                       &theta_o[kk], &err_theta_o[kk]);
if(status > 0) { 
  fprintf(stderr, "Object found: object_name=%s %s (kk=%d) but error retrieving data\n",
       object_name, comp_name, kk); 
  fprintf(stderr, ">%s<\n", in_line);
  break;
  }
else if(status < 0) { 
  fprintf(stderr, "read_series_of_measures_from_CALIB_table/Fatal error reading calibrated file \n"); 
  fprintf(stderr, ">%s<\n", in_line);
  exit(-1);
  }

 
#ifdef DEBUG
printf("Object found: object_name=%s %s (kk=%d) rho=%.3f (+/-%.3f) theta=%.2f (+/-%.2f) epoch=%.3f \n",
       object_name, comp_name, kk, rho_o[kk], err_rho_o[kk], theta_o[kk], err_theta_o[kk], 
       epoch_o[kk]);
#endif

/*If OK, validate read measures by incrementing kk : */
  kk++;
/* Assume by default that data contained in next line is not good: */
  line_is_OK = 0;

/* Go to next line and check if there are other measurements of the same object:
*/
  if(fgets(in_line,300,fp_calib_table)) {
/* Check if line starting with \idem : */
   status = latex_get_column_item(in_line, buffer, 1, verbose_if_error);
   if(status) {
    fprintf(stderr, "Fatal error reading LaTeX table: in_line=>%s<\n", in_line);
    exit(-1);
   }
   jlp_compact_string(buffer, 80);
     if(!strncmp(buffer, "\\idem", 5)) {
       line_is_OK = 5;
       } else if(isdigit(buffer[0])) {
/* Decode the name of the next line in LaTeX file: */
       read_full_name_from_CALIB_line(in_line, name2, comp_name2, 
                                   &object2_is_ADS, &comp2_is_AB);
       if(!strcmp(name1, name2)) {
       jlp_really_compact_companion(comp_name2, comp2_really_compacted, 40);
       jlp_really_compact_companion(comp_name1, comp_really_compacted, 40);
/* If component are the same, add the corresponding measurements */
       if(!strcmp(comp_name1, comp_name2) 
            || (comp_is_AB && comp2_is_AB)
            || (*comp_really_compacted != '\0' && 
               !strcmp(comp_really_compacted, comp2_really_compacted))){
       line_is_OK = 2;
       }
       } /* name1 = name2 */
    } /* not \idem case */
   } 
  } while (line_is_OK);
/* Number of measurements: */
*nmeas += kk;
return(0);
}
/****************************************************************************
* Read a series of measures with all the lines starting with \idem
* that follows the object in a LaTeX calib table
* 
****************************************************************************/
static int read_series_of_measures_gili(char *in_line, FILE *fp_calib_table, 
                                   char *object_name, 
                                   char *comp_name, char *name,
                                   int comp_is_AB, double *epoch_o, 
                                   double *rho_o, double *theta_o,
                                   double *err_rho_o, double *err_theta_o, 
                                   int *nmeas)
{
register int kk;
int comp2_is_AB, object2_is_ADS;
int line_is_OK, verbose_if_error = 0, status;
char buffer[80], name1[40], comp_name1[40], name2[40], comp_name2[40];
char comp_really_compacted[40], comp2_really_compacted[40];

strcpy(name1, name);
strcpy(comp_name1, comp_name);
jlp_compact_string(name1, 40);
jlp_compact_string(comp_name1, 40);

/* Do while line_is_OK: */
kk = *nmeas;
line_is_OK = 1;
do {
status = read_measures_from_CALIB_line_gili(in_line, &epoch_o[kk],
                                       &rho_o[kk], &err_rho_o[kk],
                                       &theta_o[kk], &err_theta_o[kk]);
if(status > 0) { 
  fprintf(stderr, "Object found: object_name=%s %s (kk=%d) but error retrieving data\n",
       object_name, comp_name, kk); 
  fprintf(stderr, ">%s<\n", in_line);
  break;
  }
else if(status < 0) { 
  fprintf(stderr, "read_series_of_measures_from_CALIB_table_gili/Fatal error reading calibrated file \n"); 
  fprintf(stderr, ">%s<\n", in_line);
  exit(-1);
  }

 
#ifdef DEBUG
printf("Object found: object_name=%s %s (kk=%d) rho=%.3f (+/-%.3f) theta=%.2f (+/-%.2f) epoch=%.3f \n",
       object_name, comp_name, kk, rho_o[kk], err_rho_o[kk], theta_o[kk], err_theta_o[kk], 
       epoch_o[kk]);
#endif

/*If OK, validate read measures by incrementing kk : */
  kk++;
/* Assume by default that data contained in next line is not good: */
  line_is_OK = 0;

/* Go to next line and check if there are other measurements of the same object:
*/
  if(fgets(in_line,300,fp_calib_table)) {
/* Check if line starting with \idem : */
   status = latex_get_column_item(in_line, buffer, 1, verbose_if_error);
   if(status) {
    fprintf(stderr, "Fatal error reading LaTeX table: in_line=>%s<\n", in_line);
    exit(-1);
   }
   jlp_compact_string(buffer, 80);
     if(!strncmp(buffer, "\\idem", 5)) {
       line_is_OK = 5;
       } else if(isdigit(buffer[0])) {
/* Decode the name of the next line in LaTeX file: */
       read_full_name_from_CALIB_line_gili(in_line, name2, comp_name2, 
                                           &comp2_is_AB);
// JLP2020: reduce the name if comp_is_AB is true:
       if(comp2_is_AB == 1) remove_AB_from_object_name(name2);

       if(!strcmp(name1, name2)) {
       jlp_really_compact_companion(comp_name2, comp2_really_compacted, 40);
       jlp_really_compact_companion(comp_name1, comp_really_compacted, 40);
/* If component are the same, add the corresponding measurements */
       if(!strcmp(comp_name1, comp_name2) 
            || (comp_is_AB && comp2_is_AB)
            || (*comp_really_compacted != '\0' && 
               !strcmp(comp_really_compacted, comp2_really_compacted))){
       line_is_OK = 2;
       }
       } /* name1 = name2 */
    } /* not \idem case */
   } 
  } while (line_is_OK);
/* Number of measurements: */
*nmeas += kk;
return(0);
}
/***************************************************************************
* Read data from a line of the calibrated Latex table 
* that was generated by "latex_calib.c"
*
* Example of syntax: 
* 18384+0850 & HU 18 & 115 & 2007.690 & R & 20 & 0.483 & 0.013 & 128.1 & 0.4 & 1 &   \\
***************************************************************************/
int read_measures_from_CALIB_line(char *in_line, double *epoch, 
                                   double *rho, double *err_rho,
                                   double *theta, double *err_theta)
{
char buffer[120];
int status, verbose_if_error = 1;

*epoch = 0.; *rho = 0.; *err_rho = 0.; 
*theta = 0.; *err_theta = 0.;

/* Read epoch in column 4: */
status = latex_get_column_item(in_line, buffer, 4, verbose_if_error);
  if(status || (sscanf(buffer,"%lf", epoch) != 1)) {
  fprintf(stderr,"read_measures_from_CALIB_line/Error reading epoch: %s\n", in_line);
  return(-1); 
  } 
/* Read rho in column 7: */
status = latex_get_column_item(in_line, buffer, 7, verbose_if_error);
  if(status) {
  fprintf(stderr,"read_measures_from_CALIB_line/Error reading rho: %s\n", in_line);
  return(-1); 
  } else {
   jlp_trim_string(buffer, 120);
   if(!strncmp(buffer,"\\nodata", 7)) {
    fprintf(stderr,"read_measures_from_CALIB_line/WARNING: Missing rho measurement: %s\n", in_line);
    return(+1); 
   } else if(sscanf(buffer,"%lf", rho) != 1) {
    fprintf(stderr,"read_measures_from_CALIB_line/Error reading rho: %s\n", in_line);
    return(-1); 
   }
  }

/* Read err_rho in column 8: */
status = latex_get_column_item(in_line, buffer, 8, verbose_if_error);
  if(status || (sscanf(buffer,"%lf", err_rho) != 1)) {
  fprintf(stderr,"read_measures_from_CALIB_line/Error reading err_rho: %s\n", in_line);
  return(-1); 
  }
/* Read theta in column 9: */
status = latex_get_column_item(in_line, buffer, 9, verbose_if_error);
  if(status || (sscanf(buffer,"%lf", theta) != 1)) {
  fprintf(stderr,"read_measures_from_CALIB_line/Error reading theta: %s\n", in_line);
  return(-1); 
  } 
/* Read err_theta in column 10: */
status = latex_get_column_item(in_line, buffer, 10, verbose_if_error);
  if(status || (sscanf(buffer,"%lf", err_theta) != 1)) {
  fprintf(stderr,"read_measures_from_CALIB_line/Error reading err_theta: %s\n", in_line);
  return(-1); 
  } 
return(0);
} 
/***************************************************************************
* Read data from a line of the calibrated Latex table 
* that was generated by "latex_calib.c"
*
* Example of syntax: 
* 18384+0850 & HU 18 & 2007.690 & 2 & 0.483 & 0.013 & 128.1 & 0.4 & 1.23 &   \\
***************************************************************************/
int read_measures_from_CALIB_line_gili(char *in_line, double *epoch, 
                                   double *rho, double *err_rho,
                                   double *theta, double *err_theta)
{
char buffer[120];
int status, verbose_if_error = 1;

*epoch = 0.; *rho = 0.; *err_rho = 0.; 
*theta = 0.; *err_theta = 0.;

/* Read epoch in column 3: */
status = latex_get_column_item(in_line, buffer, 3, verbose_if_error);
  if(status || (sscanf(buffer,"%lf", epoch) != 1)) {
  fprintf(stderr,"read_measures_from_CALIB_line/Error reading epoch: %s\n", in_line);
  return(-1); 
  } 
/* Read rho in column 5: */
status = latex_get_column_item(in_line, buffer, 5, verbose_if_error);
  if(status) {
  fprintf(stderr,"read_measures_from_CALIB_line/Error reading rho: %s\n", in_line);
  return(-1); 
  } else {
   jlp_trim_string(buffer, 120);
   if(!strncmp(buffer,"\\nodata", 7)) {
    fprintf(stderr,"read_measures_from_CALIB_line/WARNING: Missing rho measurement: %s\n", in_line);
    return(+1); 
   } else if(sscanf(buffer,"%lf", rho) != 1) {
    fprintf(stderr,"read_measures_from_CALIB_line/Error reading rho: %s\n", in_line);
    return(-1); 
   }
  }

/* Read err_rho in column 6: */
status = latex_get_column_item(in_line, buffer, 6, verbose_if_error);
  if(status || (sscanf(buffer,"%lf", err_rho) != 1)) {
  fprintf(stderr,"read_measures_from_CALIB_line/Error reading err_rho: %s\n", in_line);
  return(-1); 
  }
/* Read theta in column 7: */
status = latex_get_column_item(in_line, buffer, 7, verbose_if_error);
  if(status || (sscanf(buffer,"%lf", theta) != 1)) {
  fprintf(stderr,"read_measures_from_CALIB_line/Error reading theta: %s\n", in_line);
  return(-1); 
  } 
/* Read err_theta in column 8: */
status = latex_get_column_item(in_line, buffer, 8, verbose_if_error);
  if(status || (sscanf(buffer,"%lf", err_theta) != 1)) {
  fprintf(stderr,"read_measures_from_CALIB_line/Error reading err_theta: %s\n", in_line);
  return(-1); 
  } 
return(0);
} 
/***************************************************************************
* Read the object name from a line of the calibrated Latex table 
* that was generated by "latex_calib.c"
*
* Example of syntax: 
* 18384+0850 & STF 138 AB & 115 & 2007.690 & R & 20 & 0.483 & 0.013 & 128.1 & 0.4 & 1 &   \\
*
* INPUT:
* in_line: line read from CALIB table
*
* OUTPUT:
* wds_name: WDS name (e.g. 00013+3351)
* ads_name: ADS 234, ADS 2375 if in ADS...)
* discov_name: STF 224, COU 34, etc
* comp_name: (AB, Aa, etc)
***************************************************************************/
int read_object_name_from_CALIB_line(char *in_line, char *wds_name, 
                                     char *discov_name, char *comp_name,
                                     char *ads_name)
{
char buffer[120], *pc;
int ads_number, verbose_if_error = 1;

wds_name[0] = '\0';
discov_name[0] = '\0';
ads_name[0] = '\0';
comp_name[0] = '\0';

/* Read WDS name in column 1: */
latex_get_column_item(in_line, wds_name, 1, verbose_if_error);

/* Read discoverer's name in column 2: */
latex_get_column_item(in_line, discov_name, 2, verbose_if_error);

/* JLP2009: I compact the name */
jlp_compact_string(discov_name, 40);

/* Extract the companion name from the discoverer's name
if present: */
strcpy(buffer, discov_name);
pc = buffer;
while(isalpha(*pc) || *pc == ' ') pc++;
while(isdigit(*pc) || *pc == ' ') pc++;
strcpy(comp_name, pc);

/* Now remove the companion name from the discoverer's name: */
pc = discov_name;
while(isalpha(*pc) || *pc == ' ') pc++;
while(isdigit(*pc)) pc++;
*pc = '\0';

/* Read ADS number in column 3: */
ads_name[0] = '\0';
latex_get_column_item(in_line, buffer, 3, verbose_if_error);
  if(sscanf(buffer,"%d", &ads_number) == 1) {
    sprintf(ads_name, "ADS %d", ads_number);
  }

return(0);
} 
/***************************************************************************
* Read the object name from a line of the calibrated Latex table 
* that was generated by "latex_calib.c" in gili's format (without ADS name)
*
* Example of syntax: 
* 18384+0850 & STF 138 AB & 115 & 2007.690 & R & 20 & 0.483 & 0.013 & 128.1 & 0.4 & 1 &   \\
*
* INPUT:
* in_line: line read from CALIB table
*
* OUTPUT:
* wds_name: WDS name (e.g. 00013+3351)
* discov_name: STF 224, COU 34, etc
* comp_name: (AB, Aa, etc)
***************************************************************************/
int read_object_name_from_CALIB_line_gili(char *in_line, char *wds_name, 
                                     char *discov_name, char *comp_name)
{
char buffer[120], *pc;
int verbose_if_error = 1;

wds_name[0] = '\0';
discov_name[0] = '\0';
comp_name[0] = '\0';

/* Read WDS name in column 1: */
latex_get_column_item(in_line, wds_name, 1, verbose_if_error);

/* Read discoverer's name in column 2: */
latex_get_column_item(in_line, discov_name, 2, verbose_if_error);

/* JLP2009: I compact the name */
jlp_compact_string(discov_name, 40);

/* Extract the companion name from the discoverer's name
if present: */
strcpy(buffer, discov_name);
pc = buffer;
while(isalpha(*pc) || *pc == ' ') pc++;
while(isdigit(*pc) || *pc == ' ') pc++;
strcpy(comp_name, pc);

/* Now remove the companion name from the discoverer's name: */
pc = discov_name;
while(isalpha(*pc) || *pc == ' ') pc++;
while(isdigit(*pc)) pc++;
*pc = '\0';

return(0);
} 
/***********************************************************************
* Retrieve the residuals for a given object at a given epoch
* from the residual Latex table 
* JLP2011: retrieve as many lines as there are orbits for this object
* in the RESID file
*
* INPUT:
* resid_fname: file name of the residual LateX table
* ref_slength: maxi length of the references
* object_name: name of object (either ADS 123 or COU 432)
* comp_name: name of the companion (i.e. AB, Aa, CD, etc) 
* epoch_o: epoch of observation 
* nmax_orbits: maximum number of orbits that are allowed for output
*
* OUTPUT:
* orbit_ref: orbit reference (ref. of the publication of the orbit)
* declared as orbit_ref[ref_length * nmax_orbits]  (character string array)
* rho_o_c (arcsec), theta_o_c (deg): residuals O-C in rho and theta 
* quadrant_discrep: "\rlap{$^Q}$" when discrepancy between measure and orbit
* norbits_found: number of orbits found
***********************************************************************/
int get_values_from_RESID_table(char *resid_fname, char *object_name,
                                char *comp_name, double epoch_o, double rho_o,
                                char *orbit_ref, int ref_slength,
                                double *rho_o_c,
                                double *theta_o_c, char *quadrant_discrep,
                                int *norbits_found, int nmax_orbits)
{
char in_line[128], object_name0[128], *pc;
double epoch_o_c, rho_val0, rho0, theta0;
int status, iline, verbose_if_error = 1, kk;
char orbit_ref0[40], comp_name0[40], comp_name1[40], object_name1[40];
char old_orbit_ref[40];
int comp_is_AB, comp0_is_AB;
FILE *fp_resid_table;

*rho_o_c = -100.;
*theta_o_c = -100.;
for(kk = 0; kk < nmax_orbits; kk++) strcpy(&orbit_ref[kk*ref_slength], "");

strcpy(quadrant_discrep, "");

/* Generate full object name:
   ADS 123 AB, STF 128 AB
*/
/* Remove all blanks: */
strcpy(object_name1, object_name);
jlp_compact_string(object_name1, 40);
strcpy(comp_name1, comp_name);
jlp_compact_string(comp_name1, 40);
if(!strcmp(comp_name1, "AB") || *comp_name1 == '\0') comp_is_AB = 1;
 else comp_is_AB = 0;


/* Open LaTeX RESID table: */
if((fp_resid_table = fopen(resid_fname, "r")) == NULL) {
  fprintf(stderr, "get_values_from_RESID_table/Error opening %s\n", 
          resid_fname);
  return(-1);
 }

/* Scan all the file looking for the object name */ 
strcpy(old_orbit_ref,"none");
iline = 0;
kk = 0;
while(!feof(fp_resid_table)){
  if(fgets(in_line, 128, fp_resid_table)) {
  iline++;
  if(in_line[0] != '%' && in_line[0] != '\\' && in_line[0] != ' '
     && in_line[0] != '&') {
    status = latex_get_column_item(in_line, object_name0, 1, verbose_if_error);
    if(status) {
      fprintf(stderr, "Fatal error/Bad syntax of %s in line %d\n",
              resid_fname, iline); 
      exit(-1);
      }
    jlp_compact_string(object_name0, 40);
/* Extract the companion and cut the object name: */
      comp_name = '\0';
       if(*object_name0) {
          pc = object_name0;
          while(*pc  && (isalpha(*pc) || *pc == ' ')) pc++;
          while(*pc  && (isdigit(*pc) || *pc == ' ')) pc++;
          strcpy(comp_name0, pc);
          *pc = '\0';
        }
      jlp_compact_string(comp_name0, 40);
      if(!strcmp(comp_name0, "AB") || *comp_name0 == '\0') comp0_is_AB = 1;
      else comp0_is_AB = 0;
 
/*
printf("object_name=%s< object_name0=%s<\n", object_name1, object_name0);
printf("comp_name=%s< comp_name0=%s<\n", comp_name1, comp_name0);
printf("comp_is_AB=%d< comp0_is_AB=%d<\n", comp_is_AB, comp0_is_AB);
*/
/* Process only the meaningful lines (skipping the header...)*/
      if( 
       (!strcmp(object_name1, object_name0) && comp0_is_AB && comp_is_AB )
       || (!strcmp(object_name1, object_name0) && !strcmp(comp_name0, comp_name1))
        ) {
         status = read_values_from_RESID_line(in_line, &epoch_o_c, 
                                             &rho_val0, &rho0, &theta0);
        if(status) {
           fprintf(stderr, "Fatal error/Bad syntax of %s in line %d\n",
                   resid_fname, iline); 
           exit(-1);
           }
        latex_get_column_item(in_line, orbit_ref0, 2, verbose_if_error);
/* Truncate to the first word only 
* Cou1973b - Couteau (1973b)    => Cou1973b
*/
        jlp_trim_string(orbit_ref0, 60);
/* If the epochs are equal, the rho are equal, and orbit_ref is new, 
*  load the values of this orbit: */
        if( (ABS(epoch_o_c - epoch_o) < 0.002 ) 
            && (ABS(rho_o - rho_val0) < 0.002) 
           && (strcmp(orbit_ref0, old_orbit_ref) != 0)) {
printf("get_values_from_RESID_table/DEBUG: object=>%s epoch_o=%f epoch_o_c=%f (old_orbit_ref=%s) rho_val0=%f\n", 
        object_name1, epoch_o, epoch_o_c, old_orbit_ref, rho_val0);

          strcpy(old_orbit_ref, orbit_ref0);
          rho_o_c[kk] = rho0;
          theta_o_c[kk] = theta0;
          if(is_in_line(in_line, "$^Q$")){
            strcpy(quadrant_discrep, "\\rlap{$^Q$}");
            }

          strcpy(&orbit_ref[kk * ref_slength], orbit_ref0);

          kk++;
/*      break; */
/* JLP2011: read subsequent lines now to allow multiple orbits... */
        }
    } /* EOF !strcmp(object_name_0)... */
  } /* EOF in_line[0] != '%' */
  } /* EOF fgets() */
} /* EOF while */

*norbits_found = kk;

fclose(fp_resid_table);
return(0);
}
/***************************************************************************
* Read data from resdiual Latex table as generated by residuals_1.c
*
* Example:
* ADS 9494 AB & Sod1999 - Soderhjelm (1999) & 2007.534 & 1.842 & 0.01 & 0.10 \\
***************************************************************************/
int read_values_from_RESID_line(char *in_line, double *epoch_o_c, 
                                double *rho_val, double *rho_o_c, 
                                double *theta_o_c)
{
char buffer[120];
int status, verbose_if_error = 1;

*epoch_o_c = 0.; *rho_o_c = 0.; *theta_o_c = 0.; 

/* Read epoch in column 3: */
status = latex_get_column_item(in_line, buffer, 3, verbose_if_error);
  if(status || (sscanf(buffer,"%lf", epoch_o_c) != 1)) {
  fprintf(stderr,"Error reading epoch: %s\n", in_line);
  return(-1); 
  } 
/* Read rho in column 4: */
status = latex_get_column_item(in_line, buffer, 4, verbose_if_error);
  if(status || (sscanf(buffer,"%lf", rho_val) != 1)) {
  fprintf(stderr,"Error reading rho value: %s\n", in_line);
  return(-1); 
  } 
/* Read rho_O-C in column 5: */
status = latex_get_column_item(in_line, buffer, 5, verbose_if_error);
  if(status || (sscanf(buffer,"%lf", rho_o_c) != 1)) {
  fprintf(stderr,"Error reading rho_o_c: %s\n", in_line);
  return(-1); 
  }
/* Read theta in column 6: */
status = latex_get_column_item(in_line, buffer, 6, verbose_if_error);
  if(status || (sscanf(buffer,"%lf", theta_o_c) != 1)) {
  fprintf(stderr,"Error reading theta: %s\n", in_line);
  return(-1); 
  } 
return(0);
} 
/************************************************************
* Get character string in the icol th column of a LaTeX table
*
* INPUT:
*  in_line: string corresponding to a full line of a Latex table
*  icol: column number
*
* OUTPUT:
*  item: character string corresponding to the column #icol
* 
************************************************************/
int latex_get_column_item(char *in_line, char *item, int icol, 
                          int verbose_if_error)
{
int ic;
char *pc, buffer[512];
 pc = in_line;

/* Get to the icol th column: */
 ic = 1;
 while(ic != icol) {
   while(*pc && *pc != '&' && strncmp(pc, "\\cr", 3) 
        && strncmp(pc, "\\\\", 2)) pc++; 
    if(*pc == '&') {
      pc++;
      ic++;
     } else {
     if(verbose_if_error) {
        fprintf(stderr, "latex_get_column_item/Error: column #%d not found in >%s<\n",
             icol, in_line);
        }
     return(-1);
     }
 }

/* icol has been reached: */
 strcpy(buffer, pc);
 pc = buffer;
/* Go to the next column or the end of line: */
 while(*pc && *pc != '&' && strncmp(pc, "\\cr", 3)
        && strncmp(pc, "\\\\", 2)) pc++; 
 *pc = '\0';
 strcpy(item, buffer);

return(0);
}
/************************************************************************
* Extract ADS name from object name
*
* INPUT:
* object_name: name in the PISCO catalog
*
* OUTPUT:
* ADS_name: string containing only the ADS number
************************************************************************/
int ADS_name_from_object_name(char *object_name, char *ADS_name)
{
char *pc, buffer[40], ADS_buffer[40];
int i;

strncpy(buffer, object_name, 40);
buffer[39] = '\0';

pc = buffer;
i = 0;
while(*pc && *pc != 'A') pc++;
/* If ADS is present, retrieve the ADS number */
if(!strncmp(pc, "ADS", 3)) {
  pc += 3;
  while(*pc) ADS_buffer[i++] = *(pc++);
  }
ADS_buffer[i] = '\0';

/* Extract and clean the string following "ADS ": */
jlp_trim_string(ADS_buffer, 40);
pc = ADS_buffer;
i = 0;
while(*pc && isdigit(*pc)) pc++;
*pc = '\0';
strncpy(ADS_name, ADS_buffer, 20);
ADS_name[19] = '\0';

return(0);
}
/******************************************************************
* Convert companion name to avoid problems when finding orbits
*
* e.g.  Bb-C  -> BC
*       Aa-Bb -> AB
*       AaBb -> AB
*       Cc-D -> CD
*       Cc,D -> CD
*       AB-C -> AC
*       AC   -> AC
*****************************************************************/
int jlp_really_compact_companion(char *name_in, char *name_out, int length)
{
char *pc;
int k, separation_found, main_found;

name_in[length-1] = '\0';

/******* Look for separation */
pc = name_in;
separation_found = 0;
k = 0;
while(*pc) { 
  if(*pc == '-' || *pc == ',') { 
  separation_found = 1;
  } else if (isupper(*pc)) {
   name_out[k++] = *pc;
  }
  pc++;
  }
name_out[k] = '\0';
/* Job done if no separation: */
if(!separation_found) return(0);

/******* Look for main components if separation: */
pc = name_in;
main_found = 0;
k = 0;
while(*pc) {
  if(*pc == '-' || *pc == ',') { 
   main_found = 0;
  } else if(!main_found && isupper(*pc)) { 
   name_out[k++] = *pc;
   main_found = 1;
   }
  pc++;
  }
name_out[k] = '\0';

return(0);
}
/********************************************************************
* Remove "AB" from object_name1
********************************************************************/
int remove_AB_from_object_name(char *object_name1)
{
char buffer[64], str_to_search[16], *position_str;
int position;

strcpy(buffer, object_name1);
strcpy(str_to_search, "AB");

position_str = strstr(buffer, str_to_search);
if(position_str != NULL) {
    position = position_str - buffer;  
    buffer[position] = '\0';
    strcpy(object_name1, buffer);
}
return 0;
}

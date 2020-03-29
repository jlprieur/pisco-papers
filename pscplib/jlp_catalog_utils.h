/************************************************************************
* "jlp_catalog_utils.h"
* Set of routines used to read catalogs: WDS, OC6, Hipparcos, 
* PISCO catalog, CALIBrated table (from latex_calib.c) and Latex tables
*
* JLP 
* Version 21/04/2009
*************************************************************************/
#ifndef _jlp_catalog_utils_h /* BOF sentry */
#define _jlp_catalog_utils_h

#include <stdio.h>
#include <stdlib.h>                  /* exit() */
#include <string.h>
#include <ctype.h>                   /* isprint... */
#include <math.h>
#include <time.h>                    /* date */

#ifndef ABS
#define ABS(a) ((a) < 0.0  ? (-(a)) : (a))
#endif
#ifndef PI
#define PI 3.14159265
#endif
#define DEGTORAD   (PI/180.00)

/* Declaring linkage specification to have "correct names"
* that can be linked with C programs */

#ifdef __cplusplus
extern "C" {
#endif

int is_in_line(char *line_buffer2, char *str1);
int latex_get_column_item(char *in_line, char *item, int icol, 
                          int verbose_if_error);
int get_values_from_RESID_table(char *resid_fname, char *object_name,
                                char *comp_name, double epoch_o, double rho_o,
                                char *orbit_ref, int ref_slength,
                                double *rho_o_c, 
                                double *theta_o_c, char *quadrant_discrep,
                                int *norbits_found, int nmax_orbits);
int read_values_from_RESID_line(char *in_line, double *epoch_o_c,
                                double *rho_val,
                                double *rho_o_c, double *theta_o_c);
int get_measures_from_CALIB_table(char *calib_fname, char *object_name,
                                  char *comp_name, double *epoch_o, 
                                  double *rho_o, double *theta_o,
                                  double *err_rho_o, double *err_theta_o, 
                                  int *nmeas);
int get_measures_from_CALIB_table_gili(char *calib_fname, char *object_name,
                                  char *comp_name, double *epoch_o, 
                                  double *rho_o, double *theta_o,
                                  double *err_rho_o, double *err_theta_o, 
                                  int *nmeas);
int read_measures_from_CALIB_line(char *line_buffer, double *epoch, 
                                   double *rho, double *err_rho,
                                   double *theta, double *err_theta);
int read_measures_from_CALIB_line_gili(char *line_buffer, double *epoch, 
                                   double *rho, double *err_rho,
                                   double *theta, double *err_theta);
int read_object_name_from_CALIB_line(char *in_line, char *wds_name,
                                     char *discov_name, char *comp_name,
                                     char *ads_name);
int read_object_name_from_CALIB_line_gili(char *in_line, char *wds_name,
                                     char *discov_name, char *comp_name);
int read_full_name_from_CALIB_line(char *in_line, char *name2,
                                          char *comp_name2, int *object2_is_ADS,
                                          int *comp2_is_AB);
int read_full_name_from_CALIB_line_gili(char *in_line, char *name2,
                                          char *comp_name2, int *comp2_is_AB);

int ADS_name_from_object_name(char *object_name, char *ADS_name);
int jlp_really_compact_companion(char *name_in, char *name_out, int length);
int remove_AB_from_object_name(char *object_name1);

#ifdef __cplusplus
}
#endif


#endif

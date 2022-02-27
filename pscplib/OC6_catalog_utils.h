/************************************************************************
* "OC6_catalog_utils.h"
* Set of routines used to read OC6 catalog
*
* JLP 
* Version 21/06/2015
*************************************************************************/
#ifndef _OC6_catalog_utils_h /* BOF sentry */
#define _OC6_catalog_utils_h

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

int get_orbit_from_OC6_list_gili(char *in_line, int iline, int is_master_file,
                       char *WDS_name, char *discov_name,
                       char *comp_name, char *object_name, char *author,
                       double *Omega_node, double *omega_peri, double *i_incl, 
                       double *e_eccent, double *T_periastron, double *Period, 
                       double *a_smaxis, double *mean_motion, 
                       double *orbit_equinox, int *orbit_grade);
int line_extraction_from_OC6_catalog_gili(char *OC6_fname, int is_master_file,
                                     char *discov_name, char *comp_name, 
                                     FILE *fp_out, int *found, 
                                     int *candidate_found, 
                                     int norbits_per_object);
int get_orbit_from_OC6_list(char *in_line, int iline, int is_master_file,
                       char *WDS_name, char *ADS_name, char *discov_name,
                       char *comp_name, char *object_name, char *author,
                       double *Omega_node, double *omega_peri, double *i_incl, 
                       double *e_eccent, double *T_periastron, double *Period, 
                       double *a_smaxis, double *mean_motion, 
                       double *orbit_equinox, int *orbit_grade);
int line_extraction_from_OC6_catalog(char *OC6_fname, int is_master_file,
                                     char *ads_name, char *discov_name,
                                     char *comp_name, FILE *fp_out, int *found, 
                                     int *candidate_found, 
                                     int norbits_per_object);
int get_name_from_OC6_line(char *in_line, char *OC6_ads_name,
                           char *OC6_discov_name, char *OC6_comp_name);
int get_name_from_OC6_line_gili(char *in_line, char *OC6_discov_name, 
                                char *OC6_comp_name);
int get_OC6_full_reference(char *object_name, char *author,
                           char *OC6_references_fname, char *refer0, 
                           char *refer1);
int save_references_for_residuals(FILE *fp_out_ref1, FILE *fp_out_ref2,
                                  char *object_name, char *author,
                                  char *refer0, char *refer1);
int sort_references(FILE *fp_out_ref1, FILE *fp_out_ref2, char *object_name, 
                    char *author, char *refer0, char *refer1, int n_names);

#ifdef __cplusplus
}
#endif


#endif

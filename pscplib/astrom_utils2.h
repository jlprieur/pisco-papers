/*************************************************************************
* Set of routines to read Latex files with astrometric measurements
* created by Xdisp1 or Wdisp1 (e.g. astrom05a.tex) 
*
* JLP
* Version 10/07/2018
*************************************************************************/
#ifndef _astrom_utils2_h  /* BOF sentry */
#define _astrom_utils2_h  

#include "astrom_utils1.h"

#ifdef __cplusplus
extern "C" {
#endif

int astrom_add_epoch_from_fits_file(FILE *fp_in, FILE *fp_out, 
                                    char *fits_directory);
int astrom_check_if_measurement(char *b_in, char *fits_filename, char *date1,  
                                int *eyepiece1, int *is_measurement);
int astrom_calib_publi(FILE *fp_in, FILE *fp_out, double *calib_scale1,
                       int *calib_eyepiece1, int n_eyepieces1, double theta0, 
                       double sign, int i_filename, int i_date, int i_filter,
                       int i_eyepiece, int i_rho, int i_drho, int i_theta, 
                       int i_dtheta, int i_notes, int comments_wanted, 
                       char *filein, int in_astrom_fmt, int out_calib_fmt);
int astrom_calib_copy(FILE *fp_in, FILE *fp_out, double *calib_scale1,
                      int *calib_eyepiece1, int n_eyepieces1, double theta0,
                      double sign, int i_date, int i_eyepiece, int i_rho,
                      int i_drho, int i_theta, int i_dtheta,
                      int comments_wanted, int input_with_header);
int astrom_add_WDS(FILE *fp_in, FILE *fp_out, char *PISCO_cat, char *WDS_cat); 
int astrom_add_new_object_with_wds_data(char *b_data, OBJECT *obj, int *nobj, 
                                        int i_notes, int in_astrom_fmt);
int astrom_add_new_object_without_wds_data(char *b_data, OBJECT *obj, 
                                           int *nobj, int in_astrom_fmt);
int astrom_read_measures(FILE *fp_in, int comments_wanted, OBJECT *obj, 
                         int *nobj, int i_filename, int i_date, 
                         int i_filter, int i_eyepiece, int i_rho, 
                         int i_drho, int i_theta, int i_dtheta, int i_notes,
                         int with_wds_data, int in_astrom_fmt);
int astrom_mean_for_paper2(OBJECT *obj, int nobj);
int astrom_mean_for_full_table(OBJECT *obj, int nobj);
int astrom_compute_mean_of_two_measures(OBJECT *obj, int io, int jm1, int jm2);
int compute_mean_for_filter(OBJECT *obj, int i, int *j0, int nj0, 
                            int quadrant_value);
int merate_get_full_directory(char *fits_directory, char *date,
                              char *full_directory);
int astrom_get_name_from_2nd_col(char *b_in, char *ads_name, char *discov_name,
                                 char *comp_name, double *year);
#ifdef __cplusplus
}
#endif

#endif /* EOF sentry */

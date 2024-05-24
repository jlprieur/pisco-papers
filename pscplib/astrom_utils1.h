/*************************************************************************
* Set of routines to read Latex files with astrometric measurements
* created by Xdisp1 or Wdisp1 (e.g. astrom05a.tex) 
*
* JLP
* Version 03/03/2020
*************************************************************************/
#ifndef _astrom_utils1_h  /* BOF sentry */
#define _astrom_utils1_h  
#include <stdio.h>
#include <stdlib.h>   /* definition of exit() */
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "astrom_def.h" 

#ifdef __cplusplus
extern "C" {
#endif

int get_calib_scale(double *calib_date1, int *calib_dd1, int *calib_mm1,
                    int *calib_year1, double *calib_scale1,
                    int *calib_eyepiece1, int *n_eyepieces1, double *theta01, 
                    int *sign1, int ncalib1, int ndim, double me_bessel_epoch, 
                    int me_eyepiece, double *scale2, double *theta02, 
                    int *sign2);
int astrom_read_new_object_line_for_wds_or_ads(char *b_data, char *wds_name, 
                            char *discov_name, char *comp_name,
                            char *ads_name, 
                            double *WR, double *WT, double *WY, int i_notes,
                            int gili_format);
int astrom_read_object_data_for_name_only(char *b_data, char *discov_name, 
                            char *comp_name, double *bessel_epoch_year,
                            int gili_format);
int astrom_check_measure(char *b_data, int i_eyepiece, int i_rho, 
                         int i_drho, int i_theta, int i_dtheta);
int astrom_calibrate_measures(OBJECT *obj, int nobj, double *calib_date1, 
                              int *calib_dd1, int *calib_mm1,
                              int *calib_year1, double *calib_scale1,
                              int *calib_eyepiece1, int *n_eyepieces1,
                              double *theta01, int *sign1, int ncalib1,
                              int ndim);
int astrom_calib_data_copy(char *b_data, char *b_out,  
                           double *calib_date1, int *calib_dd1, int *calib_mm1,
                           int *calib_year1, double *calib_scale1,
                           int *calib_eyepiece1, int *n_eyepieces1,
                           double *theta01, int *sign1, int ncalib1,
                           int ndim, int i_date, int i_eyepiece, int i_rho, 
                           int i_drho, int i_theta, int i_dtheta);
int astrom_add_new_measure(char *b_data, OBJECT *obj, int i_obj, 
                           int i_filename, int i_date, int i_filter, 
                           int i_eyepiece, int i_rho, int i_drho, int i_theta, 
                           int i_dtheta, int i_notes, int comments_wanted);
int astrom_read_new_measure(char *b_data, int i_filename, int i_date, 
                            int i_filter, int i_eyepiece, int i_rho, int i_drho,                            int i_theta, int i_dtheta, int i_notes, 
                            char *filename2, char *notes2, 
                            double *BesselEpoch2, 
                            char *filter2, int *eyepiece2, double *rho2, 
                            double *err_rho2, double *theta2, 
                            double *err_theta2);
int astrom_decode_data(char *b_data, char *date, double *BesselEpoch, 
                       double *rho, 
                       double *drho, double *theta, double *dtheta, 
                       int *eyepiece, int i_date, int i_eyepiece, int i_rho, 
                       int i_drho, int i_theta, int i_dtheta);
int astrom_compute_bessel_epoch_value(char *b_data, char *date, 
                                      double *BesselEpoch, int icol);
int astrom_ra_sort_objects(OBJECT *obj, int *index_obj, int nobj);
int astrom_name_sort_objects(OBJECT *obj, int *index_obj, int nobj);
int astrom_read_WDS_CHARA(char *notes, double *WR, double *WT, double * WY);
int astrom_read_quadrant_Q(char *notes, int *quadrant, int *dquadrant,
                           int comments_wanted);
int astrom_read_quadrant_LQ(char *notes, int *quadrant, int *dquadrant,
                            int comments_wanted);
int astrom_quadrant_is_consistent(MEASURE *me);
int astrom_compute_statistics(FILE *fp_out, OBJECT *obj, int nobj, 
                              char *filein);
int astrom_correct_theta_with_WDS_CHARA(OBJECT *ob, MEASURE *me);
int astrom_read_bessel_epoch_from_notes(char *notes, double *BesselEpoch, 
                                        int delete_epoch_item);
int astrom_read_julian_epoch_from_notes(char *notes, double *JulianEpoch, 
                                        int delete_epoch_item);
int astrom_read_dmag_from_notes(char *notes, double *dmag, double *ddmag, 
                                int comments_wanted);
int astrom_preformat_wds_name(char *obj_wds, char *wds_name);
int astrom_check_if_object_name(char *b_in, int *contains_object_name,
                                int *contains_WDS_name);
int astrom_is_record_file(char *filename);
int astrom_is_direct_file(char *filename);
int jlp_split_discov_comp(char *full_name, int fname_length,
                          char *discov_name, char *comp_name);
int check_filter(char *filter, char *b_data);

#ifdef __cplusplus
}
#endif

#endif  /* EOF sentry */

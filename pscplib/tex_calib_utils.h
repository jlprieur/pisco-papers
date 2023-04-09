/*************************************************************************
* Set of routines to read Latex files with astrometric measurements
*
* JLP
* Version 27/10/2021
*************************************************************************/
#ifndef _tex_calib_utils_h  /* BOF sentry */
#define _tex_calib_utils_h  

#include "astrom_utils1.h"

#ifdef __cplusplus
extern "C" {
#endif

int extract_companion_from_name(char *input_object_name, char *object_name0, 
                                char *comp_name0, int *object_len0);
int tex_calib_read_measures(char *filein1, OBJECT *obj1, int *nobj1,
                            int *nunres1, int in_calib_fmt);
int tex_calib_read_measures_from_line(char *b_in, OBJECT *obj, int *nobj, 
                                      int *nunres1, int in_calib_fmt);
int tex_calib_compare(double rho3, double err_rho3, double theta3,
                      double err_theta3, double rho1, double err_rho1,
                      double theta1, double err_theta1, 
                      double max_drho, double max_dtheta,
                      int *compatible_meas31);
int astrom_write_publi_table(FILE *fp_out, int comments_wanted, OBJECT *obj,
                             int *index_obj, int nobj, int tabular_only);
int astrom_write_publi_table_gili(FILE *fp_out, int comments_wanted,
                                  OBJECT *obj, int *index_obj, int nobj,
                                  int tabular_only);
int tex_calib_write_miniheader(FILE *fp_out);
int tex_calib_get_meas_from_object(OBJECT *obj1, int nobj1, 
                                   char *WDSName0, char *ObjectName0,
                                   double epoch0,
                                   char *WDSName1, char *ObjectName1,
                                   double *epoch1, char *filt1, int *eyep1, 
                                   double *rho1, double *err_rho1, 
                                   double *theta1, double *err_theta1);
int tex_calib_get_meas_from_csv_object(OBJECT *obj1, int nobj1, 
                                   char *ObjectName0, double epoch0,
                                   char *ObjectName1, double *epoch1,
                                   int *eyep1, double *rho1, double *err_rho1, 
                                   double *theta1, double *err_theta1,
                                   double *dmag1);

#ifdef __cplusplus
}
#endif

#endif /* EOF sentry */

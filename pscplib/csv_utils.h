/*************************************************************************
* Set of routines to read csv files with astrometric measurements
*
* JLP
* Version 15/01/2023
*************************************************************************/
#ifndef _csv_utils_h  /* BOF sentry */
#define _csv_utils_h  

#include "astrom_utils1.h"

#ifdef __cplusplus
extern "C" {
#endif

int csv_read_gili_measures(char *filein1, OBJECT *obj1, int *nobj1,
                           int nobj_maxi, double scale_mini);
int csv_read_gili_measures_from_line(char *b_data, OBJECT *obj1, int *nobj1,
                                     double scale_mini);
int csv_read_object_from_filename(char *b_data, int i_column, char *out_string);
int csv_read_string(char *b_data, int i_column, char *out_string);
int csv_read_dvalue(char *b_data, int i_column, double *dvalue);
int csv_read_epoch_from_notes(char *b_data, int i_notes, double *epoch1);
int csv_read_dmag_from_notes(char *b_data, int i_notes, double *dmag1);

#ifdef __cplusplus
}
#endif

#endif /* EOF sentry */

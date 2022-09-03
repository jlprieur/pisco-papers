/************************************************************************
* "latex_utils.h"
*
* JLP 
* Version 21/08/2020
*************************************************************************/
#ifndef _latex_utils_h /* BOF sentry */
#define _latex_h

#include "astrom_def.h"

/* Declaring linkage specification to have "correct names"
* that can be linked with C programs */

#ifdef __cplusplus
extern "C" {
#endif

int latex_read_ivalue(char *b_data, int *value, int icol);
int latex_read_dvalue(char *b_data, double *value, int icol, int verbose);
int latex_read_fvalue(char *b_data, float *value, int icol, int verbose);
int latex_read_svalue(char *b_data, char *value, int icol);
int latex_write_dvalue(char *b_data, char *b_out, double value, int icol,
                       int nber_of_decimals);

int latex_get_column_item(char *in_line, char *item, int icol, 
                          int verbose_if_error);
int latex_set_column_item(char *in_line, int in_line_length,
                          char *new_item, int new_item_len, 
                          int icol, int verbose_if_error);
int latex_remove_column(char *in_line, int i_col, int ilen);
int latex_remove_dollar(char *in_string, char *out_string, int len_string);
int latex_add_emptycols_to_ncols(char *in_string, char *out_string, 
                                 int len_string_max, int ncols_max);

#ifdef __cplusplus
}
#endif


#endif

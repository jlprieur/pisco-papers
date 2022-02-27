/************************************************************************
* "WDS_catalog_utils.h"
* To retrieve data from the WDS catalog
*
* JLP 
* Version 12/02/2012
*************************************************************************/
#ifndef __WDS_catalog_utils   /* BOF sentry */
#define __WDS_catalog_utils   

int search_discov_name_in_WDS_catalog(char *WDS_catalog, char *discov_name,
                                      char *WDS_name, int *found);
int get_data_from_WDS_catalog(char *WDS_catalog, char *discov_name,
                              char *comp_name,
                              char *wds_name, double *WdsLastYear,
                              double *WdsLastRho, double *WdsLastTheta,
                              double *WdsMagA, double *WdsMagB,
                              char *WdsSpectralType,
                              int *found);
int read_coordinates_from_WDS_catalog(char *WDS_name, char *WDS_catalog, 
                                      double *alpha, double *delta,
                                      double *equinox, int *found);
int get_data_from_WDS_and_HIP_catalogs(char *WDS_catalog, 
                                       char *HIC_catalog, char *HIP_catalog,
                                       char *discov_name, char *comp_name, 
                                       char *wds_name,
                                       double *magV, double *B_V_index,
                                       double *paral, double *err_paral,
                                       double *magV_A, double *magV_B,
                                       char *spectral_type, int *found_in_WDS,
                                       int *found_in_Hip_cat);

#endif   /* EOF sentry */

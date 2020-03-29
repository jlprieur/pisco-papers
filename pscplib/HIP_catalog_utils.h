/************************************************************************
* "HIP_catalog_utils.h"
* Initially to add the WDS/HIP/CCDM names to all objects 
* of the PISCO catalog 
*
* JLP 
* Version 12/02/2012
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>                  /* exit() */
#include <string.h>
#include <ctype.h>                   /* isprint... */
#include <math.h>
#include <time.h>                    /* date */

/*
#define DEBUG
#define DEBUG_1
*/

int HIP_name_from_HIP_HDS_WDS_cross(char *WDS_name, char *HIP_HDS_WDS_cross, 
                                    char *HIP_name, int *found);
int search_object_in_HIC_catalog(char *HIC_catalog, double alpha, double delta, 
                                 double equinox, char *HIP_name, 
                                 char *CCDM_name,
                                 double *V_mag, double *B_V_index, 
                                 double D_tolerance, int *found);
int read_data_in_HIP_catalog(char *HIP_catalog, char *HIP_name, 
                             double *paral, double *err_paral, int *found);
int check_consistency_coord_WDSname(double alpha_Pcat, double delta_Pcat, 
                                    double equinox_Pcat,
                                    char *object_name, char *WDS_name, 
                                    int *is_OK);
int check_consistency_ADSname(char *ADS_name, char *WDS_name, 
                              char *discov_name, char *ADS_WDS_cross, 
                              int *is_OK);

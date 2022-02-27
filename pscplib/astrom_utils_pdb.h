/*************************************************************************
* Set of routines to read Latex files with astrometric measurements
* created by Xdisp1 or Wdisp1 (e.g. astrom05a.tex) 
*
* JLP
* Version 19/09/2020
*************************************************************************/
#ifndef _astrom_utils_pdb_h  /* BOF sentry */
#define _astrom_utils_pdb_h  

#ifdef __cplusplus
extern "C" {
#endif

int astrom_add_WDS_from_discov(FILE *fp_in, FILE *fp_out, char *WDS_catalog); 
int astrom_add_discov_and_WDS(FILE *fp_in, FILE *fp_out, char *WDS_catalog, 
                              char *ADS_WDS_cross);
int astrom_decode_new_object_name(char *b_in, char *ads_name,
                                  char *discov_name, char *comp_name,
                                  char *WDS_name, double *year);
int get_discov_from_ads_wds_crossref(char *ADS_WDS_cross, char *ads_name1,
                                     char *discov_name1, char *comp_name1,
                                     char *wds_name1);
int get_ads_from_ads_wds_crossref(char *ADS_WDS_cross, char *discov_name1,
                                  char *ads_name1);


#ifdef __cplusplus
}
#endif


#endif /* EOF sentry */

/*************************************************************************
* Set of routines to read Latex files with astrometric measurements
* created by Xdisp1 or Wdisp1 (e.g. astrom05a.tex) 
*
* To convert a Latex array with measurements in pixels to array in arcseconds
* (and conversion of angles from XY axis to Noth-East referential)
*
* Format of input files:
\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}
& Name  & Epoch & Filter & Eyep. & $\rho$ & $\sigma_\rho$ & $\theta$ &
$\sigma_\theta$ & comments \\
 &       &       &        & (mm)
& (pixels) & (pixels) & (deg.) & (deg.) & \\
\hline % -----------------------------------------------------------------------
%%
16564+6502 = STF 2118 AB & ADS 10279 & 2004. & & & & & & & orb \\
or 16564+6502 = STF 2118 AB & 2004. & & & & & & & \\ (gili's format, old_Gdpisco_fmt=1)
%%
%% 090904_ads10279_Vd_a
%% Thresholds: -1,8
%% rho_center=14.62 (+/-1.87) theta_center=156.71 (+/-7.32)or theta=-23.29
%% rho_center=14.72 theta_center=156.70  or theta = -23.30
%% rho_center=14.37 (+/-1.89) theta_center=-23.64 (+/-7.54)
%% rho_center=14.34 theta_center=-23.73
%% mean: rho = 14.50 \pm 0.17  theta = -23.61 \pm 0.27 (too small -> 0.3)(n=8)
& 090904\_ads10279\_Vd\ & 09/09/2004 & V & 20 & 14.50 & 0.17 & -23.61 & 0.3 & \\
....
*
Keywords in the notes:
EP=2004.133 (BesselEpoch)
WR=1.2      (Rho, WDS)
WT=127      (Theta, WDS)
WY=2002     (Year of WDS observation) 
Q=2         (Quadrant with restricted triple-correplation)
LQ=3        (Quadrant with long integration)
*
* JLP
* Version 02/05/2010
*************************************************************************/
#include <stdio.h> 
#include <stdlib.h> // exit() 
#include <ctype.h>  // isprint 
#include <string.h> 
#include "astrom_utils1.h" 
#include "jlp_numeric.h" // JLP_QSORT, MINI, MAXI, PI, etc
#include "jlp_fitsio.h" // JLP_besselian_epoch 
#include "jlp_string.h"
#include "latex_utils.h"  // latex_read_svalue...

#define DEBUG
/*
*/

/*
int astrom_read_new_object_line_for_wds_or_ads(char *b_data, char *wds_name, 
                            char *discov_name, char *ads_name, 
                            double *WR, double *WT, double *WY, int i_notes,
                            int in_astrom_fmt);
int astrom_read_object_data_for_name_only(char *b_data, char *star_name,  
                            int *epoch_year, int in_astrom_fmt);
int astrom_check_measure(char *b_data, int i_eyepiece, int i_rho, 
                         int i_drho, int i_theta, int i_dtheta);
int astrom_calibrate_measures(OBJECT *obj, int nobj, double *calib_date1, 
                              int *calib_dd1, int *calib_mm1,
                              int *calib_year1, double *calib_scale1,
                              int *calib_eyepiece1, int *n_eyepieces1,
                              double *theta01, int *sign1, int ncalib1, int ndim);
int astrom_calib_data_copy(char *b_data, char *b_out,
                           double *calib_date1, int *calib_dd1, int *calib_mm1,
                           int *calib_year1, double *calib_scale1,
                           int *calib_eyepiece1, int *n_eyepieces1,
                           double *theta01, int *sign1, int ncalib1, int ndim,
                           int i_date, int i_eyepiece, int i_rho,
                           int i_drho, int i_theta, int i_dtheta);
int astrom_add_new_measure(char *b_data, OBJECT *obj, int i_obj, int i_filename, 
                           int i_date, int i_filter, int i_eyepiece, int i_rho,
                           int i_drho, int i_theta, int i_dtheta, int i_notes,
                           int comments_wanted);
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
int astrom_read_dmag_from_notes(char *notes, double *dmag, double *ddmag,
                                int comments_wanted);
int astrom_preformat_wds_name(char *obj_wds, char *wds_name);
int astrom_check_if_object_name(char *b_in, int *contains_object_name);
int astrom_ra_sort_objects(OBJECT *obj, int *index_obj, int nobj);
int astrom_name_sort_objects(OBJECT *obj, int *index_obj, int nobj);
int check_filter(char *filter, char *b_data);
*/


/*************************************************************************
*
* INPUT:
*  ncalib1 : number of calibrations
* OUTPUT:
*  scale2, theta02, sign2 (if status == 0)
*************************************************************************/
int get_calib_scale(double *calib_date1, int *calib_dd1, int *calib_mm1,
                    int *calib_year1, double *calib_scale1,
                    int *calib_eyepiece1, int *n_eyepieces1, double *theta01, 
                    int *sign1, int ncalib1, int ndim, double me_bessel_epoch, 
                    int me_eyepiece, double *scale2, double *theta02, 
                    int *sign2)
{
int i, k, icalib, eyepiece1, status = -1;
*scale2 = 0.;
*theta02 = 0.;
*sign2 = 1;

icalib = -1;
for(i = 0; i < ncalib1; i++) {
   if(me_bessel_epoch < calib_date1[i]) {
     icalib = i - 1;
     break;
     }
   } 

/* Scale according to eyepiece: */
 if(icalib >= 0) {
    *theta02 = theta01[icalib];
    *sign2 = sign1[icalib];
// printf("QQQ: me_bessel_epoch=%.4f icalib=%d calib_date1[icalib]=%.4f \n",
//         me_bessel_epoch, icalib, calib_date1[icalib]);
    for(k = 0; k < n_eyepieces1[icalib]; k++) {
      if(calib_eyepiece1[icalib * ndim + k] == me_eyepiece) { 
                   *scale2 = calib_scale1[icalib * ndim + k];
                   status = 0;
                   break;
                   }
   }
// printf("QQQ: scale=%.5f\n", *scale2);
 }
 if(status != 0) {
    fprintf(stderr, "Error: icalib=%d me_bessel_epoch=%.4f icalib=%d calib_date1[icalib]=%.4f \n",
         icalib, me_bessel_epoch, icalib, calib_date1[icalib]);
   }
return(status);
}
/*************************************************************************
*
 22388+4419 = HO 295 AB & ADS 16138 & 2004. & & & & & & & orb \\
*
* For 2004 data, WDS name was not present...
HO 295 AB & ADS 16138 & 2004. & & & & & & & orb \\
*
* INPUT:
* i_notes: column nber containing the notes
* in_astrom_fmt : =1 (gili_format) if empty ADS column, no orb info
* 22388+4419 = HO 295 AB &  & 2004. & & & & & & & \\
*                 =2 (calern_format) if no ADS column, no orb info
* 22388+4419 = HO 295 AB & 2004. & & & & & & & \\
*************************************************************************/
int astrom_read_new_object_line_for_wds_or_ads(char *b_data, 
                            char *wds_name, char *discov_name, char *comp_name,
                            char *ads_name, double *WR, double *WT, 
                            double *WY, int i_notes, int in_astrom_fmt)
{
char *pc, buff[64], wds_null[64], wds_not_found[64];
int istat, WDS_status, ADS_status, icol, ncol, new_object_line = 0;
int status;

WDS_status = 1;
ADS_status = 1;
ads_name[0] = '\0';
discov_name[0] = '\0';
comp_name[0] = '\0';
wds_name[0] = '\0';

strcpy(wds_null,"00000+0000");
strcpy(wds_not_found, "WDS_NAME_NOT_FOUND");

/* 2009/NEW/Preliminary step: check is syntax is OK: */
pc = b_data;
ncol = 1;
while(*pc) {
 if(*pc == '&') ncol++;
 pc++;
 }
// DEBUG 
printf("4555/in_astrom_fmt=%d ncol=%d\n", in_astrom_fmt, ncol);
/*
* in_astrom_fmt : =1 if gili's format (empty ADS column, no orb info)
*                 =2 if gili's format (no ADS column, no orb info)
* Gili's format:
* 14595+1753 = COU188 &  & 2013 & & & & & & & WY=2010 WT=227 WR=0.3 \\
* Calern format:
* 00226+5417 = A907   & 2017    & & & & & & & WY=2016 WT=211 WR=0.8 \\
*/
// & cou188\_a & 05/06/2013 &  & 1 & 3.50 & 0.10 & -43.21 & 0.6 & EP=2013.4292 Q=3 \\

// JLP 2021: now format has been changed and gili's format is ok with 10 columns
 new_object_line = 0;
 if(in_astrom_fmt == 1) {
  if((ncol == 10) && (b_data[0] != '&')) { 
     new_object_line = 1; 
   } else {
#ifdef DEBUG
  fprintf(stderr, "read_new_object_line/Bad new object syntax in following line:\n %s\n", b_data); 
  fprintf(stderr, "read_new_object_line/Error ncol=%d (all lines should have 10 columns!), in_astrom_fmt=%d\n", 
          ncol, in_astrom_fmt);
#endif
  return(-1); 
  }
 } else {
// JLP2022: I add the possibility of 8 columns
//  if((b_data[0] != '&') && ((ncol == 9) || (ncol == 8))) { 
  if((b_data[0] != '&') && (b_data[0] != '\\') && (b_data[0] != '%')) { 
printf("4555/OK1\n");
     new_object_line = 1; 
   } else {
printf("4555/NOK1\n");
#ifdef DEBUG
  fprintf(stderr, "read_new_object_line/Bad new object syntax in following line:\n %s\n", b_data); 
  fprintf(stderr, "read_new_object_line/Error ncol=%d (should be 9 for new object, 10 for measures !), in_astrom_fmt=%d\n", 
          ncol, in_astrom_fmt);
#endif
   return(-1); 
  }
 }

/* First step: decodes WDS name (not available for 2004a and 2004b...): */
icol = 1;
istat = latex_read_svalue(b_data, buff, icol); 

/* Look and check whether WDS name is present */
if(!istat && *buff) {
pc = buff;
buff[59]='\0';
/* Look for the second item (discov_name) starting after a "=" symbol: */
while(*pc && *pc != '=') pc++;
  if(*pc == '=') {
  pc++;
  strcpy(discov_name, pc);
  WDS_status = 0;
  }

/* Look for the first item ending with a "=" symbol: */
status = -1;
if(!WDS_status) {
  pc = buff;
  while(*pc && *pc != '=') pc++;
    if(*pc == '=') {
    *pc = '\0';
    strcpy(wds_name, buff);
    status = 0;
    }
} else if (strstr(b_data, wds_not_found) != NULL) {
// 16564+6502 = STF 2118 AB & ADS 10279 & 2004. & & & & & & & orb \\
   strcpy(wds_name, wds_null);
   strncpy(wds_name, wds_null, strlen(wds_null));
   strcpy(discov_name, buff);
#ifdef DEBUG
printf("wds_name=%s< len=%d wds_null=%s len2=%d\n", 
           wds_name, strlen(wds_name), wds_null, strlen(wds_null));
printf("wds_name=%s discov_name=%s\n", wds_name, discov_name);
#endif
   status = 0;
}

 ADS_status = 0;
/* Reads ADS name in 2nd column: */
if(in_astrom_fmt == 1) {
 ADS_status = 0;
 strcpy(ads_name, "");
/* Stop here if no discover's name was found either: */
 if(WDS_status) return(2);
} else {
icol = 2;
ads_name[0] = '\0';
istat = latex_read_svalue(b_data, buff, icol); 
if(!istat){
 ADS_status = 0;
 pc = buff;
 while(*pc && strncmp(pc,"ADS",3)) pc++;
 if(!strncmp(pc,"ADS",3)){
    pc +=3;
    strcpy(ads_name, pc);
    }
 else {
   strcpy(ads_name,"\\nodata");
  }
 }
}

/* Read WR, WT and WY values if present */
   latex_read_svalue(b_data, buff, i_notes);
   *WR = -1.; *WT = -1.; *WY = -1.;
   astrom_read_WDS_CHARA(buff, WR, WT, WY);

} /* EOF !status case */

return(status);
}
/*************************************************************************
*
 22388+4419 = HO 295 AB & ADS 16138 & 2004. & & & & & & & orb \\
or:
& A1053f200S & 2011 & & & & & & & \\
or:
HO209AB & 2015 & & & & & & & \\
*
* INPUT:
*  b_data
*  in_astrom_fmt: flag set to 1 if gili's format (no ADS column and no orb info)
*
* OUTPUT:
*  discov_name, comp_name, bessel_epoch_year
*************************************************************************/
int astrom_read_object_data_for_name_only(char *b_data, char *discov_name,
                                          char *comp_name, 
                                          double *bessel_epoch_year,
                                          int in_astrom_fmt)
{
int status = -1, ival;
char *pc, buff[60], full_name[64];
int istat, star_name_status, icol;

star_name_status = 1;
discov_name[0] = '\0';
comp_name[0] = '\0';
*bessel_epoch_year = 0.;

if(b_data[0] == '&') {
/* Read discov_name in 2nd column: */
icol = 2;
istat = latex_read_svalue(b_data, buff, icol); 
if((istat == 0) && (buff[0] != '\0')){
 star_name_status = 0;
 strcpy(full_name, buff);

/* Then try to read bessel_epoch_year in 3rd column: */
  star_name_status = -1;
  icol = 3;
  istat = latex_read_ivalue(b_data, &ival, icol); 
  if(istat == 0){
    if((ival > 2000) && (ival < 2100)) {
      *bessel_epoch_year = (double)ival;
      star_name_status = 0;
    }
  }
}
// EOF &
} else {
/* Read discov_name in 1st column: */
icol = 1;
istat = latex_read_svalue(b_data, buff, icol); 
if((istat == 0) && (buff[0] != '\0')){
 star_name_status = 0;
 strcpy(full_name, buff);
}
/* Then try to read bessel_epoch_year in 2nd column: */
  star_name_status = -1;
  icol = 2;
  istat = latex_read_ivalue(b_data, &ival, icol); 
  if(istat == 0){
    if((ival > 2000) && (ival < 2100)) {
      *bessel_epoch_year = (double)ival;
      star_name_status = 0;
    }
  }
}
// Split discov_name and comp_name:
if(star_name_status == 0) {
 jlp_split_discov_comp(full_name, 64, discov_name, comp_name);
 }

return(star_name_status);
}
/*************************************************************************
*
* Example: COU1142, A231AB, HSL1Aa,Ac
*************************************************************************/
int jlp_split_discov_comp(char *full_name, int fname_length,
                         char *discov_name, char *comp_name)
{
char *pc;

comp_name[0] = '\0';

/* Remove all the blanks: */
  jlp_compact_string(full_name, fname_length);

// copy to discov_name for output:
strcpy(discov_name, full_name);

// First scan the letters
pc = discov_name;
while(*pc) {
  if(isalpha(*pc) != 0) 
    pc++;
  else 
    break;
  }

// Then scan the numbers: 
while(*pc) {
  if(isdigit(*pc) != 0) 
    pc++;
  else 
    break;
  }

// copy to comp_name for output:
if(*pc) {
 strcpy(comp_name, pc);
 }

// Cut here the string for discov_name:
*pc = '\0';

/* DEBUG
printf("full=%s< disc=%s< comp=%s<\n", full_name, discov_name, comp_name);
*/

return(0);
}
/*************************************************************************
* Check whether input line is a measurement
*
* INPUT:
* i_eyepiece: column number of eyepiece focal length information
* i_rho: column number with rho values
* i_drho: column number with drho values
* i_theta: column number with theta values
* i_dtheta: column number with dtheta values
*
* return 0 if successful
*************************************************************************/
int astrom_check_measure(char *b_data, int i_eyepiece, int i_rho, 
                         int i_drho, int i_theta, int i_dtheta)
{
int status;
double ww; 
int iverbose = 0;

status = latex_read_dvalue(b_data, &ww, i_rho, iverbose); 
if(!status) status = latex_read_dvalue(b_data, &ww, i_drho, iverbose); 
if(!status) status = latex_read_dvalue(b_data, &ww, i_theta, iverbose); 
if(!status) status = latex_read_dvalue(b_data, &ww, i_dtheta, iverbose); 
if(!status) status = latex_read_dvalue(b_data, &ww, i_eyepiece, iverbose); 

return(status);
}
/*************************************************************************
*
* INPUT:
* calib_scale1 : scale values, scale_10, scale_20, scale_32
* calib_eyepiece1 : focal lengths, 10, 20, 32
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
* wds_name, discov_name, ads_name : object designation in various catalogues  
*
* OUPUT:
* obj.meas updated with calibrated values
*************************************************************************/
int astrom_calibrate_measures(OBJECT *obj, int nobj, double *calib_date1, 
                              int *calib_dd1, int *calib_mm1,
                              int *calib_year1, double *calib_scale1,
                              int *calib_eyepiece1, int *n_eyepieces1,
                              double *theta01, int *sign1, int ncalib1,
                              int ndim)
{
MEASURE *me;
double scale, theta02;
int nm, sign2, status; 
int i, j, k;

for(i = 0; i < nobj; i++) {
  nm = (obj[i]).nmeas;
  for(j = 0; j < nm; j++) {
    me = &(obj[i]).meas[j];
    status = get_calib_scale(calib_date1, calib_dd1, calib_mm1, calib_year1,
                             calib_scale1, calib_eyepiece1, n_eyepieces1,
                             theta01, sign1, ncalib1, ndim, me->bessel_epoch, 
                             me->eyepiece, &scale, &theta02, &sign2);
     if(status != 0) {
        fprintf(stderr, "astrom_calibrate_measures/Fatal error: obj[%d].discov=%s eyepiece=%d\n",
                i, (obj[i]).discov_name, me->eyepiece);
        exit(-1);
        break;
      } 

   if(me->rho != NO_DATA) { 
    me->rho *= scale;
    me->drho *= scale;
    me->theta = me->theta * sign2 + theta02;
    while(me->theta < 0.) me->theta += 360.;
    while(me->theta >= 360.) me->theta -= 360.;
    }
   } /* EOF loop on j */
} /* EOF loop on i */

return(0);
}
/*************************************************************************
*
* INPUT:
* b_data: line with raw measurement
* calib_scale1 : scale values, scale_10, scale_20, scale_32
* calib_eyepiece1 : focal lengths, 10, 20, 32
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
* wds_name, discov_name, ads_name : object designation in various catalogues  
*
* OUTPUT:
* b_out: line with calibrated measurement
*************************************************************************/
int astrom_calib_data_copy(char *b_data, char *b_out,
                           double *calib_date1, int *calib_dd1, int *calib_mm1,
                           int *calib_year1, double *calib_scale1,
                           int *calib_eyepiece1, int *n_eyepieces1,
                           double *theta01, int *sign1, int ncalib1, 
                           int ndim, int i_date, int i_eyepiece, int i_rho,
                           int i_drho, int i_theta, int i_dtheta)
{
int status, eyepiece, k, sign;
double BesselEpoch, rho, drho, theta, dtheta, scale, theta0;
char date[30];

/* By default simply copy input to output: */
 strcpy(b_out, b_data);

 status = astrom_decode_data(b_data, date, &BesselEpoch, &rho, &drho, &theta, 
                             &dtheta, &eyepiece, i_date, i_eyepiece, i_rho, 
                             i_drho, i_theta, i_dtheta);
// If not a good line with measurements, return from here
 if(status) {
  return(status);
  }

if(rho != NO_DATA) { 
  status = get_calib_scale(calib_date1, calib_dd1, calib_mm1, calib_year1,
                           calib_scale1, calib_eyepiece1, n_eyepieces1,
                           theta01, sign1, ncalib1, ndim, BesselEpoch, 
                           eyepiece, &scale, &theta0, &sign);
   if(status != 0) {
      fprintf(stderr, "astrom_calib_data_copy/get_calib_scale/Fatal error: BesselEpoch=%.3f eyepiece=%d\n",
              BesselEpoch, eyepiece);
      exit(-1);
   } 
 rho *= scale;
 drho *= scale;
 theta = theta * sign + theta0;
 while(theta < 0.) theta += 360.;
 while(theta >= 360.) theta -= 360.;
 dtheta = dtheta;

/* rho and drho with 3 decimals */
status = latex_write_dvalue(b_data, b_out, rho, i_rho, 3); 
if(!status) {
   strcpy(b_data, b_out);
   status = latex_write_dvalue(b_data, b_out, drho, i_drho, 3); 
   }
/* theta and dtheta with 1 decimal */
if(!status) {
   strcpy(b_data, b_out);
   status = latex_write_dvalue(b_data, b_out, theta, i_theta, 1); 
   }
if(!status) {
   strcpy(b_data, b_out);
   status = latex_write_dvalue(b_data, b_out, dtheta, i_dtheta, 1); 
   }
} /* EOF !NO_DATA */

/* BesselEpoch with 3 decimals (October 2008)*/
if((BesselEpoch > 0) && (!status)) {
   strcpy(b_data, b_out);
   status = latex_write_dvalue(b_data, b_out, BesselEpoch, i_date, 3); 
   }
if(status) {
 printf("calib_data/Fatal error, updating array!\n");
 return(-1);
 }

return(0);
}
/**************************************************************************
* Add a new measure for an existing object
*
* INPUT:
* i_obj: index of current object in OBJECT *obj
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
*
* OUTPUT:
* rho, drho, theta, dtheta
* eyepiece
* BesselEpoch
*
**************************************************************************/
int astrom_add_new_measure(char *b_data, OBJECT *obj, int i_obj, 
                           int i_filename, int i_date, int i_filter, 
                           int i_eyepiece, int i_rho, int i_drho, int i_theta, 
                           int i_dtheta, int i_notes, int comments_wanted)
{
MEASURE *me;
double BesselEpoch, JulianEpoch, rho, drho, theta, dtheta, ww, exact_epoch;
double dmag, ddmag;
char notes[64], filter[20], date[40], filename[60], nd_string[16], *pc;
int eyepiece, quadrant, dquadrant, status, nm;
int iverbose = 0, is_newdouble = 0, delete_epoch_item;

/* Return if no object has been entered yet: */
if(i_obj < 0) return(-1);

#ifdef DEBUG
printf("astrom_add_new_measure/add new measure from discov_name=%s (i=%d) nm=%d\n", 
       (obj[i_obj]).discov_name, i_obj, (obj[i_obj]).nmeas);
#endif

eyepiece = 0;
quadrant = 0;
dquadrant = 0;
rho = NO_DATA;
drho = NO_DATA;
theta = NO_DATA;
dtheta = NO_DATA;
BesselEpoch = -1.;
JulianEpoch = -1.;
dmag = -1;
ddmag = 0.;
filter[0] = '\0';
notes[0] = '\0';

/* Read date and compute approximative Bessel epoch from the date: */
   status = astrom_compute_bessel_epoch_value(b_data, date, &BesselEpoch, 
                                              i_date); 
#ifdef DEBUG
   printf("compute epoch with date: BesselEpoch=%.2f i_date=%d, status=%d\n", 
           BesselEpoch, i_date, status);
#endif

   status = latex_read_svalue(b_data, notes, i_notes);
   if(comments_wanted == 1) {
     delete_epoch_item = 0;
   } else {
     delete_epoch_item = 1;
   }
/* Check if EP= is available in the notes: */
   status = astrom_read_bessel_epoch_from_notes(notes, &exact_epoch, 
                                                delete_epoch_item);
   if(status == 0) BesselEpoch = exact_epoch;

/* Check if EJUL= is available in the notes: */
   status = astrom_read_julian_epoch_from_notes(notes, &JulianEpoch, 
                                                delete_epoch_item);
// Remove trailing and extra blanks:
   jlp_trim_string(notes, 64);

status = latex_read_dvalue(b_data, &rho, i_rho, iverbose); 
#ifdef DEBUG
printf("rho=%.2f i_rho=%d, status=%d\n", rho, i_rho, status);
#endif
status = latex_read_dvalue(b_data, &drho, i_drho, iverbose); 
status = latex_read_dvalue(b_data, &theta, i_theta, iverbose); 
status = latex_read_dvalue(b_data, &dtheta, i_dtheta, iverbose); 
status = latex_read_dvalue(b_data, &ww, i_eyepiece, iverbose); 
if(status == 0) eyepiece = (int)(ww+0.5);

/* Read non-compulsory parameters */
   latex_read_svalue(b_data, filename, i_filename);
   jlp_compact_string(filename, 60);
   latex_read_svalue(b_data, filter, i_filter);
/* JLP2008: I change sf to W filter: */
/* JLP 2008 I remove the first blank characters: */
   pc = filter; 
   while(*pc && *pc == ' ') pc++; 
   strcpy(filter, pc);
   pc = filter; 
/* Conversion to upper case characters (JLP2010) 
*  e.g., Rl becomes RL
*/
   while(*pc && isprint(*pc)) {*pc = toupper(*pc); pc++;} 
   *pc = '\0';
   if(!strncmp(filter,"SF",2)) {
/*
     printf("******** filter sf changed to >%s< \n",filter);
*/
     strcpy(filter, "W");
     } 
   status = check_filter(filter, b_data);

/* Check if Dm= is available in the notes: */
   astrom_read_dmag_from_notes(notes, &dmag, &ddmag, comments_wanted);

/* Check if ND is present in the notes: */
   strcpy(nd_string, "ND");
   if(strstr(notes, nd_string) != NULL) 
     is_newdouble = 1;
   else
     is_newdouble = 0;
//   printf(" notes=%s nd_string=%s new_double=%d \n", notes, nd_string, is_newdouble);

/* Read the quadrant value if present, and 
* removes "Q=" from notes if comments_wanted == 0 */
   astrom_read_quadrant_Q(notes, &quadrant, &dquadrant, comments_wanted);
/* 2008: attempt with LQ if problem with Q: */
/*
Q=2         (Quadrant with restricted triple-correplation)
LQ=3        (Quadrant with long integration)
*/
   if(quadrant == -1 || dquadrant == 1) 
      astrom_read_quadrant_LQ(notes, &quadrant, &dquadrant, comments_wanted);

/* Store data and increase nber of measurements for this object */
if(!status) {
   nm = (obj[i_obj]).nmeas;
   if(nm >= NMEAS) {
      fprintf(stderr, "astrom_new_measure/Fatal error: i_obj=%d nm=%d maximum number of measurements (=%d) is reached for %s\n", 
              i_obj, nm, NMEAS, filename); 
      exit(-1);
      }
   me = &(obj[i_obj]).meas[nm];

   me->rho = rho; 
/* Minimum value for rho error: 0.1 pixel or 0.5% */
   drho = MAXI(drho, 0.1); 
   me->drho = MAXI(drho, rho*0.005); 
   me->theta = theta; 
/* Minimum value for theta error: 0.3 degree */
   me->dtheta = MAXI(dtheta, 0.3); 
   me->quadrant = quadrant; 
   me->dquadrant = dquadrant; 
   me->dmag = dmag; 
   me->is_new_double = is_newdouble;
   me->flagged_out = 0;
   me->eyepiece = eyepiece; 
   me->bessel_epoch = BesselEpoch; 
   me->from_recorded_file = astrom_is_record_file(filename);
   strcpy(me->filename, filename);
   strcpy(me->filter, filter);
   strcpy(me->notes, notes);
   strcpy(me->date, date);
   ((obj[i_obj]).nmeas)++;

#ifdef DEBUG
     printf("astrom_add_new_measure ******************\n");
     printf(" i_obj=%d, new measure successfully added (nm=%d)\n", i_obj,
              (obj[i_obj]).nmeas);
     printf(" rho=%.2f drho=%.2f theta=%.2f dtheta=%.2f eyep.=%d Q=%d dQ=%d record=%d notes=>%s<\n", 
                   me->rho, me->drho, me->theta, me->dtheta, me->eyepiece,
                   me->quadrant, me->dquadrant, me->from_recorded_file, me->notes);
     printf("******************\n");
#endif
   }
#ifdef DEBUG
else
  printf("add_new_measure/Failure adding new measurement \n");
#endif

return(status);
}
/**************************************************************************
*
**************************************************************************/
int check_filter(char *filter, char *b_data)
{
int status = -1;

// Remove all blanks:
  jlp_compact_string(filter, 64);
  if(!strcmp(filter, "B") || !strcmp(filter, "V")
     || !strcmp(filter, "R") || !strcmp(filter, "RL")
     || !strcmp(filter, "I") || !strcmp(filter, "W") ) status = 0;

  if(status == -1) {
    fprintf(stderr, " Check_filter/Fatal error: filter=%s\n", filter);
    fprintf(stderr, " b_data=%s\n", b_data);
    exit(-1);
   }

return(status);
}
/**************************************************************************
* Read a new measure 
*
* INPUT:
* i_obj: index of current object in OBJECT *obj
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
*
* OUTPUT:
*   filename2, notes2
*   epoch2, filter2, eyepiece2
*   rho2, drho2, theta2, dtheta2
*
**************************************************************************/
int astrom_read_new_measure(char *b_data, int i_filename, int i_date, 
                            int i_filter, int i_eyepiece, int i_rho, int i_drho,                            int i_theta, int i_dtheta, int i_notes, 
                            char *filename2, char *notes2, 
                            double *BesselEpoch2, 
                            char *filter2, int *eyepiece2, double *rho2, 
                            double *drho2, double *theta2, double *dtheta2)
{
int comments_wanted2 = 1;
double ww, exact_epoch2, dmag2, ddmag2, JulianEpoch2;
char date2[40], *pc;
int quadrant2, dquadrant2, status, nm, iverbose = 0, delete_epoch_item;

*eyepiece2 = 0;
quadrant2 = 0;
dquadrant2 = 0;
*rho2 = NO_DATA;
*drho2 = NO_DATA;
*theta2 = NO_DATA;
*dtheta2 = NO_DATA;
*BesselEpoch2 = -1.;
JulianEpoch2 = -1.;
dmag2 = -1;
ddmag2 = 0.;
filter2[0] = '\0';
notes2[0] = '\0';

/* Read date and compute approximative Bessel epoch from the date: */
    status = astrom_compute_bessel_epoch_value(b_data, date2, BesselEpoch2, 
                                               i_date); 
#ifdef DEBUG
    printf("compute BesselEpoch with date: BesselEpoch=%.2f i_date=%d, status=%d\n", 
            *BesselEpoch2, i_date, status);
#endif
   if(comments_wanted2 == 1) {
     delete_epoch_item = 0;
   } else {
     delete_epoch_item = 1;
   }
/* Check if EP= is available in the notes: */
   status = astrom_read_bessel_epoch_from_notes(notes2, &exact_epoch2, 
                                                delete_epoch_item);
   if(status == 0) *BesselEpoch2 = exact_epoch2;
/* Check if EJUL= is available in the notes: */
   status = astrom_read_julian_epoch_from_notes(notes2, &JulianEpoch2, 
                                                delete_epoch_item);

status = latex_read_dvalue(b_data, rho2, i_rho, iverbose); 
#ifdef DEBUG
printf("rho=%.2f i_rho=%d, status=%d\n", *rho2, i_rho, status);
#endif
status = latex_read_dvalue(b_data, drho2, i_drho, iverbose); 
status = latex_read_dvalue(b_data, theta2, i_theta, iverbose); 
status = latex_read_dvalue(b_data, dtheta2, i_dtheta, iverbose); 
status = latex_read_dvalue(b_data, &ww, i_eyepiece, iverbose); 
if(status == 0) *eyepiece2 = (int)(ww+0.5);

/* Read non-compulsory parameters */
   latex_read_svalue(b_data, filename2, i_filename);
   latex_read_svalue(b_data, filter2, i_filter);
/* JLP2008: I change sf to W filter: */
/* JLP 2008 I remove the first blank characters: */
   pc = filter2; 
   while(*pc && *pc == ' ') pc++; 
   strcpy(filter2, pc);
/* Conversion to upper case characters (JLP2010) 
*  e.g., Rl becomes RL
*/
   pc = filter2; 
   while(*pc && isprint(*pc)) {*pc = toupper(*pc); pc++;} 
   *pc = '\0';
   if(!strncmp(filter2,"SF",2)) {
/*
     printf("******** filter sf changed to >%s< \n",filter2);
*/
     strcpy(filter2, "W");
     } 
   status = latex_read_svalue(b_data, notes2, i_notes);

/* Check if Dm= is available in the notes: */
   astrom_read_dmag_from_notes(notes2, &dmag2, &ddmag2, comments_wanted2);

/* Read the quadrant value if present, and 
* removes "Q=" from notes if comments_wanted == 0 */
   astrom_read_quadrant_Q(notes2, &quadrant2, &dquadrant2, comments_wanted2);
/* 2008: attempt with LQ if problem with Q: */
/*
Q=2         (Quadrant with restricted triple-correplation)
LQ=3        (Quadrant with long integration)
*/
   if(quadrant2 == -1 || dquadrant2 == 1) 
      astrom_read_quadrant_LQ(notes2, &quadrant2, &dquadrant2, comments_wanted2);

/* Minimum value for rho error: 0.1 pixel or 0.5% */
   *drho2 = MAXI(*drho2, 0.1); 
   *drho2 = MAXI(*drho2, (*rho2)*0.005); 
/* Minimum value for theta error: 0.3 degree */
   *dtheta2 = MAXI(*dtheta2, 0.3); 

return(status);
}
/**************************************************************************
* Decode line with results of measurements
* Example:
* & 041107\_ads15547ab\_Rd\_8\_a & 04/11/2007 & R & 20 & 11.55 & 0.37 & 49.12 
* & 0.7 & Q=4 EP=2007.8437 \\
*
* INPUT:
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
*
* OUTPUT:
* rho, drho, theta, dtheta
* eyepiece
* BesselEpoch
*
**************************************************************************/
int astrom_decode_data(char *b_data, char *date, double *BesselEpoch, 
                       double *rho, 
                       double *drho, double *theta, double *dtheta, 
                       int *eyepiece, int i_date, int i_eyepiece, int i_rho, 
                       int i_drho, int i_theta, int i_dtheta)
{
double ww;
int status, ncol, iverbose = 0;
char *pc;

/* 2009/NEW/Preliminary step: check is syntax is OK: */
pc = b_data;
ncol = 1;
while(*pc) {
 if(*pc == '&') ncol++;
 pc++;
 }
if(ncol != 10) {
  fprintf(stderr, "astrom_decode_data/Bad syntax in the following line:\n %s\n", b_data); 
  fprintf(stderr, "astrom_decode_data/Error ncol=%d (should be 10!)\n", 
          ncol);
  fprintf(stderr, "i_date=%d i_eyepiece=%d i_rho=%d i_drho=%d i_theta=%d i_dtheta=%d\n", i_date, i_eyepiece, i_rho, i_drho, i_theta, i_dtheta); 
  fprintf(stderr, "line: %s\n", b_data); 
  return(-1); 
  }

*eyepiece = 0;
*rho = *drho = *theta = *dtheta = 0;

/* Read date: */
*BesselEpoch = 0.;
status = astrom_compute_bessel_epoch_value(b_data, date, BesselEpoch, i_date); 

status = latex_read_dvalue(b_data, rho, i_rho, iverbose); 
if(!status) status = latex_read_dvalue(b_data, drho, i_drho, iverbose); 
if(!status) status = latex_read_dvalue(b_data, theta, i_theta, iverbose); 
if(!status) status = latex_read_dvalue(b_data, dtheta, i_dtheta, iverbose); 
if(!status) {
  status = latex_read_dvalue(b_data, &ww, i_eyepiece, iverbose); 
  *eyepiece = (int)(ww+0.5);
  }

#ifdef DEBUG
if(!status) printf(" rho=%.2f drho=%.2f theta=%.2f dtheta=%.2f eyepiece=%d\n", 
                   *rho, *drho, *theta, *dtheta, *eyepiece);
#endif

return(status);
}
/**************************************************************************
* Read BesselEpoch value from date information in column #icol from b_data string
*
* Input format: dd/mm/yy, e.g. 12/2/2004 or 01/12/1998 or 31/06/2002
* Output: BesselEpoch, as a fraction of year, e.g. 2004.234
*
**************************************************************************/
int astrom_compute_bessel_epoch_value(char *b_data, char *date, 
                                      double *BesselEpoch, int icol) 
{
int ival, status, dd, mm, iyy;
double yy, time;

/* Read date: */
date[0] = '\0';
status = latex_read_svalue(b_data, date, icol); 

/* Removes the blanks in "date" string: */
sscanf(date, "%s", date);

if(!status) { 
   ival = sscanf(date, "%d/%d/%d", &dd, &mm, &iyy);
/*
printf("astrom_compute_bessel_epoch_value/date=>%s< dd=%d mm=%d iyy=%d ival=%d\n", 
        date, dd, mm, iyy, ival);
*/
   if(ival < 3) status = 1;
  }

if(!status) { 
/* Assume observations at 10:30 pm, local time, i.e., 20:30 (U.T) in summer */
/* Assume observations at 9:30 pm, local time, i.e., 20:30 (U.T) in winter */
time = 20.5;
yy = (double)iyy;
/* In "jlp_fits_utils.c": */
status = JLP_besselian_epoch(yy, mm, dd, time, BesselEpoch);

/*
printf("astrom_compute_bessel_epoch_value/ BesselEpoch=%f\n", *BesselEpoch); 
*/
}

return(status);
}
/***************************************************************************
* astrom_name_sort_objects
* Calling routine  JLP_QSORT_INDX_CHAR()
*
* INPUT:
*  array[nn]: array to be sorted
*
* OUTPUT:
*  array[nn]: sorted array
*  index[nn]: array giving the index of the input array,
*             to sort other arrays in the same way if necessary
*             (f.i. array2[i] := array2[index[i]])
****************************************************************************/
int astrom_name_sort_objects(OBJECT *obj, int *index_obj, int nobj)
{
int k, length, j, status;
char buffer[60*10000];

if(nobj > 10000) 
   {
   printf("astrom_name_sort_objects/Fatal error: nobj=%d > %d\n",
           nobj, 10000);
   exit(-1);
   }
/* When length = 1, sort elements according to first letter only */
length = 60;

for(j = 0; j < nobj; j++) 
  sprintf(&buffer[j*60], "%s%s", obj[j].discov_name, obj[j].comp_name);

/* Sort "buffer" array by alphabetic order: */
JLP_QSORT_INDX_CHAR(buffer, &length, index_obj, &nobj);

return(0);
}
/***************************************************************************
* astrom_ra_sort_objects
* Calling routine  int JLP_QSORT_INDX_DBLE(double *array, int *index, int *nn)
* INPUT:
*  array[nn]: array to be sorted
*
* OUTPUT:
*  array[nn]: sorted array
*  index[nn]: array giving the index of the input array,
*             to sort other arrays in the same way if necessary
*             (f.i. array2[i] := array2[index[i]])
****************************************************************************/
int astrom_ra_sort_objects(OBJECT *obj, int *index_obj, int nobj)
{
double *ra;
int i, j1, j2, nswap;

ra = (double *)malloc((nobj) * sizeof(double));

for(i = 0; i < nobj; i++) ra[i] = obj[i].ra;

JLP_QSORT_INDX_DBLE(ra, index_obj, &nobj);

/* Final check on declination criterion too: */
nswap = 1;
while(nswap > 0) {
nswap = 0;
for(i = 0; i < nobj-1; i++) {
 j1 = index_obj[i];
 j2 = index_obj[i+1];
/* Swap indices if declination is not sorted out: */
 if((obj[j1].ra == obj[j2].ra) && (obj[j1].dec > obj[j2].dec)) {
      index_obj[i] = j2;
      index_obj[i+1] = j1;
      nswap++;
#ifdef DEBUG
      printf("ra1=%f ra2=%f dec1=%d dec2=%d\n", obj[j1].ra, obj[j2].ra, 
              obj[j1].dec, obj[j2].dec);
#endif
   } 
 } /* EOF loop on i */
#ifdef DEBUG
printf("Sorting declination now: nswap=%d for this iteration\n",nswap);
#endif
} /* EOF while loop */

free(ra);
return(0);
}
/***************************************************************************
* Read WDS_CHARA data from keywords (WR, WT, WY) if present in the note column
*
* INPUT:
*  notes: string from notes column of input file 
*
* OUTPUT
* WR, WT, WY  (-1 if absent)
*
***************************************************************************/
int astrom_read_WDS_CHARA(char *notes, double *WR, double *WT, double * WY)
{
char *pc;
int status;
status = 0;

*WR=-1.;
*WT=-1;
*WY=-1;

pc = notes;
while(*pc && strncmp(pc,"WR=",3)) pc++;

if(!strncmp(pc,"WR=",2)){
 pc +=3;
 sscanf(pc,"%lf",WR);
   }

pc = notes;
while(*pc && strncmp(pc,"WY=",3)) pc++;

if(!strncmp(pc,"WY=",2)){
 pc +=3;
 sscanf(pc,"%lf",WY);
   }

pc = notes;
while(*pc && strncmp(pc,"WT=",3)) pc++;

if(!strncmp(pc,"WT=",2)){
 pc +=3;
 sscanf(pc,"%lf",WT);
   }

return(status);
}
/***************************************************************************
* Read the quadrant value if present, 
* and removes "Q=*" from notes (if comments_wanted == 0) 
*
* INPUT:
*  notes: string from notes column of input file 
*
* OUTPUT
* quadrant: 0 if Q was not present
*           -1 if Q=?
*           1, 2, 3 or 4 if Q=1, Q=2, Q=3 or Q=4
*  notes: same as input string, but without "Q=*" (if comments_wanted == 0) 
*
***************************************************************************/
int astrom_read_quadrant_Q(char *notes, int *quadrant, int *dquadrant,
                           int comments_wanted)
{
char *pc, *pc1, buff1[256], buff2[256];
int status, k;
status = 0;

// Default values:
*quadrant=-1;
*dquadrant=1;

pc = strstr(notes, "Q=");

/* Set quadrant indetermination to 1 
* if is a question marks such as Q=? Q=1?, Q=2?, Q=3? or Q=4? */

if(pc != NULL) {
 strcpy(buff1, pc);
 pc1 = buff1;
 while(*pc1 && (isalpha(*pc1) || *pc1 != ' ' || isdigit(*pc1))) pc1++; 
 *pc1 = '\0';
 if(!strcmp(buff1, "Q=1")) {
    *quadrant = 1;
    *dquadrant = 0;
 } else if(!strcmp(buff1, "Q=2")) {
    *quadrant = 2;
    *dquadrant = 0;
 } else if(!strcmp(buff1, "Q=3")) {
    *quadrant = 3;
    *dquadrant = 0;
 } else if(!strcmp(buff1, "Q=4")) {
    *quadrant = 4;
    *dquadrant = 0;
 }

#ifdef DEBUG
printf("DEBUG before notes=%s< buff1=%s<  quad=%d dquad=%d ", notes,
                          buff1, *quadrant, *dquadrant);
#endif

/* Update notes, after removing quadrant indication: */
if(comments_wanted == 0) {
// buff2 is the second part after Q=2?
  pc = strstr(notes, "Q=");
  pc +=2;
  while(*pc && ((*pc == '-') || (*pc == '?') || isdigit(*pc))) pc++; 
  if(*pc) strcpy(buff2, pc);
  else strcpy(buff2, "");
// buff1 is the first part after Q=2?
  strcpy(buff1, notes);
  pc = strstr(buff1, "Q=");
  *pc = '\0'; 
  sprintf(notes, "%s %s", buff1, buff2);
#ifdef DEBUG
printf("after notes=%s<", notes);
#endif
  }
#ifdef DEBUG
printf("\n");
#endif

} // EOF pc != NULL

return(status);
}
/*********************************************************************
* Nearly the same as read_quadrant_Q but with LQ keyword
* No longer used ?
*
* Example:
* Q=2         (Quadrant with restricted triple-correplation)
* LQ=3        (Quadrant with long integration)
*
*********************************************************************/
int astrom_read_quadrant_LQ(char *notes, int *quadrant, int *dquadrant,
                            int comments_wanted)
{
char *pc, buffer[80];
int status, k;
status = 0;

/* Compared to read_quadrant_Q,
* read_quadrant_LQ should not initialize quadrant and dquadrant! */

pc = notes;
k = 0;
for(k = 0; k < 68; k++) 
  if(strncmp(pc,"LQ=",3)) {
    buffer[k] = *pc; pc++;
  } else {
    break;
  }

if(!strncmp(pc,"LQ=",3)){
 pc +=3;
 switch(*pc) {
   case '1':
    *quadrant=1;
    break;
   case '2':
    *quadrant=2;
    break;
   case '3':
    *quadrant=3;
    break;
   case '4':
    *quadrant=4;
    break;
   default:
   case '?':
    *quadrant=-1;
    break;
   }

/* Set quadrant indetermination by checking whether there
* is a question marks such as Q=1?, Q=2?, Q=3? or Q=4? */
 if(*pc) pc++;
 if(*pc == '?') {pc++; *dquadrant = 1;}
 else *dquadrant = 0;

/* Normal syntax is f.i. "LQ=4," */
 pc++;
 while(*pc && k < 79) {pc++; buffer[k++] = *pc;}
 buffer[k++] = '\0';

/* Update notes, after removing quadrant indication: */
if(!comments_wanted) strcpy(notes,buffer);
}

return(status);
}
/*******************************************************************************
* Check if theta value is compatible with the quadrant value
*
*******************************************************************************/
int astrom_quadrant_is_consistent(MEASURE *me)
{
/* tol = tolerance for accepting the quadrant (5 degrees here)
*/
double theta, tol = 5.;
int quad, good_quad;

quad = me->quadrant;
theta = me->theta;

switch(quad)
  {
   case 1:
     if(theta > (0. - tol) && theta < (90. + tol)) good_quad = 1;
     else good_quad = -1;
     break;
   case 2:
     if(theta > (90. - tol) && theta < (180. + tol)) good_quad = 1;
     else good_quad = -1;
     break;
   case 3:
     if(theta > (180. - tol) && theta < (270. + tol)) good_quad = 1;
     else good_quad = -1;
     break;
   case 4:
     if(theta > (270. - tol) && ((theta < (360. + tol)) || (theta < tol))) 
           good_quad = 1;
     else good_quad = -1;
     break;
/* Case Q=? */
   case -1:
   default:
     good_quad = 0;
     break;
  }

return(good_quad);
}
/*********************************************************************
*
*********************************************************************/
int astrom_compute_statistics(FILE *fp_out, OBJECT *obj, int nobj, char *filein)
{
MEASURE *me;
double delta_theta, delta_rho;
char *pc, latex_filein[128];
int nm, nmeas, no_detected, no_data;
int nquad, nquad_uncert, is_quad, nbad_quad, k;
int i, j;

/* Compute number of measurements: */
nmeas = 0;
no_detected = 0;
nquad = 0;
nquad_uncert = 0;
nbad_quad = 0;

for(i = 0; i < nobj; i++) {
 nm = (obj[i]).nmeas;
   for(j = 0; j < nm; j++) {
   me = &(obj[i]).meas[j];
     if(!me->flagged_out) {
     nmeas ++;
     no_data = ((me->rho == NO_DATA) || (me->theta == NO_DATA)) ? 1 : 0; 
     no_detected += no_data;
     is_quad = (me->quadrant > 0) ? 1 : 0;
     nquad += is_quad;
     nquad_uncert += me->dquadrant;
// Compare quadrants if WDS measures more recent than 1980:
     if((me->dquadrant != 1) && (me->theta != NO_DATA) && (obj[i]).WY > 1980. && is_quad) {
        delta_theta = me->theta - (obj[i]).WT;
        if(delta_theta > 180.) delta_theta -= 360.;
        if(delta_theta < -180.) delta_theta += 360.;
        delta_theta = ABS(delta_theta);
/* Check if difference is not too large: */
        if(delta_theta > 90.) {
          nbad_quad ++;
          fprintf(fp_out,"Inconsistent quadrant for %s/%s%s %s (Q=%d, mtheta=%.1f WT=%.1f Dt=%.1f WY=%d)\n",
                (obj[i]).wds, (obj[i]).discov_name, (obj[i]).comp_name, 
                me->filename, me->quadrant, me->theta, 
                (obj[i]).WT, delta_theta, (int)(obj[i]).WY);
          }
       }
// New JLP2018, also compare the rho measurements:
// Compare rho if WDS measures more recent than 1980:
      if((me->rho != NO_DATA) && ((obj[i]).WR > 0. && me->rho > 0.) && ((obj[i]).WY > 1980.)) {
        delta_rho = me->rho - (obj[i]).WR;
        delta_rho = ABS(delta_rho) / me->rho;
        if(delta_rho > 0.5) {
        fprintf(fp_out,"Inconsistent rho for %s/%s%s (mrho=%.2f WR=%.2f Dr/r=%.2f WY=%d)\n",
                (obj[i]).wds, (obj[i]).discov_name, (obj[i]).comp_name, 
                me->rho, (obj[i]).WR, delta_rho, (int)(obj[i]).WY);
         }
        }
     } // !me->flagged_out
   }
 }

pc = filein;
k = 0;
while(*pc) {
/* Change "_" to "\_" for LateX : */
 if(*pc == '_') {latex_filein[k++] = '\\'; latex_filein[k++] = '_';}
 else {latex_filein[k++] = *pc;}
 pc++;
 }

fprintf(fp_out, " Input file: %s \n \n", latex_filein);
printf(" Number of objects: %d \n", nobj);
fprintf(fp_out, " Number of objects: %d \n \n", nobj);

printf(" Number of observations: %d with %d measurements and %d cases of no detection \n ", 
        nmeas, nmeas - no_detected, no_detected);
fprintf(fp_out, " Number of observations: %d with %d measurements and %d cases of no detection \n \n ", 
        nmeas, nmeas - no_detected, no_detected);

printf("Quadrant was determined for %d measurements (rejecting %d uncertain determinations)\n", 
        nquad, nquad_uncert);
fprintf(fp_out, "Quadrant was determined for %d measurements (rejecting %d uncertain determinations)\n \n", 
        nquad, nquad_uncert);

printf("Warning: %d quadrant values are inconsistent with CHARA theta last measurements! \n",
         nbad_quad);
fprintf(fp_out,"Warning: %d quadrant values are inconsistent with CHARA theta last measurements! \n \n",
         nbad_quad);
printf("OK: %d quadrant values are consistent with CHARA theta last measurements! \n",
         nquad - nbad_quad);
fprintf(fp_out,"OK: %d quadrant values are consistent with CHARA theta last measurements! \n \n",
         nquad - nbad_quad);

fprintf(fp_out, " In column 9,  $*$ indicates that $\\theta$ was determined with our quadrant value (or with the long integration)\n \n");
fprintf(fp_out, " In column 9,  $!$ indicates that $\\theta$ could not be determined neither with this value, nor with WDS CHARA last measurement\n \n");

return(0);
}
/******************************************************************************
*
* Good value if WDS_CHARA measurement is in the interval 
*          [me->theta - 90, me->theta + 90]
* Otherwise, correct our measurement by adding +180 degrees.
******************************************************************************/
int astrom_correct_theta_with_WDS_CHARA(OBJECT *ob, MEASURE *me)
{
double delta_theta;
int status;
status = 0;

if(ob->WY == -1. || ob->WT == -1. || me->theta == NO_DATA) {
   status = -1;
  }
else {
  delta_theta = ABS(me->theta - ob->WT);
  if(delta_theta > 90.) {
    me->theta += 180.;
    if(me->theta > 360.) me->theta -= 360.;
    }
  }
return(status);
}
/***************************************************************************
* Read the Besselian epoch value if present, 
* and removes "EP=*" from notes (if delete_epoch_item != 0) 
*
* INPUT:
*  notes: string from notes column of input file 
*
* OUTPUT
*  BesselEpoch: -1 if EP was not present
*  notes: same as input string, but without "EJUL=*" (if delete_epoch_item != 0)
*
***************************************************************************/
int astrom_read_bessel_epoch_from_notes(char *notes, double *BesselEpoch, 
                                        int delete_epoch_item)
{
char *pc, buffer[80];
int status, k;
status = -1;

*BesselEpoch=-1.;

pc = notes;
k = 0;
while(*pc && strncmp(pc,"EP=",3)) {buffer[k++] = *pc; pc++;}

/* Syntax is EP=2005.2341 */
if(!strncmp(pc,"EP=",3)){
 pc +=3;
 sscanf(pc,"%lf", BesselEpoch);
 status = 0;

/* Normal syntax is f.i. "EP=2005.2341 " */

/* Look for comma, blank or end of notes: */
 while(*pc && *pc != ',' && *pc != ' ') pc++;

/* Copy end of notes: */
 while(*pc && k < 79) {pc++; buffer[k++] = *pc;}
 buffer[k++] = '\0';

/* Update notes, after removing epoch data: */
if(delete_epoch_item != 0) strcpy(notes,buffer);
}

return(status);
}
/***************************************************************************
* Read the Julian epoch value if present, 
* and removes "EJUL=*" from notes (if delete_epoch_item != 0) 
*
* INPUT:
*  notes: string from notes column of input file 
*
* OUTPUT
*  JulianEpoch: -1 if EJUL was not present
*  notes: same as input string, but without "EJUL=*" (if delete_epoch_item != 0)
*
***************************************************************************/
int astrom_read_julian_epoch_from_notes(char *notes, double *JulianEpoch, 
                                        int delete_epoch_item)
{
char *pc, buffer[80];
int status, k;
status = -1;

*JulianEpoch=-1.;

pc = notes;
k = 0;
while(*pc && strncmp(pc,"EJUL=",5)) {buffer[k++] = *pc; pc++;}

/* Syntax is EP=2005.2341 */
if(!strncmp(pc,"EJUL=",5)){
 pc +=5;
 sscanf(pc,"%lf", JulianEpoch);
 status = 0;

/* Normal syntax is f.i. "EJUL=2005.2341 " */

/* Look for comma, blank or end of notes: */
 while(*pc && *pc != ',' && *pc != ' ') pc++;

/* Copy end of notes: */
 while(*pc && k < 79) {pc++; buffer[k++] = *pc;}
 buffer[k++] = '\0';

/* Update notes, after removing epoch data: */
if(delete_epoch_item != 0) strcpy(notes,buffer);
}

return(status);
}
/***************************************************************************
* Read Dmag if present, 
* and removes "Dm=*" from notes (if comments_wanted == 0) 
*
* INPUT:
*  notes: string from notes column of input file 
*
* OUTPUT
*  dmag: -1 if Dm was not present
*  ddmag: -1 if Dm was not present
*  notes: same as input string, but without "Dm=*" (if comments_wanted == 0)
*
***************************************************************************/
int astrom_read_dmag_from_notes(char *notes, double *dmag, double *ddmag,
                                int comments_wanted)
{
char *pc, buffer[80];
int status, k;
status = 0;

*dmag=-1.;
*ddmag=-1.;

// First go with "dm"
pc = notes;
k = 0;
while(*pc && strncmp(pc,"dm=",2)) {buffer[k++] = *pc; pc++;}

/* Syntax is dm1.23 */
if(!strncmp(pc,"dm",2)){
 pc +=2;
 sscanf(pc,"%lf", dmag);

/* Normal syntax is f.i. "dm1.23 " */

/* Look for question mark, blank or end of notes: */
 while(*pc && *pc != '?' && *pc != ' ') pc++;

/* Copy end of notes: */
 while(*pc && k < 79) {pc++; buffer[k++] = *pc;}
 buffer[k++] = '\0';

/* Update notes, after removing dmag data: */
if(!comments_wanted) strcpy(notes,buffer);
}

// Second  go with "Dm="
pc = notes;
k = 0;
while(*pc && strncmp(pc,"Dm=",3)) {buffer[k++] = *pc; pc++;}

/* Syntax is Dm=1.23 */
if(!strncmp(pc,"Dm=",3)){
 pc +=3;
 sscanf(pc,"%lf", dmag);

/* Normal syntax is f.i. "Dm=1.23 " */

/* Look for comma, blank or end of notes: */
 while(*pc && *pc != ',' && *pc != ' ') pc++;

/* Copy end of notes: */
 while(*pc && k < 79) {pc++; buffer[k++] = *pc;}
 buffer[k++] = '\0';

/* Update notes, after removing dmag data: */
if(!comments_wanted) strcpy(notes,buffer);
}

return(status);
}
/*************************************************************************
* Conversion of WDS name for Latex output, so that
* WDS names with negative declination will be aligned with positive ones 
*
* INPUT:
* obj_wds: e.g., 17304-0104
*
* OUTPUT:
* wds_name: e.g., 17304$-$0104
*
*************************************************************************/
int astrom_preformat_wds_name(char *obj_wds, char *wds_name)
{
char *pc;
int k;

pc = obj_wds;
k = 0;
while(*pc && *pc != '-' && k < 40) wds_name[k++] = *pc++;
if(*pc == '-') {
  strcpy(&wds_name[k], "$-$");
  k += 3;
  pc++;
  while(*pc) wds_name[k++] = *pc++;
 }
wds_name[k] = '\0';

#ifdef DEBUGG
printf("preformat_wds_name/input=%s output=%s\n", obj_wds, wds_name);
#endif

return(0);
}
/*********************************************************************
* Check if current line is compatible with the syntax of a new object
* Example:
* 16564+6502 = STF 2118 AB & ADS 10279 & 2004. & & & & & & & orb \\
* Or:
*  & ADS 10279 & 2004. & & & & & & & \\
* Or:
  COU 2118 &  & 2004. & & & & & & & orb \\
*
* OUTPUT:
* contains_object_name: flag set to one if ADS or discoverer's name was found
* contains_WDS_name: flag set to one if WDS name was found
*********************************************************************/
int astrom_check_if_object_name(char *b_in, int *contains_object_name,
                                int *contains_WDS_name)
{
int ncol, icol, istat, iverbose = 0;
double year;
char *pc;

*contains_object_name = 0;
*contains_WDS_name = 0;

/* If comment, exit from here: */
if(b_in[0] == '%' || b_in[1] == '%' || b_in[0] == '\\') return(0);

/* First check that there are 10 columns: */
pc = b_in;
ncol = 1;
while(*pc) {
  if(*pc == '&') ncol++;
  pc++;
 }
if(ncol != 10) return(-2);

/* Then check that there is a "+" or "-" (WDS name)
* or  a "=" in the first column: */
pc = b_in;
ncol = 1;
while(*pc) {
  if(isalpha(*pc)) *contains_object_name = 1;
  if(*pc == '+' || *pc == '-' || *pc == '=') {
    *contains_WDS_name = 1;
    *contains_object_name = 1;
    return(0);
    }
  if(*pc == '&') {ncol++; break;}
  pc++;
 }
 
/* Then check that there is a date in the 3rd column: */
icol = 3;
year = 0.;
istat = latex_read_dvalue(b_in, &year, icol, iverbose); 
if(!istat && year > 2000. && year < 2020.){
  *contains_object_name = 1;
  }
 
return(0);
}
/************************************************************************
* Check if recorded file,
* i.e. filename with _Br _Rr _Vr or _sfr
*
************************************************************************/
int astrom_is_record_file(char *filename)
{
char *pc;
int is_record_file;

is_record_file = 0;

/* Check if recorded file: */
pc = filename;
while(*pc && strncmp(pc,"_Vr",3) && strncmp(pc,"_Rr",3)
      && strncmp(pc,"_sfr",4) && strncmp(pc,"_Br",3)) pc++;
if(!strncmp(pc,"_Vr",3) || !strncmp(pc,"_Rr",3)
      || !strncmp(pc,"_sfr",4) || !strncmp(pc,"_Br",3)) is_record_file = 1;

return(is_record_file);
}
/************************************************************************
* Check if direct file,
* i.e. filename with _Bd _Rd _Vd or _sfd
*
************************************************************************/
int astrom_is_direct_file(char *filename)
{
char *pc;
int is_direct_file;

is_direct_file = 0;

/* Check if direct file: */
pc = filename;
while(*pc && strncmp(pc,"_Vd",3) && strncmp(pc,"_Rd",3)
      && strncmp(pc,"_Rld",4) && strncmp(pc,"_RLd",4)
      && strncmp(pc,"_sfd",4) && strncmp(pc,"_Bd",3)) pc++;
if(!strncmp(pc,"_Vd",3) || !strncmp(pc,"_Rd",3)
      || !strncmp(pc,"_Rld",4) || !strncmp(pc,"_RLd",4)
      || !strncmp(pc,"_sfd",4) || !strncmp(pc,"_Bd",3)) is_direct_file = 1;

return(is_direct_file);
}

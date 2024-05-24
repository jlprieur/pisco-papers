/************************************************************************
* "HIP_catalog_utils.c"
* Initially to add the WDS/HIP/CCDM names to all objects 
* of the PISCO catalog 
*
* JLP 
* Version 20/04/2009
*************************************************************************/
#include "HIP_catalog_utils.h" 
#include "jlp_string.h"

#define ABS(a) ((a) < 0.0  ? (-(a)) : (a))
#ifndef PI
#define PI 3.14159265
#endif
#define DEGTORAD   (PI/180.00)


/*
#define DEBUG
#define DEBUG_1
*/

/* Contained here:

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
*/
/*****************************************************************************
* Check the consistency between the WDS name and the coordinates
*
* INPUT:
* alpha_Pcat, delta_Pcat, equinox_Pcat: coordinates of the object 
*                                       from the PISCO catalog
*             (alpha in hours and delta in degrees)
* 
* object_name: name in the PISCO catalog
* WDS_name: WDS name found in the WDS catalog
*
* OUTPUT:
* is_OK: 1 if consistent, 0 otherwise
*****************************************************************************/
int check_consistency_coord_WDSname(double alpha_Pcat, double delta_Pcat, 
                                    double equinox_Pcat, char *object_name, 
                                    char *WDS_name, int *is_OK)
{
double alpha_WDS, delta_WDS, D_alpha, D_delta;
int a1, a2, d1, d2, nval;
char sign[1];

*is_OK = 0;

/* Decode the WDS name
Example:
23244+6917
corresponds to:
23 H 24.4 min +69 deg 17' 
*/
nval = sscanf(WDS_name, "%2d%3d%c%2d%2d", &a1, &a2, sign, &d1, &d2);
if(nval != 5) {
 fprintf(stderr, 
         "check_consistency_coord_WDSname/Fatal error reading WDS name: >%s<\n",
         WDS_name);
 exit(-1);
 }
alpha_WDS = (double)a1 + ((double)a2)/600.;
delta_WDS = (double)d1 + ((double)d2)/60.;

if(sign[0] == '-') delta_WDS *= -1.;

/* Difference of alpha in hours and of delta in degrees */
D_alpha = ABS(alpha_WDS - alpha_Pcat);
D_delta = ABS(delta_WDS - delta_Pcat);

/* Tolerance: 0.1 degrees */
if(D_alpha > 0.1/15. || D_delta > 0.1) {
  *is_OK = 0.;
   fprintf(stderr, "check_consistency/Error: %s WDS name %s (%.3f %.3f) is not consistent with the coordinates of your catalog\n", 
         object_name, WDS_name, alpha_WDS, delta_WDS); 
   fprintf(stderr, "(RA_Pcat=%.3f DEC_Pcat=%.3f Equinox=%d) Difference=%.3f H %.3f deg\n", 
         alpha_Pcat, delta_Pcat, (int)equinox_Pcat, D_alpha, D_delta); 
  } else {
  *is_OK = 1;
  }

return(0);
}
/************************************************************************
* Check the consistency between the WDS name and the object name 
*
* INPUT:
* ADS_name: ADS name in the PISCO catalog
* WDS_name: WDS name corresponding to the discoverer's name 
*
* OUTPUT:
* is_OK: 1 if consistent, 0 otherwise
************************************************************************/
int check_consistency_ADSname(char *ADS_name, char *WDS_name, 
                              char *discov_name, char *ADS_WDS_cross, 
                              int *is_OK)
{
char in_line[80], ADS_name0[20], WDS_name0[20], discov_name0[20];
int i, iline, found;
FILE *fp_ADS_WDS_cross;

/* If ADS is absent from object name (i.e. only discoverer's name)
* return from here:
*/
if(!*ADS_name) {
  *is_OK = 1;
  return(0);
  }

jlp_trim_string(ADS_name, 20);

/* Open input file containing the ADS_WDS cross-references */
if((fp_ADS_WDS_cross = fopen(ADS_WDS_cross, "r")) == NULL) {
   fprintf(stderr, "update_PISCO_catalog_main/Fatal error opening ADS/WDS cross-ref.: %s\n",
           ADS_WDS_cross);
   exit(-1);
  }

*is_OK = 0;

/* Scan the cross-reference file: */
/* Syntax:
BDS         ADS     Discovr Comp         WDS       Omit? 

13663     17158     A  1248           00000+7530
12704         1     STF3053  AB       00026+6606
13664     17180     BGH   1  AB-C     00024+1047

ADS numbers in fields 11 to 15
Discov names in fields 21 to 27
Components in fields 30 to 36
WDS numbers in fields 39 to 48
*/
iline = 0;
found = 0;
while(!feof(fp_ADS_WDS_cross)) {
  if(fgets(in_line, 120, fp_ADS_WDS_cross)) {
    iline++;
    if(in_line[0] != '%') {
/* Warning: C arrays start at 0, hence should remove one from field number: */
      for(i = 0; i < 6; i++) ADS_name0[i] = in_line[10+i];
      ADS_name0[6] = '\0';
      jlp_trim_string(ADS_name0, 20);
/* Warning: C arrays start at 0, hence should remove one from field number: */
      for(i = 0; i < 7; i++) discov_name0[i] = in_line[20+i];
      discov_name0[7] = '\0';
/* Warning: C arrays start at 0, hence should remove one from field number: */
      for(i = 0; i < 10; i++) WDS_name0[i] = in_line[38+i];
      WDS_name0[10] = '\0';
/* Compare ADS and WDS: */
      if(!strcmp(ADS_name, ADS_name0)) {
        if(!strcmp(WDS_name, WDS_name0)) {
           *is_OK = 1;
         } else {
           fprintf(stderr, "check_consistency_ADS/Error: ADS=%s corresponds to WDS=%s not to %s as indicated in your catalog!\n",
                 ADS_name, WDS_name0, WDS_name);
           fprintf(stderr, "(in your catalog: ADS%s = %s, in WDS cross-ref: ADS%s = %s)\n",
                 ADS_name, discov_name, ADS_name0, discov_name0);
           *is_OK = 0;
         }
        found = 1;
        break;
        }
      } /* EOF inline[0] == '%" */
    } /* EOF fgets */
} /* EOF while */

/* Close cross-reference file */
fclose(fp_ADS_WDS_cross);
return(0);
}
/*************************************************************************
* Search for the Hipparcos object located at the coordinates (alpha, delta)
* in the Hipparcos Input Catalog
* WARNING: this catalog contains parallaxes 
*          obtained BEFORE Hipparcos measurements!
*
   1-  6  I6     ---     HIC      [1/120313]+ Hipparcos Input Catalogue
                                    running number.
   8- 11  A4     ---     Comp     Component(s) considered in this entry
      13  A1     ---     Target  *[A-Hjg] Satellite target in case
                                    of joint entry
  15- 16  I2     h       RAh      Right ascension J2000 (hours), at Epoch
  18- 19  I2     min     RAm      Right ascension (minutes)
  21- 26  F6.3   s       RAs      Right ascension (seconds)
      28  A1     ---     DE-      Declination J2000 (sign)
  29- 30  I2     deg     DEd      Declination (degrees)
  32- 33  I2     arcmin  DEm      Declination (minutes)
  35- 39  F5.2   arcsec  DEs      Declination (seconds)
  41- 44  I4     a       Epoch   *Epoch for the position, generally 2000
 
* INPUT:
*  HIC_catalog: name of the Hipparcos Input catalog
*  alpha, delta, equinox: coordinates of the object to be searched for
*                         (alpha in hours and delta in degrees)
*
* OUTPUT:
*  HIP_name: Hipparcos name corresponding to (alpha, delta) 
*  CCDM_name: CCDM name corresponding to (alpha, delta) 
*  V_mag: V magnitude of the object
*  B_V_index: (B-V) index of the object
*  D_tolerance: tolerance of D_alpha/D_delta in degrees
*  found: 1 is object was found, 0 otherwise
*************************************************************************/
int search_object_in_HIC_catalog(char *HIC_catalog, double alpha, double delta, 
                                 double equinox, char *HIP_name, 
                                 char *CCDM_name, double *V_mag, 
                                 double *B_V_index, double D_tolerance, 
                                 int *found)
{
char cat_line0[512], c_sign[0], HIP_name0[20], buffer[512];
int ilen, iline, a1, a2, d1, d2, i_equinox0, nval;
double a3, alpha0, delta0, D_alpha, D_delta, fw;
float d3;
FILE *fp_HIC_cat;

/* Initialization: */
HIP_name[0] = '\0';
CCDM_name[0] = '\0';
*V_mag = 100.;
*B_V_index = 100.;

/* Open input file containing the HIC catalog */
if((fp_HIC_cat = fopen(HIC_catalog, "r")) == NULL) {
   fprintf(stderr, "search_discov_name_in_HIC_catalog/Fatal error opening HIC catalog: %s\n",
           HIC_catalog);
   return(-1);
  }

#ifdef DEBUG
printf("search_object_in_HIC_catalog: alpha=%f delta=%f equinox=%f D_tolerance=%f\n",
        alpha, delta, equinox, D_tolerance);
#endif

/* Look for the data concerning this object: */
*found = 0;
iline = 0;
while(!feof(fp_HIC_cat)) {
 if(fgets(cat_line0, 512, fp_HIC_cat)) {
   iline++;
   if(cat_line0[0] != '%') {
/*
   1-  6  I6     ---     HIC      [1/120313]+ Hipparcos Input Catalogue
                                    running number.
   8- 11  A4     ---     Comp     Component(s) considered in this entry
      13  A1     ---     Target  *[A-Hjg] Satellite target in case
                                    of joint entry
  15- 16  I2     h       RAh      Right ascension J2000 (hours), at Epoch
  18- 19  I2     min     RAm      Right ascension (minutes)
  21- 26  F6.3   s       RAs      Right ascension (seconds)
      28  A1     ---     DE-      Declination J2000 (sign)
  29- 30  I2     deg     DEd      Declination (degrees)
  32- 33  I2     arcmin  DEm      Declination (minutes)
  35- 39  F5.2   arcsec  DEs      Declination (seconds)
  41- 44  I4     a       Epoch   *Epoch for the position, generally 2000
*/
/* Warning: C arrays start at 0, hence should remove one from field number: */
      strncpy(HIP_name0, &cat_line0[0], 6);
      HIP_name0[6] = '\0';
/* Other values: */
      ilen = 44 - 15 + 1;
      strncpy(buffer, &cat_line0[14], ilen),
      buffer[ilen] = '\0';
      nval = sscanf(buffer, "%2d %2d %6lf %c %2d %2d %5f %4d",
                    &a1, &a2, &a3, c_sign, &d1, &d2, &d3, &i_equinox0);
      if(nval != 8) {
         fprintf(stderr, "Fatal error reading coordinates in Hipparcos input catalog/iline=%d nval=%d\n buffer=>%s<\n",
                 iline, nval, buffer);
         exit(-1);
         }
/* alpha in hours and delta in degrees */
      alpha0 = (double)a1 + ((double)a2)/60. + a3/3600.;
      delta0 = (double)d1 + ((double)d2)/60. + d3/3600.;
      if(c_sign[0] == '-') delta0 *= -1.;
/* DEBUG:
      if(iline < 10) {
      printf("%.80s \n", cat_line0);
      printf("%i %i %.3f >%c< %i %i %.3f %d (RA=%.3f DEC=%.3f)\n", 
             a1, a2, a3, c_sign[0], d1, d2, d3, i_equinox0, alpha0, delta0);
      }
EOF DEBUG*/
      D_alpha = ABS(alpha0 - alpha);
      D_delta = ABS(delta0 - delta);
/* D_tolerance in degrees (0.1 arcmin is a good value) */
      if(D_alpha < D_tolerance/15. && D_delta < D_tolerance) {
        strcpy(HIP_name, HIP_name0);
        jlp_trim_string(HIP_name, 7);
        *found = 1;
/* Read the other parameters: 
 191-196  F6.3   mag     Vmag     V magnitude
 198-201  F4.3   mag     e_Vmag   Error of V magnitude
 203-208  F6.3   mag     B-V      ?B-V colour index
 210-213  F4.3   mag     e_B-V    ?Error of B-V
     215  A1     ---     r_B-V   *Source of photometry (see Table B2)
 217-227  A11    ---     Sp      *Spectral type and luminosity class
     229  A1     ---     r_Sp    *Source of the spectral type (see Table B3)
 231-235  I5     mas     Plx      ?Parallax (milli-arcsec) *** JLP: before 1991!
 237-238  I2     mas     e_Plx    ?Probable error of parallax (milli-arcsec)
                                   *** JLP: before 1991 ! 
*/
       ilen = 6;
       strncpy(buffer, &cat_line0[190], ilen),
       buffer[ilen] = '\0';
       if((nval = sscanf(buffer, "%6lf", &fw)) == 1) *V_mag = fw; 
       ilen = 6;
       strncpy(buffer, &cat_line0[202], ilen),
       buffer[ilen] = '\0';
       if((nval = sscanf(buffer, "%6lf", &fw)) == 1) *B_V_index = fw; 
/*
 286-295  A10    ---     CCDM    *CCDM number (details in annex1)
 297-298  A2     ---     CCDMcomp Components considered
 300-302  A3     deg     PA      *[ 0-9.EFPNS] Position angle components
 304-309  F6.1   arcsec  Sep      ?Separation between the components considered.
                                  *** JLP: before 1991!
 311-314  F4.1   mag     Dmag     ?Magnitude difference between the components
                                  *** JLP: before 1991!
     316  A1     ---     n_CCDM  *[OA*] Information on orbital systems
                                  *** JLP: before 1991!
*/
       ilen = 13;
       strncpy(CCDM_name, &cat_line0[285], ilen),
       CCDM_name[ilen] = '\0';
#ifdef DEBUG
printf("search_object_in_HIC_catalog: OK found (HIP%s and CCDM%s) !\n", 
       HIP_name, CCDM_name);
#endif
// When found, break here:
       break;
       } // EOF found
   } /* EOF cat_line[0] != '%' */
 } /* EOF fgets... */
} /* EOF while */

/* Close Hipparcos input catalog: */
fclose(fp_HIC_cat);
return(0);
}
/*************************************************************************
* Read data from the Hipparcos/Tycho main Catalog
*
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
       1  A1    ---     Catalog   [H] Catalogue (H=Hipparcos)               (H0)
   9- 14  I6    ---     HIP       Identifier (HIP number)                   (H1)
      16  A1    ---     Proxy    *[HT] Proximity flag                       (H2)
  18- 28  A11   ---     RAhms     Right ascension in h m s, ICRS (J1991.25) (H3)
  30- 40  A11   ---     DEdms     Declination in deg ' ", ICRS (J1991.25)   (H4)
  42- 46  F5.2  mag     Vmag      ? Magnitude in Johnson V                  (H5)
      48  I1    ---     VarFlag  *[1,3]? Coarse variability flag            (H6)
      50  A1    ---   r_Vmag     *[GHT] Source of magnitude                 (H7)
  52- 63  F12.8 deg     RAdeg    *? alpha, degrees (ICRS, Epoch=J1991.25)   (H8)
  65- 76  F12.8 deg     DEdeg    *? delta, degrees (ICRS, Epoch=J1991.25)   (H9)
      78  A1    ---     AstroRef *[*+A-Z] Reference flag for astrometry    (H10)
  80- 86  F7.2  mas     Plx       ? Trigonometric parallax                 (H11)
  88- 95  F8.2 mas/yr   pmRA     *? Proper motion mu_alpha.cos(delta), ICRS(H12)
  97-104  F8.2 mas/yr   pmDE     *? Proper motion mu_delta, ICRS           (H13)
 106-111  F6.2  mas   e_RAdeg    *? Standard error in RA*cos(DEdeg)        (H14)
 113-118  F6.2  mas   e_DEdeg    *? Standard error in DE                   (H15)
 120-125  F6.2  mas   e_Plx       ? Standard error in Plx                  (H16)
  
* INPUT:
*  HIP_catalog: name of the Hipparcos/Tycho main catalog
*  HIP_name: Hipparcos name of the object 
*
* OUTPUT:
*  paral, err_paral: parallax and error on this parallax
*************************************************************************/
int read_data_in_HIP_catalog(char *HIP_catalog, char *HIP_name, 
                             double *paral, double *err_paral, int *found)
{
char cat_line0[512], HIP_name0[20], buffer[512];
int ilen, iline, nval;
double fw;
FILE *fp_HIP_cat;

/* Initialization: */
*paral = -1.;
*err_paral = -1.;

/* Open input file containing the HIP catalog */
if((fp_HIP_cat = fopen(HIP_catalog, "r")) == NULL) {
   fprintf(stderr, "read_data_in_HIP_catalog/Fatal error opening HIP catalog: %s\n",
           HIP_catalog);
   return(-1);
  }

/* Look for the data concerning this object: */
*found = 0;
iline = 0;
while(!feof(fp_HIP_cat)) {
 if(fgets(cat_line0, 512, fp_HIP_cat)) {
   iline++;
   if(cat_line0[0] != '%') {
/*
*/
/* Warning: C arrays start at 0, hence should remove one from field number: 
   9- 14  I6    ---     HIP       Identifier (HIP number)                   (H1)
*/
    strncpy(HIP_name0, &cat_line0[8], 6);
    HIP_name0[6] = '\0';
    jlp_compact_string(HIP_name0, 7);
    if(!strcmp(HIP_name0, HIP_name)) {
      *found = 1;
/* 80- 86  F7.2  mas    Plx      ? Trigonometric parallax                 (H11)
 120-125  F6.2  mas   e_Plx      ? Standard error in Plx                  (H16)
*/
      ilen = 7;
      strncpy(buffer, &cat_line0[79], ilen),
      buffer[ilen] = '\0';
      if((nval = sscanf(buffer, "%lf", &fw)) == 1) *paral = fw; 
      ilen = 6;
      strncpy(buffer, &cat_line0[119], ilen),
      buffer[ilen] = '\0';
      if((nval = sscanf(buffer, "%lf", &fw)) == 1) *err_paral = (double) fw; 
      break;
      }
   } /* EOF cat_line[0] != '%' */
 } /* EOF fgets... */
} /* EOF while */

/* Close Hipparcos main catalog: */
fclose(fp_HIP_cat);
return(0);
}
/************************************************************************
* Get the Hipparcos name from the WDS name 
* by reading the HIP/HDS/WDS cross-reference file from the WDS website
* (only for debugging, since there are very few objects in this list !)
*
* INPUT:
* WDS_name: WDS name 
* HIP_HDS_WDS_cross: name of the cross-catalog between HIP/HDS/WDS names
*
* OUTPUT:
* HIP_name: Hipparcos name corresponding to WDS_name
* founs: 1 if found, 0 otherwise
************************************************************************/
int HIP_name_from_HIP_HDS_WDS_cross(char *WDS_name, char *HIP_HDS_WDS_cross, 
                                    char *HIP_name, int *found)
{
int iline;
char line_buffer[80], WDS_name0[20];
FILE *fp_HIP_HDS_WDS_cross;

/* Initialization: */
HIP_name[0] = '\0';

/* Open HIP/WDS cross reference file: */
if((fp_HIP_HDS_WDS_cross = fopen(HIP_HDS_WDS_cross, "r")) == NULL) {
  fprintf(stderr, "HIP_name_from_HIP_HDS_WDS_cross/Fatal error opening %s\n",
          HIP_HDS_WDS_cross);
  exit(-1);
 }

/* Scan all the file looking for the object name */
*found = 0;
iline = 0;
while(!feof(fp_HIP_HDS_WDS_cross) && !(*found)){
  if(fgets(line_buffer, 80, fp_HIP_HDS_WDS_cross)) {
  iline++;
  if(line_buffer[0] != '%') {
    strncpy(WDS_name0, line_buffer, 10);
    WDS_name0[10] = '\0';
    if(!strncmp(WDS_name0, WDS_name,10)) {
      *found = 1;
/* Get Hipparcos number: */
      strncpy(HIP_name, &line_buffer[23], 6);
      HIP_name[6] = '\0';
      jlp_trim_string(HIP_name, 7);
      break;
      }
    } /* EOF inline[0] == '%" */
  } /* EOF fgets */
} /* EOF while */

/* Close cross-reference file */
fclose(fp_HIP_HDS_WDS_cross);
return(0);
}

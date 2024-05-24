/************************************************************************
* "WDS_catalog_utils.c"
* To retrieve data from the WDS catalog
*
* JLP 
* Version 20/11/2009
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>                /* exit() */
#include <string.h>
#include <ctype.h>                 /* isprint, isdigit ... */
#include <math.h>
#include <time.h>                  /* date */
#include "jlp_catalog_utils.h"     // Routines used to read catalogs:

#include "WDS_catalog_utils.h"     /* Prototypes of the routines defined here */ 
#include "HIP_catalog_utils.h"     
#include "jlp_string.h"

#define ABS(a) ((a) < 0.0  ? (-(a)) : (a))
#ifndef PI
#define PI 3.14159265
#endif
#define DEGTORAD   (PI/180.00)


/*
#define DEBUG_1
#define DEBUG
*/

/* Defined here:
int search_discov_name_in_WDS_catalog(char *WDS_catalog, char *discov_name,
                                      char *wds_name, int *found);
int get_data_from_WDS_catalog(char *WDS_catalog, char *discov_name,
                              char *comp_name,
                              char *wds_name, double *WdsLastYear, 
                              double *WdsLastRho, double *WdsLastTheta, 
                              double *WdsMagA, double *WdsMagB,
                              char *WdsSpectralType,
                              int *found);
int read_coordinates_from_WDS_catalog(char *wds_name, char *WDS_catalog, 
                                      char *str_alpha, char *str_delta,
                                      double *alpha, double *delta,
                                      double *equinox, int *found);
*/

/*************************************************************************
* Format of WDS catalog (version of 2009-2012):
*
  COLUMN     Format                     DATA
  --------   ------         ----------------------------
  1  -  10   A10             2000 Coordinates
  11 -  17   A7              Discoverer & Number
  18 -  22   A5              Components
  24 -  27   I4              Date (first)
  29 -  32   I4              Date (last)
  34 -  37   I4              Number of Observations (up to 9999)
  39 -  41   I3              Position Angle (first - XXX)
  43 -  45   I3              Position Angle (last  - XXX)
  47 -  51   F5.1            Separation (first)
  53 -  57   F5.1            Separation (last)
  59 -  63   F5.2            Magnitude of First Component
  65 -  69   F5.2            Magnitude of Second Component
  71 -  79   A9              Spectral Type (Primary/Secondary)
  81 -  84   I4              Primary Proper Motion (RA)
  85 -  88   I4              Primary Proper Motion (Dec)
  90 -  93   I4              Secondary Proper Motion (RA)
  94 -  97   I4              Secondary Proper Motion (Dec)
  99 - 106   A8              Durchmusterung Number
 108 - 111   A4              Notes
 113 - 130   A18             2000 arcsecond coordinates
*
* INPUT:
*  WDS_catalog: name of the WDS catalog
*  discov_name: discoverer's name of the object to be searched for
*
* OUTPUT:
*  wds_name: WDS name corresponding to discov_name
*  found: 1 is object was found, 0 otherwise
*************************************************************************/
int search_discov_name_in_WDS_catalog(char *WDS_catalog, char *discov_name,
                                      char *comp_name,
                                      char *wds_name, char *wds_discov_name,
                                      char *wds_comp_name, int *found)
{
FILE *fp_WDS_cat;
char cat_line0[256], discov_name0[20], comp_name0[20], wds_name0[20];
int iline, same_comp;

*found = 0;

/* Removes all the blanks since 7 characters for WDS, and 8 characters 
* for Marco's file */
jlp_compact_string(discov_name, 20);
jlp_compact_string(comp_name, 20);

/* Open input file containing the WDS catalog */
if((fp_WDS_cat = fopen(WDS_catalog, "r")) == NULL) {
   fprintf(stderr, "search_discov_name_in_WDS_catalog/Fatal error opening WDS catalog: %s\n",
           WDS_catalog);
   return(-1);
  }

/* Look for the data concerning this object: */
*found = 0;
iline = 0;
while(!feof(fp_WDS_cat) && (*found != 1)) {
 if(fgets(cat_line0, 256, fp_WDS_cat)) {
   iline++;
   if(cat_line0[0] != '%') {
/*   
  1  -  10   A10             2000 Coordinates
  11 -  17   A7              Discoverer & Number
  18 -  22   A5              Components
*/
   strncpy(wds_name0, &cat_line0[0], 10);
   wds_name0[10] = '\0';
/* 7 characters for WDS, 8 characters for Marco's file */
   strncpy(discov_name0, &cat_line0[10], 7);
   discov_name0[7] = '\0';
   strncpy(comp_name0, &cat_line0[17], 5);
   comp_name0[5] = '\0';

/* Removes all the blanks since 7 characters for WDS, and 8 characters 
* for Marco's file */
   jlp_compact_string(discov_name0, 20);
   jlp_compact_string(comp_name0, 20);
   same_comp = 0; 
   if(!strcmp(comp_name, comp_name0)) same_comp = 1; 
   if((comp_name[0] == '\0') && !strcmp(comp_name0, "AB")) same_comp = 1; 
   if((comp_name0[0] == '\0') && !strcmp(comp_name, "AB")) same_comp = 1; 
   if(!strcmp(discov_name, discov_name0)) {
     strcpy(wds_name, wds_name0);
     strcpy(wds_discov_name, discov_name0);
     strcpy(wds_comp_name, comp_name0);
     *found = 1;
   if(same_comp == 1) {
     *found = 2;
     break;
     }
#ifdef DEBUG_1
     printf("search_discov_name_in_WDS_catalog/Object found (discov_name=%s< comp_name=%s<) found=%d", 
             discov_name, comp_name,  *found);
     printf(" WDS=%s discov0=%s< comp0=%s<\n", wds_name0, discov_name0, comp_name0);
#endif
     }
   } /* EOF cat_line[0] != '%' */
  } /* EOF fgets... */
} /* EOF while */

fclose(fp_WDS_cat);
return(0);
}
/***********************************************************************
* Get miscellaneous data from the WDS catalog from discov_name and comp_name
*
* INPUT:
*  WDS_catalog: name of the WDS catalog
*  discov_name1: discoverer's name of the object to be searched for
*  comp_name1: companion name (if any) of thye object to be searched for
*
* OUTPUT:
*  wds_name: WDS name corresponding to discov_name
*  WdsLastYear, WdsLastRho, WdsLastTheta: year, rho, theta of the last 
*           observation reported in the WDS
*  WdsMagA, WdsMagB : magnitude of primary and secondary components
*  WdsSpectralType : spectral type of A and B (when known)
*  wds_meas_found: 1 is measures for this object has been found, 0 otherwise
***********************************************************************/
int get_data_from_WDS_catalog(char *WDS_catalog, char *discov_name1,
                              char *comp_name1,
                              char *wds_name, double *WdsLastYear, 
                              double *WdsLastRho, double *WdsLastTheta, 
                              double *WdsMagA, double *WdsMagB,
                              char *WdsSpectralType, int *wds_meas_found)
{
FILE *fp_WDS_cat;
char cat_line0[256], discov_name0[20], comp_name0[20], wds_name0[20];
char cvalue[64], notes0[64], orbit0, *pc, full_discov_name[64];
char full_discov_name0[64];
double dvalue;
int iline, ivalue;

// Copy input discov_name and comp name: 
// Handle case of AB companion 
jlp_compact_string(comp_name1, 20);
if(!strcmp(comp_name1, "AB")) {
  sprintf(full_discov_name, "%s", discov_name1);
  } else {
  sprintf(full_discov_name, "%s%s", discov_name1, comp_name1);
  }

/* DEBUG
printf("get_data_from_WDS_catalog: looking for discov=%s comp=%s\n", 
       discov_name1, comp_name1);
printf("get_data_from_WDS_catalog: looking for full discov=%s\n", full_discov_name);
*/

*WdsLastYear = 0.;
*WdsLastRho = 0.;
*WdsLastTheta = 0.;
*WdsMagA = 0.;
*WdsMagB = 0.;
WdsSpectralType[0] = '\0';
wds_name[0] = '\0';

/* Removes all the blanks since 7 characters for WDS, and 8 characters 
* for Marco's file */
jlp_compact_string(full_discov_name, 20);

/* Open input file containing the WDS catalog */
if((fp_WDS_cat = fopen(WDS_catalog, "r")) == NULL) {
   fprintf(stderr, "get_data_from_WDS_catalog/Fatal error opening WDS catalog: %s\n",
           WDS_catalog);
   return(-1);
  }

/* Look for the data concerning this object: */
*wds_meas_found = 0;
iline = 0;
while(!feof(fp_WDS_cat)) {
 if(fgets(cat_line0, 256, fp_WDS_cat)) {
   iline++;
   if(cat_line0[0] != '%') {
/*   
  1  -  10   A10             2000 Coordinates
  11 -  17   A7              Discoverer & Number
  18 -  22   A5              Components
  59 -  63   F5.2            Magnitude of First Component
  65 -  69   F5.2            Magnitude of Second Component
  71 -  79   A9              Spectral Type (Primary/Secondary)
  108 - 111  A4              Notes
*/
//** 1. WDS name: "1  -  10   A10"
   strncpy(wds_name0, &cat_line0[0], 10);
   wds_name0[10] = '\0';
//** 2. Discoverer & Number: "11 -  17   A7" 
// 7 characters for WDS, 8 characters for Marco's file
   strncpy(discov_name0, &cat_line0[10], 7);
   discov_name0[7] = '\0';
//** 3. Components: "18 -  22   A5" 
   strncpy(comp_name0, &cat_line0[17], 5);
   comp_name0[5] = '\0';
//** 4. Magnitude of First Component: "59 -  63   F5.2" 
   strncpy(cvalue, &cat_line0[58], 5);
   ivalue = sscanf(cvalue, "%lf", &dvalue);
   if(ivalue == 1) *WdsMagA = dvalue;
//** 5. Magnitude of Second Component: "65 -  69   F5.2" 
   strncpy(cvalue, &cat_line0[64], 5);
   ivalue = sscanf(cvalue, "%lf", &dvalue);
   if(ivalue == 1) *WdsMagB = dvalue;
//** 6. Spectral Type (Primary/Secondary): "71 -  79   A9" 
   strncpy(WdsSpectralType, &cat_line0[70], 9);
   WdsSpectralType[9] = '\0';
//** 7. Orbit
// O: Orbit. A published orbit exists
// Orbit, briefly described in WDSNOT MEMO and has entry in Orbit Catalog
   strncpy(notes0, &cat_line0[107], 4);
   notes0[4] = '\0';
   pc = strchr(notes0, 'O');
   if(pc != NULL) orbit0 = 1;
   else orbit0 = 0;

/* Removes all the blanks since 7 characters for WDS, and 8 characters 
* for Marco's file */
   jlp_compact_string(discov_name0, 20);
// Handle case of AB companion 
   jlp_compact_string(comp_name0, 20);
   if(!strcmp(comp_name0, "AB")) {
     sprintf(full_discov_name0, "%s", discov_name0);
   } else {
     sprintf(full_discov_name0, "%s%s", discov_name0, comp_name0);
   }
   jlp_compact_string(full_discov_name0, 64);
   if(!strcmp(full_discov_name, full_discov_name0)) {
// Copy to output value of wds_name:
     strcpy(wds_name, wds_name0);
#ifdef DEBUG_1
     printf("get_data_from_WDS_catalog/Full object name found in WDS catalog (discov_name =%s)\n", full_discov_name);
     printf("WDS=%s discov=%s comp=%s\n", wds_name0, discov_name0, comp_name0);
     printf("WDS_cat/notes0=%s pc_notes=%s orbit=%d\n", notes0, pc, orbit0);
#endif
/* Read WdsLastYear, WdsLastRho, WdsLastTheta: 
  24 -  27   I4              Date (first)
  29 -  32   I4              Date (last)
  34 -  37   I4              Number of Observations (up to 9999)
  39 -  41   I3              Position Angle (first - XXX)
  43 -  45   I3              Position Angle (last  - XXX)
  47 -  51   F5.1            Separation (first)
  53 -  57   F5.1            Separation (last)
*/
/* (last) year */
     strncpy(cvalue, &cat_line0[28], 4);
     cvalue[4] = '\0';
     if(sscanf(cvalue, "%d", &ivalue) == 1) *WdsLastYear = ivalue;

/* (last) theta */
     strncpy(cvalue, &cat_line0[42], 3);
     cvalue[3] = '\0';
     if(sscanf(cvalue, "%d", &ivalue) == 1) *WdsLastTheta = ivalue;

/* (last) rho */
     strncpy(cvalue, &cat_line0[52], 5);
     cvalue[5] = '\0';
     if(sscanf(cvalue, "%lf", &dvalue) == 1) *WdsLastRho = dvalue;

     *wds_meas_found = 1;
#ifdef DEBUG
     printf("last rho=%f theta=%f epoch=%f found=%d\n",
             *WdsLastRho, *WdsLastTheta, *WdsLastYear, *wds_meas_found);
#endif
     break;
     }
   } /* EOF cat_line[0] != '%' */
  } /* EOF fgets... */
} /* EOF while */

fclose(fp_WDS_cat);
return(0);
}
/***********************************************************************
* Look for accurate coordinates in WDS catalog
*
Columns   1- 10:    The  hours, minutes, and tenths of minutes of Right 
                    Ascension for 2000, followed by the degrees and minutes of
                    Declination for 2000, with + and - indicating north and
                    south declinations. The positions given represent our best
                    estimates of these values. Where possible, these are based
                    on the ACRS and PPM data, with proper motion incorporated.

Columns 113-130:    The hours, minutes, seconds and tenths of seconds (when
                    known) of Right Ascension for 2000, followed by the degrees,
                    minutes, and seconds of Declination for 2000, with + and - 
                    indicating north and south declinations. The positions given
                    represent our best estimates of these values. Where 
                    possible, these are based on the Hipparcos and Tycho data, 
                    with proper motion incorporated. While the arcminute 
                    coordinate (columns 1-10) refer to the primary of a multiple
                    system, the arcsecond coordinate (columns 113-130) refer to
                    the primary of the subsystem. For example, while the BC pair
                    of an A-BC multiple will have the same 10 digit WDS 
                    coordinate, the arcsecond coordinate of the BC pair will be
                    at the "B" position.
* Example:
06019+6052         in 1-10
060156.93+605244.8 in 113-130
*
* INPUT:
*  wds_name: WDS name of the object 
*  WDS_catalog: name of the WDS catalog
*
* OUTPUT:
* alpha, delta, equinox: coordinates of the object
*            (alpha in hours and delta in degrees)
************************************************************************/
int read_coordinates_from_WDS_catalog(char *wds_name, char *WDS_catalog, 
                                      char *str_alpha, char *str_delta,
                                      double *alpha, double *delta,
                                      double *equinox, int *found)
{
FILE *fp_WDS_cat;
char cat_line0[256], wds_name0[40], cvalue[64], sign[1];
int iline, hh, hm, hs, hss, dd, dm, ds, dss;

*alpha = 0.;
*delta = 0.;
*equinox = 2000.;

/* Removes all the blanks since 10 characters for WDS */ 
jlp_compact_string(wds_name, 40);

/* Open input file containing the WDS catalog */
if((fp_WDS_cat = fopen(WDS_catalog, "r")) == NULL) {
   fprintf(stderr, "read_coordinates_in_WDS_catalog/Fatal error opening WDS catalog: %s\n",
           WDS_catalog);
   return(-1);
  }

/* Look for the data concerning this object: */
*found = 0;
iline = 0;
while(!feof(fp_WDS_cat)) {
 if(fgets(cat_line0, 256, fp_WDS_cat)) {
   iline++;
   if(cat_line0[0] != '%') {
/*   
  1  -  10   A10             2000 Coordinates (= WDS name)
  11 -  17   A7              Discoverer & Number
  18 -  22   A5              Components
*/
   strncpy(wds_name0, &cat_line0[0], 10);
   wds_name0[10] = '\0';

   if(!strcmp(wds_name, wds_name0)) {
     *found = 1;
#ifdef DEBUG
     printf("read_coordinates_from_WDS_catalog/Object found in WDS catalog (wds_name =%s)\n", wds_name);
#endif

/* Read WdsLastYear, WdsLastRho, WdsLastTheta: 
  113 - 130   A18            2000 precise coordinates
Example:
060156.93+605244.8 in 113-130
*/
     strncpy(cvalue, &cat_line0[112], 18);
     strncpy(str_alpha, &cat_line0[112], 9);
     str_alpha[9] = '\0';
     strncpy(str_delta, &cat_line0[121], 9);
     str_delta[9] = '\0';
     cvalue[18] = '\0';
     if(sscanf(cvalue, "%02d%02d%02d.%02d%c%02d%02d%02d.%d", 
        &hh, &hm, &hs, &hss, sign, &dd, &dm, &ds, &dss) == 9) {
        *alpha = (double)hh + ((double)hm)/60. 
                + ((double)hs + (double)hss/10.)/3600.;
        *delta = (double)dd + ((double)dm)/60. 
                + ((double)ds + (double)dss/10.)/3600.;
        if(sign[0] == '-') *delta *= -1.;
        else if(sign[0] != '+') {
           fprintf(stderr,"read_coordinates/Fatal error: sign=%s\n", sign);
           exit(-1);
           }
#ifdef DEBUG_1
     printf("WDS%s : %02d%02d%02d.%02d%s%02d%02d%02d.%d\n",
            wds_name, hh, hm, hs, hss, sign, dd, dm, ds, dss); 
     printf("WDS%s : alpha=%f delta=%f \n", wds_name, *alpha, *delta);
#endif
     } else {
       fprintf(stderr,"read_coordinates_from_WDS_catalog/Error in catalog\n");
       fprintf(stderr,"Warning: error reading coordinates: >%s< of WDS%s in line%d\n",
               cvalue, wds_name, iline);
     } 

     *found = 1;
     break;
     }
   } /* EOF cat_line[0] != '%' */
  } /* EOF fgets... */
} /* EOF while */

fclose(fp_WDS_cat);
return(0);
}
/***********************************************************************
* Get miscellaneous data from the WDS catalog
*
* INPUT:
*  WDS_catalog: name of the WDS catalog
*  discov_name: discoverer's name of the object to be searched for
*
* OUTPUT:
*  wds_name: WDS name corresponding to discov_name
*  V_mag, B_V_index : magnitude and color index  from Hipparcos
*  paral, err_paral : parallax and error from Hipparcos
*  magV_A, magV_B : magnitude of primary and secondary components from WDS 
*  spectral_type : spectral types of primary and secondary components from WDS 
*  found_in_WDS: 1 is object was found in WDS, 0 otherwise
***********************************************************************/
int get_data_from_WDS_and_HIP_catalogs(char *WDS_catalog, 
                                       char *HIC_catalog, char *HIP_catalog, 
                                       char *discov_name, char *comp_name, 
                                       char *wds_name, 
                                       double *V_mag, double *B_V_index,
                                       double *paral, double *err_paral,
                                       double *magV_A, double *magV_B,
                                       char *spectral_type, int *found_in_WDS,
                                       int *found_in_Hip_cat)
{
char WDS_name[40], HIP_name[40], CCDM_name[40];
char *pc, buffer[128], str_alpha[64], str_delta[64];
double alpha_wds, delta_wds, equinox_wds, D_tolerance;
double WdsLastYear, WdsLastRho, WdsLastTheta;
int found, is_OK, status;

wds_name[0] = '\0';
spectral_type[0] = '\0';
*V_mag = 100.;
*B_V_index = 100.;
*paral = -1.;
*err_paral = 0.;
*magV_A = 100.;
*magV_B = 100.;
*found_in_WDS = 0;
*found_in_Hip_cat = 0;

/* Search for WDS number in WDS catalog using discov_name */
get_data_from_WDS_catalog(WDS_catalog, discov_name, comp_name,
                          wds_name, &WdsLastYear,
                          &WdsLastRho, &WdsLastTheta, magV_A, magV_B,
                          spectral_type, found_in_WDS);

if(*found_in_WDS != 0) {
 read_coordinates_from_WDS_catalog(wds_name, WDS_catalog, str_alpha,
                                   str_delta, &alpha_wds, 
                                   &delta_wds, &equinox_wds, &found);
 if(found != 1) {
   fprintf(stderr, "get_data_from_WDS_and_HIP_catalogs/Error: coorrds. of %s not found\n",
          wds_name);
   exit(-1);
   }

#ifdef DEBUG
printf("get_data_from_WDS_and_HIP_catalogs/discov_name=%s comp_name=%s wds_name=%s alpha_wds=%f delta_wds=%f equinox_wds=%f\n",
       discov_name, comp_name, wds_name, alpha_wds, delta_wds, equinox_wds);
#endif

/* Search for WDS number in HIPHIP//HDS/WDS cross-file
* (only for debugging, since there are very few objects in this list)
* NO longer used !
  HIP_name_from_HIP_HDS_WDS_cross(wds_name, HIP_HDS_WDS_cross, HIP_name_cross,
                                  &found_in_cross);
*/

/* Now only look for HIP object close to WDS coordinates
*/
/* Tolerance of the coordinates in degrees
* used for searching for Hipparcos names.
* 3 arcminutes is a good value with the coordinates 
* derived from the WDS names (hour-min+/-deg-min) */
// eg: 22070+3605 22 h 07 min, 36 deg 05 arcmin
     D_tolerance = 3./60.;

/* DEBUG: search in Hipparcos catalog */
    search_object_in_HIC_catalog(HIC_catalog, alpha_wds, delta_wds,
                                 equinox_wds, HIP_name, CCDM_name,
                                 V_mag, B_V_index, D_tolerance,
                                 found_in_Hip_cat);
#ifdef DEBUG
    printf("get_data_from_WDS_and_HIP_catalogs/DEBUG found_in_Hip_cat=%d\n",
            found_in_Hip_cat);
#endif

    if(*found_in_Hip_cat) {
#ifdef DEBUG
    printf("%s%s=%s was found in Hipparcos catalog (=HIP%s=CCDM%s)\n",
            discov_name, comp_name, wds_name, HIP_name, CCDM_name);
    printf("V=%.3f B-V=%.3f\n", *V_mag, *B_V_index);
#endif
    read_data_in_HIP_catalog(HIP_catalog, HIP_name, paral, err_paral,
                             &found);
    if(found) {
#ifdef DEBUG
       printf("Paral=%.2f+/-%.2f\n", *paral, *err_paral);
#endif
      } else {
       fprintf(stderr, "Fatal error: HIP=%s not in Hipparcos main catalog!\n", HIP_name);
       exit(-1);
      }
 // EOF if found_in_HIP_cat
#ifdef DEBUG
    } else {
      printf("discov_name=%s not found in HIC catalog!\n",
              discov_name);
#endif
    }
  }

return(status);
}
/*************************************************************************
*
* Example: 
01022+2705
19022$-$2705
**************************************************************************/
int decode_WDS_name(char *WDS_name, double *WDS_alpha, double *WDS_delta)
{
char *pc, wds_alpha[32], wds_delta[32];
int status, isign, ia1, ia2, id1, id2, nval;
*WDS_alpha = 0.;
*WDS_delta = 0.;

// First part:
strcpy(wds_alpha, WDS_name);
pc = wds_alpha;
while(*pc && isdigit(*pc) ) pc++;
strcpy(wds_delta, pc);
*pc = '\0';
nval = sscanf(wds_alpha, "%02d%03d", &ia1, &ia2);
 if(nval == 2) *WDS_alpha = (double)ia1 + (double)ia2 /1000.;

// Second part:
// + or $-$
pc = wds_delta;
if(*pc == '+') { 
    isign = +1;
    pc++;
  } else {
    if(*pc == '$') {
      pc++;
      if(*pc == '-') {
        pc++;
        isign = -1;
        }
        pc++;
     } else if(*pc == '-') {
        pc++;
        isign = -1;
     } else {
       printf("Decode_WDS_name/Fatal error decoding >%s< \n", WDS_name);
       exit(-1);
     }
  }

// Last part:
strcpy(wds_delta, pc);
  while(*pc) pc++;
  *pc = '\0';
nval = sscanf(wds_delta, "%02d%02d", &id1, &id2);
 if(nval == 2) *WDS_delta = isign *((double)id1 + (double)id2 /100.);

printf("wds_name: >%s< = >%s< + >%s< sign=%d\n", 
       WDS_name, wds_alpha, wds_delta, isign);

printf("WDS_alpha=%f WDS_delta=%f\n", *WDS_alpha, *WDS_delta); 

return(status);
}


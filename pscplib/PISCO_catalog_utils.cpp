/************************************************************************
* "PISCO_catalog_utils.c"
* Set of routines used to read PISCO catalog
*
* JLP 
* Version 05/06/2015
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>    // exit()
#include <string.h>    // strcpy()
#include "jlp_numeric.h" // JLP_QSORT, MINI, MAXI, PI, etc
#include "jlp_string.h"  // jlp_compact_string

/* The prototypes of routines included here
* are defined in "jlp_catalog_utils.h":
*/
#include "PISCO_catalog_utils.h"
#include "jlp_catalog_utils.h" // is_in_line()
#include "latex_utils.h" // latex_get_column_item() 

/*
#define DEBUG 
#define DEBUG_1 
*/
/************************************************************
* Search for the coordinates of an object 
* in the PISCO catalog containing the list of objects (used by TAV1.EXE)
*
* INPUT : 
* NameInPiscoCatalog: name of object
* PISCO_catalog_name: name of PISCO catalog ("zeiss_doppie.cat")
*
* OUTPUT : 
* alpha, delta: coordinates of object (in radians)
* epoch_o: epoch of observation
* coord_equinox: equinox corresponding to the coordinates 
*
*************************************************************/
int get_coordinates_from_PISCO_catalog(char *PISCO_catalog_name, 
                                       char *NameInPiscoCatalog,
                                       double *alpha, double *delta, 
                                       double *coord_equinox)
{
char in_line1[80], in_line2[80], *pc;
FILE *fp_cat;
int found, digit_found, object_len, status;
char compacted_object_name[40], object_in_catalog[40];

strcpy(compacted_object_name, NameInPiscoCatalog);
jlp_compact_string(compacted_object_name, 40);

if((fp_cat = fopen(PISCO_catalog_name, "r")) == NULL) {
  fprintf(stderr, "get_coordinates_from_PISCO_catalog/Error opening %s\n", 
          PISCO_catalog_name);
  return(-1);
 }

/* Compute the useful length of the compacted_object_name: 
* ADS123AB should be 6 only 
* COU12AB should be 5 only 
* ADS123AC should also be 8, but the discovery component should be AC ! 
* (this is not useful for coordinates...)
*/
pc = compacted_object_name;
object_len = 0;
digit_found = 0;
while(*pc) {
  if(isdigit(*pc)) digit_found = 1;
  if(!isdigit(*pc) && digit_found) break; 
  pc++; 
  object_len++;
  }

/* Look for the data concerning this object: */
found = 0;
while(!feof(fp_cat)) {
  if(fgets(in_line1,80,fp_cat)) {
  if(in_line1[0] != '%') {
    strncpy(object_in_catalog, in_line1, 20); 
    object_in_catalog[19] = '\0';
    pc = object_in_catalog;

digit_found = 0;
while(*pc) {
  if(isdigit(*pc)) digit_found = 1;
  if(!isdigit(*pc) && digit_found) break;
  pc++; 
  }
*pc = '\0';
jlp_compact_string(object_in_catalog, 40);

/* ADS 213 is sometimes written as: ADS 213AB
* but not as ADS 2136
*/
    if(!strncmp(object_in_catalog, compacted_object_name, object_len)){
       if(fgets(in_line2,80,fp_cat)) {
         found = 1;
         break;
         } else {
         fprintf(stderr, "Error reading second line of object data\n");
         break;
         }
     } 
   } /* EOF in_line1 != % */ 
 } /* EOF fgets */
} /* EOF while */

if(found) {
  status = read_coordinates_from_PISCO_catalog(in_line1, alpha, delta, 
                                               coord_equinox);
  if(status) {
    fprintf(stderr, "Error reading PISCO_catalog at line: >%s< (status=%d)\n", 
            in_line1, status);
    } else {
#ifdef DEBUG
    printf("object=%s object_len=%d alpha=%f delta=%f \n %s\n %s\n", 
           NameInPiscoCatalog, object_len, *alpha, *delta, 
           in_line1, in_line2);
#endif
    status = 0;
    }
  } else {
  status = -1;
  }

/* Close PISCO_catalog */
fclose(fp_cat);
return(status);
}
/************************************************************
* Search for photometric data of an object 
* in the PISCO catalog containing the list of objects (used by TAV1.EXE)
*
* INPUT : 
* NameInPiscoCatalog: name of object (ADS or discoverer's name)
* PISCO_catalog_name: name of PISCO catalog ("zeiss_doppie.cat")
*
* OUTPUT : 
*
*************************************************************/
int get_data_from_PISCO_catalog(char *PISCO_catalog_name, 
                                char *NameInPiscoCatalog, 
                                double *magV, double *B_V, 
                                double *paral, double *err_paral, 
                                double *magV_A, double *magV_B, 
                                char *spectral_type, char *discov_name, 
                                char *comp_name,
                                char *ads_name, char *WDS_name) 
{
char in_line1[128], in_line2[128], compacted_line[128], *pc;
char compacted_object[64], compacted_discov[64], new_ads_name[80];
FILE *fp_cat;
int iline, found; 
int object_len, status, verbose_if_error;

*magV_A = 100.;
*magV_B = 100.;
spectral_type[0] = '\0';
discov_name[0] = '\0';
comp_name[0] = '\0';
*magV = 100.;
*paral = 0.;
*err_paral = 0.;

if((fp_cat = fopen(PISCO_catalog_name, "r")) == NULL) {
  fprintf(stderr, "get_data_from_PISCO_catalog/Error opening %s\n", 
          PISCO_catalog_name);
  return(-1);
 }

/* Compute the useful length of the object_name: 
* ADS123AB should be 6 only 
* COU12Ac should be 5 only 
* But ADS123AC should also be 8, but the discovery component should be AC !
*/
strcpy(compacted_object, NameInPiscoCatalog);
jlp_compact_string(compacted_object, 64);
pc = compacted_object;
object_len = 0;
while(*pc) {
/* To upper case: not used since it modifies Ab-B components for example 
*  if(islower(*pc)) *pc = toupper(*pc);
*/
  pc++; 
  object_len++;
  }

/* Look for the data concerning this object: 
ADS 16198 (2.50, 58) & 22 42 29.9 & 39 18 39 & 2000.0
& A: 9.5, B: 9.9      & BU   176,       & WDS22425+3917 \cr
*/
found = 0;
iline = 0;

while(!feof(fp_cat)) {
  if(fgets(in_line1,128,fp_cat)) {
  iline++;
  if(in_line1[0] != '%' && in_line1[0] != '&') {
/* Detect ADS name and store it,
 (needed in case filename corresponds to discoverer's name): */
    new_ads_name[0] = '\0';
    verbose_if_error = 0;
    latex_get_column_item(in_line1, new_ads_name, 1, verbose_if_error);
    pc = new_ads_name;
/* Remove the heading blanks */
    while(*pc && *pc != ' ') pc++;
    pc++;
/* Stop at the blanks encontered in the middle */
    while(*pc && *pc != ' ') pc++;
    *pc = '\0';
/* ADS 213 is sometimes written as: ADS 213AB
* but not as ADS 2136
*/
    strcpy(compacted_line, in_line1);
    jlp_compact_string(compacted_line,128);
    if(!strncmp(compacted_line, compacted_object, object_len) 
       && !isdigit(compacted_line[object_len]) ){
         strcpy(ads_name, new_ads_name);
         found = 1;
/* Read 2nd line if object found in first line: */
         if(!fgets(in_line2,128,fp_cat)) {
           fprintf(stderr, "Error reading second line of object data\n");
/* Close PISCO_catalog */
           fclose(fp_cat);
           return(-1);
           }
/* Look for the discoverer's name the 2nd line */
         PISCO_catalog_read_line2_discov(in_line2, discov_name, comp_name);
         strcpy(compacted_discov, discov_name);
         jlp_compact_string(compacted_discov, 64);
         break;
     } /* EOF !strncmp(compacted_line, compacted_object, object_len) */
/* Case of second line: */
   } else if(in_line1[0] == '&') {
    strcpy(in_line2, in_line1);
/* Look for the discoverer's name the 2nd line */
    PISCO_catalog_read_line2_discov(in_line2, discov_name, comp_name);
    strcpy(compacted_discov, discov_name);
    jlp_compact_string(compacted_discov, 64);
    if(!strncmp(compacted_object, compacted_discov, object_len)
        && !isdigit(compacted_discov[object_len])) {
      fprintf(stderr, "get_data_from_PISCO_catalog/Object %s found in 2nd line!\n",
              NameInPiscoCatalog);
      if(new_ads_name[0] != '\0') {
              fprintf(stderr, "Hence the object name %s is corrected to: %s\n",
              NameInPiscoCatalog, new_ads_name);
         strcpy(NameInPiscoCatalog, new_ads_name);
         }
      strcpy(ads_name, new_ads_name);
      found = 1;
      break;
    }
  } /* EOF in_line1[0] == '&' */
 } /* EOF fgets */
} /* EOF while */


if(found) {
  status = PISCO_catalog_read_line2(in_line2, magV_A, magV_B, spectral_type); 
  if (status) {
      fprintf(stderr, "From read_line2/Fatal error/Bad syntax in line #%d\n", iline);
      exit(-1);
      }

  status = PISCO_catalog_read_line2_WDS(in_line2, WDS_name);
  if (status) {
      fprintf(stderr, "From read_line2_WDS/Fatal error/Object=%s Bad syntax in line #%d\n", NameInPiscoCatalog, iline);
      exit(-1);
      }

  status = PISCO_catalog_read_line2_Hip(in_line2, magV, B_V, paral, err_paral);
/*
  if (status) {
      fprintf(stderr, "From read_line2_Hip/Fatal error/Bad syntax in line #%d\n", iline);
      exit(-1);
      }
*/

  status = 0;
 
/*
  status = read_coordinates_from_PISCO_catalog(in_line1, alpha, delta, 
                                               coord_equinox);
  if(status) {
    fprintf(stderr, "Error reading PISCO_catalog at line: >%s< (status=%d)\n", 
            in_line1, status);
    } 
*/
/* EOF case: object was found */
} else  {
/* Case: object was not found */
 status = -1;
}

/* Close PISCO_catalog */
fclose(fp_cat);
return(status);
}
/*********************************************************************
* Read coordinates in a line extracted from 
* the PISCO catalog containing the list of objects (used by TAV1.EXE)
*
* Syntax:
* name & 23 12 03.5 & -23 34 21 & 2000 
*
* INPUT:
* in_line: line containing the data
*
* OUTPUT:
* alpha, delta: coordinates of object
* coord_equinox: corresponding equinox
*********************************************************************/
int read_coordinates_from_PISCO_catalog(char *in_line, double *alpha, 
                                        double *delta, double *coord_equinox)
{
int a1, a2, d1, d2, nval, isign;
double a3, d3;
char *pc;

*alpha = 0.;
*delta = 0.;
*coord_equinox = 0.;

/* Skip the first column: */
pc = in_line;
while(*pc != '&' && *pc) pc++;
if(*pc == '\0') return(-1);

/* Read alpha in the 2nd column: */
pc++;
nval = sscanf(pc, "%d %d %lf", &a1, &a2, &a3);
if(nval != 3) return(-2);

/* Conversion from hours to radians: */
*alpha = (a1 + (double)a2 / 60. + a3 / 3600.) * PI / 12.;

/* Moves to the 3rd column: */
while(*pc != '&' && *pc) pc++;
if(*pc == '\0') return(-3);

/* Read delta in the 3rd column: */

/* First detect the minus sign if present: */
pc++;
while(*pc == ' ') pc++;
isign = +1;
if(*pc == '-') {
  isign = -1;
  pc++;
} else if (*pc == '+') pc++;

nval = sscanf(pc, "%d %d %lf", &d1, &d2, &d3);
if(nval != 3) return(-4);

/* Conversion from degrees to radians: */
*delta = isign * (d1 + (double)d2 / 60. + d3 / 3600.) * DEGTORAD;

/* Moves to the 4th column: */
while(*pc != '&' && *pc) pc++;
if(*pc == '\0') return(-5);

/* Read equinox in the 4th column: */
pc++;
nval = sscanf(pc, "%lf", coord_equinox);
if(nval != 1) return(-6);

#ifdef DEBUG_1
printf("OK: alpha = %d %d %.2f (=%f radians)\n", a1, a2, a3, *alpha);
printf("delta = %d %d %.2f sign=%d (=%f radians)\n", d1, d2, d3, isign, *delta);
printf("equinox= %.2f\n", *coord_equinox);
#endif

return(0);
}
/***************************************************************************
* Read the first line of PISCO catalog (used by TAV1.EXE)
*
* Example:
* ADS 326AB (3.40, 340) &  0 23 59.8 &-03 28 31 & 2000.0
* or:
* COU 653  (0.42, 257) &  0 26 58.0 & 30 57 47 & 2000.0
*
* INPUT:
* in_line: character string to be read
*
* OUTPUT:
* discov_name: name of this object as given by its discoverer 
*              using the WDS syntax (with 8 characters)
* object_name: object name (ADS 326AB, COU 653 for instance) 
* comp_name: name of the components (e.g., Aa, AB, AC, etc)
* alpha_Pcat, delta_Pcat, equinox_Pcat: coordinates of the object
*
* RETURN: -1 if error
*         10 if reference star
***************************************************************************/
int PISCO_catalog_read_line1(char *in_line, char *object_name, char *comp_name,
                             char *ADS_name, double *alpha_Pcat, 
                             double *delta_Pcat, double *equinox_Pcat)
{
char buffer[120], *pc;
int status, bad_syntax, nval, i1, i2, i3, verbose_if_error = 1;
double d3, sign;

object_name[0] = '\0';
comp_name[0] = '\0';

/* Extract the object name  (e.g., ADS 1234AB, COU 128, ADS 234, or HIP 12322) */
strcpy(buffer, in_line);
pc = buffer;
bad_syntax = 1;
while(*pc && *pc != ' ') pc++;
/* First blank: */
if(*pc == ' ') {
  pc++;
/* Second blank: */
  while(*pc && *pc != ' ' && (isalnum(*pc) || *pc == '-')) pc++;
  *pc = '\0';
  bad_syntax = 0;
  }
if(bad_syntax) {
  fprintf(stderr, "PISCO_catalog_read_line1/Bad syntax >%s<\n", in_line);
  return(-1);
  } 

/* Load object_name: */
strcpy(object_name, buffer);

/* If HIP xxxx, assumes it is a standard star used for focusing*/
if(!strncmp(object_name, "HIP", 3)) {
  return(10);
  }

/* Retrieve the ADS name: */ 
ADS_name_from_object_name(object_name, ADS_name);

/* Extract the companion name  Aa-C, AB, etc*/
strcpy(buffer, in_line);
pc = buffer;
/* Skip the first blank: */
while(*pc && *pc != ' ') pc++;
if(*pc == ' ') {
 pc++;
 while(*pc && isdigit(*pc)) pc++;
 if(*pc != ' ') {
  strcpy(buffer, pc);
  pc = buffer;
  while(*pc && *pc != ' ' && (*pc == '-' || isalpha(*pc))) pc++;
  *pc = '\0';
  strcpy(comp_name, buffer);
  }
}

/* Read coordinates: 
* Example: 
ADS 5042 (0.56, 259) &  6 25 34.2 & 22 27 28 & 2000.0
*/
*alpha_Pcat = 0.;
*delta_Pcat = 0.;
*equinox_Pcat = 0.;

/* Right ascension in 2nd column: */
status = latex_get_column_item(in_line, buffer, 2, verbose_if_error);
if(status) {
  fprintf(stderr, "PISCO_catalog_read_line1/Error: missing RA column\n");
  return(-1);
  }
nval = sscanf(buffer,"%d %d %lf", &i1, &i2, &d3);
if(nval != 3) {
  fprintf(stderr, "PISCO_catalog_read_line1/Error: bad right ascension!\n");
  return(-1);
  }
*alpha_Pcat = (double)i1 + ((double)i2)/60. + d3/3600.;

/* Declination in 3rd column: */
status = latex_get_column_item(in_line, buffer, 3, verbose_if_error);
if(status) {
  fprintf(stderr, "PISCO_catalog_read_line1/Error: missing DEC column\n");
  return(-1);
  }
/* Detect the sign first 
Example:
 -00 45 32
*/
pc = buffer;
while(*pc && *pc != '-') pc++;
if(*pc == '-') sign = -1.;
else sign = +1.;

nval = sscanf(buffer,"%d %d %d", &i1, &i2, &i3);
if(nval != 3) {
  fprintf(stderr, "PISCO_catalog_read_line1/Error: bad declination!\n");
  return(-1);
  }
/* Use only absolute values, 
* since the sign is taken into account with the value of "sign"*/
if(i1 < 0) i1 *= -1;
*delta_Pcat = sign * ((double)i1 + ((double)i2)/60. + ((double)i3)/3600.);

/* Equinox in 4th column: */
status = latex_get_column_item(in_line, buffer, 4, verbose_if_error);
if(status) {
  fprintf(stderr, "PISCO_catalog_read_line1/Error: missing Equinox column\n");
  return(-1);
  }
nval = sscanf(buffer,"%lf", &d3);
if(nval != 1) {
  fprintf(stderr, "PISCO_catalog_read_line1/Error: bad equinox!\n");
  return(-1);
  }
*equinox_Pcat = d3;

return(0);
}
/***************************************************************************
* PISCO catalog (list of targets used by TAV1.EXE)
* Read the first parameters contained in the second line 
* (not the Hipparcos parameters)
*
* Example:
* & A: 5.3, B: 6.6 F2   & COU   14, orbit & \cr
* or:
* & A: 6.0, B: 8.2 B2   & STT  457,       & \cr
* or:
* & A: 9.1, B: 9.8 K0   & HLD   60, orbit & WDS00014+3937 HIP110 CCDM00014+3937AB V=8.61 B-V=0.79 PI=20.42(1.91) \cr
*
* INPUT:
* in_line: character string to be read
*
* OUTPUT:
* magV_A, magV_B: V magnitude of A and B components
* spectral_type: (global) spectral type of the binary
***************************************************************************/
int PISCO_catalog_read_line2(char *in_line, double *magV_A, double *magV_B, 
                             char *spectral_type)
{
char buffer[120], buffer0[120], *pc, *pc0;
double val0;
int status, nval, verbose_if_error = 1;

/* Initial values: */
spectral_type[0] = '\0';
*magV_A = 100.;
*magV_B = 100.;

/* Read magV_A, magV_B, spectral_type 
* Example: & A: 9.5, B: 9.5 F0   & STF  261,    & \cr
* or: & A: 9.5, B: 9.5 F0   & STF  261, orbit   & \cr
*/
status = latex_get_column_item(in_line, buffer, 2, verbose_if_error);
if(status) {
  fprintf(stderr, "PISCO_catalog_read_line2/Error: missing column\n");
  return(-1);
  }
pc = buffer;
while(strncmp(pc, "A:", 2) && *pc) pc++;
if(!strncmp(pc, "A:", 2)) {
 pc += 2;
/* Save following characters to buffer0: */
 strcpy(buffer0, pc);
 pc0 = buffer0;
/* Cut the string when the first "," or alpha character is found: */
 while(*pc0 !=',' && !isalpha(*pc0) && *pc0) pc0++;
 *pc0 = '\0';
 if((nval = sscanf(buffer0, "%lf", &val0)) == 1) *magV_A = val0;

/* Decode magV_B now */
 while(strncmp(pc, "B:", 2) && *pc) pc++;
 if(!strncmp(pc, "B:", 2)) {
   pc += 2;
/* Remove trailing blanks: */
   while(*pc == ' ' && *pc) pc++;
/* Save following characters to buffer0: */
   strcpy(buffer0, pc);
   pc0 = buffer0;
   while(*pc0 != ' ' && *pc0) pc0++;
   *pc0 = '\0';
   if((nval = sscanf(buffer0, "%lf", &val0)) == 1) *magV_B = val0;
/* Extract spectral type: */
   while(!isalpha(*pc) && *pc) pc++;
   strcpy(spectral_type, pc);
   } /* EOF !strncmp(pc, "B:", 2) */
 } /* EOF !strncmp(pc, "A:", 2) */

return(0);
}
/***************************************************************************
* PISCO catalog (list of targets used by TAV1.EXE)
* Read the discoverer's name contained in the second line 
*
* Example:
* & A: 5.3, B: 6.6 F2   & COU   14, orbit & \cr
* or:
* & A: 6.0, B: 8.2 B2   & STT  457,       & \cr
* or:
* & A: 9.1, B: 9.8 K0   & HLD   60, orbit & WDS00014+3937 HIP110 CCDM00014+3937AB V=8.61 B-V=0.79 PI=20.42(1.91) \cr
*
* INPUT:
* in_line: character string to be read
*
* OUTPUT:
* discov_name: name of this object as given by its discoverer 
*              using the WDS syntax (with 8 characters)
***************************************************************************/
int PISCO_catalog_read_line2_discov(char *in_line, char *discov_name,
                                    char *comp_name)
{
char buffer[120], *pc;
int i, status, verbose_if_error = 0;

/* Initial values: */
discov_name[0] = '\0';
comp_name[0] = '\0';

/* Read discoverer's name: */
status = latex_get_column_item(in_line, buffer, 3, verbose_if_error);
if(status) {
/*
  fprintf(stderr, "PISCO_catalog_read_line2_discov/Error: missing column\n");
*/
  return(-1);
  }

/* Extract the discoverer's name (with a comma at the end) 
* Example: & A: 9.5, B: 9.5 F0   & STF  261,    & \cr
* or: & A: 9.5, B: 9.5 F0   & STF  261, orbit   & \cr
*/
pc = buffer;
/* Removes the leading blanks: */
while(*pc && *pc == ' ') pc++;

/* Stops when ',' or '&' is encountered, with a maximum of 10 characters: */
i = 0;
while(i < 10 && *pc && *pc != ',') discov_name[i++] = *(pc++); 
discov_name[i] = '\0';

return(0);
}
/***************************************************************************
* PISCO catalog (list of targets used by TAV1.EXE)
* Read the WDS name contained in the second line 
*
* Example:
* & A: 9.1, B: 9.8 K0   & HLD   60, orbit & WDS00014+3937 HIP110 CCDM00014+3937AB V=8.61 B-V=0.79 PI=20.42(1.91) \cr
*
* INPUT:
* in_line: character string to be read
*
* OUTPUT:
* WDS_name: name of this object in Washington Double Star data base 
*            (coordinates at equinox 2000)
***************************************************************************/
int PISCO_catalog_read_line2_WDS(char *in_line, char *WDS_name)
{
char buffer[120], *pc;
int status, verbose_if_error = 0;

/* Initial values: */
WDS_name[0] = '\0';

/* Extract the WDS name from the 4th column:  
& A: 5.4, B: 8.4 K1 & KUI 66,  &  WDS14148+1006 HIP69612 CCDM14148+1006AB V=5.29 B-V=1.01
*/
status = latex_get_column_item(in_line, buffer, 4, verbose_if_error);
if(status) {
/*
  fprintf(stderr, "PISCO_catalog_read_line2_WDS/Error: missing column\n");
*/
  return(-1);
  }
pc = buffer;
/* Removes the leading blanks: */
while(*pc && *pc == ' ') pc++;

/* Look for 'WDS': */
while(*pc && strncmp(pc, "WDS", 3)) pc++; 
if(!strncmp(pc, "WDS", 3)) {
  pc += 3;
  strncpy(WDS_name, pc, 10);
  WDS_name[10] = '\0';
  }

return(0);
}
/***************************************************************************
* PISCO catalog (list of targets used by TAV1.EXE)
* Read the Hipparcos parameters contained in the second line 
*
* Example:
* & A: 5.8, B: 8.9 G7   & BU   733, orbit & WDS00022+2705 HIP171 CCDM00021+2706AB V=5.80 B-V=0.69 PI=80.63(3.03) \cr
*
* INPUT:
* in_line: character string to be read
*
* OUTPUT:
* magV, B_V, paral, err_paral
***************************************************************************/
int PISCO_catalog_read_line2_Hip(char *in_line, double *magV, double *B_V, 
                                 double *paral, double *err_paral)
{
char buffer[120], buffer0[120], buffer1[120], *pc, *pc0, *pc1;
double val0;
int status, nval, verbose_if_error = 0;

/* Initial values: */
*magV = 100.;
*B_V = 100.;
*paral = 0.;
*err_paral = 0.;

status = latex_get_column_item(in_line, buffer, 4, verbose_if_error);
if(status) {
/*
  fprintf(stderr, "PISCO_catalog_read_line2_Hip/Error: missing column\n");
*/
  return(-1);
  }

/* Decode V= */
pc = buffer;
while(strncmp(pc, " V=", 3) && *pc) pc++;
if(!strncmp(pc, " V=", 3)) {
 pc += 3;
/* Save following characters to buffer0: */
 strcpy(buffer0, pc);
 pc0 = buffer0;
/* Cut the string when the first " " or alpha character is found: */
 while(*pc0 !=' ' && !isalpha(*pc0) && *pc0) pc0++;
 *pc0 = '\0';
 if((nval = sscanf(buffer0, "%lf", &val0)) == 1) *magV = val0;
}

/* Decode B-V= */
pc = buffer;
while(strncmp(pc, "B-V=", 4) && *pc) pc++;
if(!strncmp(pc, "B-V=", 4)) {
 pc += 4;
/* Save following characters to buffer0: */
 strcpy(buffer0, pc);
 pc0 = buffer0;
/* Cut the string when the first " " or alpha character is found: */
 while(*pc0 !=' ' && !isalpha(*pc0) && *pc0) pc0++;
 *pc0 = '\0';
 if((nval = sscanf(buffer0, "%lf", &val0)) == 1) *B_V = val0;
}

/* Decode PI= */
pc = buffer;
while(strncmp(pc, "PI=", 3) && *pc) pc++;
if(!strncmp(pc, "PI=", 3)) {
 pc += 3;
/* Save following characters to buffer0: */
 strcpy(buffer0, pc);
 pc0 = buffer0;
/* Cut the string when the first " " or alpha character is found: */
 while(*pc0 !=' ' && *pc0 !='(' && !isalpha(*pc0) && *pc0) pc0++;
/* Save the error to buffer1 if '('present: */
 if(*pc0 == '(') { 
   strcpy(buffer1, pc0+1);
   pc1 = buffer1;
   while(*pc1 !=' ' && *pc0 !=')' && !isalpha(*pc1) && *pc1) pc1++;
   *pc1 = '\0';
   if((nval = sscanf(buffer1, "%lf", &val0)) == 1) *err_paral = val0;
   }
/* Decode parallax: */
 *pc0 = '\0';
 if((nval = sscanf(buffer0, "%lf", &val0)) == 1) *paral = val0;
}

/* DEBUG: 
printf("%s\n magV=%.2f B_V=%.2f paral=%.4f +/- %.4f\n",
        in_line, *magV, *B_V, *paral, *err_paral);
*/

return(0);
}

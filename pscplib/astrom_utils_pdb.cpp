/*************************************************************************
* astrom_utils_pdb.cpp
* JLP
* Version 20/09/2020
*************************************************************************/
#include "ctype.h"            /* isdigit() */
#include "string.h"            /* strcpy() */
#include "astrom_utils1.h"    // astrom_read_object_data_for_name_only()
#include "astrom_utils2.h"
#include "jlp_catalog_utils.h" /* get_data_from_PISCO_catalog */
#include "WDS_catalog_utils.h" /* get_data_from_WDS_catalog */
#include "jlp_fitsio.h"  // get_epoch_from_fits_file()
#include "jlp_numeric.h" // JLP_QSORT, MINI, MAXI, PI, etc
#include "jlp_string.h"
#include "latex_utils.h"  // latex_read_svalue...

#include "astrom_utils_pdb.h"
/*
int astrom_add_WDS_from_discov(FILE *fp_in, FILE *fp_out, char *WDS_catalog);
int astrom_add_discov_and_WDS(FILE *fp_in, FILE *fp_out, char *WDS_catalog, 
                              char *ADS_WDS_cross);
int astrom_decode_new_object_name(char *b_in, char *ads_name, 
                                  char *discov_name, char *comp_name, 
                                  char *wds_name, double *year);
int get_discov_from_ads_wds_crossref(char *ADS_WDS_cross, char *ads_name1,
                                     char *discov_name1, char *comp_name1,
                                     char *wds_name1);
int get_ads_from_ads_wds_crossref(char *ADS_WDS_cross, char *discov_name1,
                                  char *ads_name1);
*/

#define DEBUG
/*
*/

/*************************************************************************
* Scan the astrom file and add various parameters 
* WDS number,
* WY, WT, WR (year, theta and rho of the last observation) 
* are read from WDS catalog
*
* INPUT:
* fp_in: input old astrom file pointer
* fp_out: output new astrom file pointer
* WDS_cat: wdsweb_summ.txt (used to obtain the last observations)
*
**************************************************************************/
int astrom_add_WDS_from_discov(FILE *fp_in, FILE *fp_out, char *WDS_catalog) 
{
double magV, B_V, paral, err_paral, magV_A, magV_B; 
double year, WdsLastYear, WdsLastTheta, WdsLastRho, mag_A, mag_B;
int iline, status, gili_format;
int found, found1;
char b_in[256], discov_name[64], *pc;
char wds_name[64], wds_discov_name[64], wds_comp_name[64];
char comp_name[64], spectral_type[64]; 

year = 0.;
discov_name[0] = '\0';
wds_name[0] = '\0';

iline = 0;
while(!feof(fp_in))
{
/* Maximum length for a line will be 170 characters: */
  if(fgets(b_in,170,fp_in))
  {
  iline++;
  b_in[255] = '\0';
/* Remove ^M (Carriage Return) if present: */
  pc = b_in;
  while(*pc) {
  if(*pc == '\r') *pc = ' ';
  pc++;
  }
 
  if((b_in [0] != '\\') && (b_in[0] != '%')) {
/* Check if input LateX line contains the name of the object:
*/
  gili_format = 0;
  status = astrom_read_object_data_for_name_only(b_in, discov_name, 
                                                 comp_name, &year,
                                                 gili_format);
#ifdef DEBUG
 printf("from astrom_read_object_data_for_name_only with %s status=%d\n", 
          b_in, status);
#endif
/* Remove the extra-blanks: */
  jlp_trim_string(discov_name, 64);
  jlp_trim_string(comp_name, 64);

  if(status == 0) {

#ifdef DEBUG
    printf("DEBUG: %s From decode: discov=%s< comp_name=%s< year=%.2f\n", 
            b_in, discov_name, comp_name, year);
#endif

/*
if(comp_name[0] != '\0') printf("DEBUG123: %s comp_name=%s -> ", discov_name, comp_name);
*/

// Change some companion names:
  if(!strcmp(comp_name, "ab")) strcpy(comp_name, "AB");
  if(!strcmp(comp_name, "bc")) strcpy(comp_name, "BC");
  if(!strcmp(comp_name, "cd")) strcpy(comp_name, "CD");
  if(!strcmp(comp_name, "ac")) strcpy(comp_name, "AC");
  if(!strcmp(comp_name, "ab-c")) strcpy(comp_name, "AB-C");
/*
if(comp_name[0] != '\0') printf("DEBUG1234/comp_name=%s\n ", comp_name);
*/

        get_data_from_WDS_catalog(WDS_catalog, discov_name, comp_name,
                                  wds_name,
                                  &WdsLastYear, &WdsLastRho, &WdsLastTheta, 
                                  &mag_A, &mag_B, spectral_type, 
                                  &found);
#ifdef DEBUG
printf("(discov_name=%s< wds_name=>%s< stat=%d)\n", 
       discov_name, wds_name, status);
#endif
         if(found != 0) {
           sprintf(b_in, "%s = %s%s & %d & & & & & & & WY=%d WT=%d WR=%3.1f \\\\\n",
                   wds_name, discov_name, comp_name, 
                   (int)year, (int)WdsLastYear, 
                   (int)WdsLastTheta, WdsLastRho);
        } else {
            search_discov_name_in_WDS_catalog(WDS_catalog, discov_name,
                                              comp_name,
                                              wds_name, wds_discov_name,
                                              wds_comp_name, &found1);

          if(found1 != 0) {
          fprintf(stderr, "astrom_add_WDS_from_discov/Error: observations not found for: wds=%s discov=%s comp=%s\n",  
                   wds_name, discov_name, comp_name);
           sprintf(b_in, "%s = %s%s & %d & & & & & & \\\\\n",
                   wds_name, discov_name, comp_name, (int)year);
           } else {
          fprintf(stderr, "astrom_add_WDS_from_discov/Error: wds_name not found for discov=%s\n",  
                   discov_name);
           sprintf(b_in, "%s%s & %d & & & & & & \\\\\n",
                   discov_name, comp_name, (int)year);
           }
        } 
    } /* EOF status == 0 (new object name) */

// Copy current line to file:
        fputs(b_in, fp_out);
   } /* EOF if "&" */
  } /* EOF if fgets() */
} /* EOF while loop */

return(0);
}
/*************************************************************************
* Scan the file and add various parameters 
* WDS number, discoverer's name from WDS catalog
* WY, WT, WR (year, theta and rho of the last observation) 
* are read from WDS catalog
*
* INPUT:
* WDS_cat: wdsweb_summ.txt (used to obtain the last observations)
* ADS_WDS_cross : bsd_ads_wds_cross reference file
*
**************************************************************************/
int astrom_add_discov_and_WDS(FILE *fp_in, FILE *fp_out, char *WDS_catalog, 
                              char *ADS_WDS_cross)
{
double magV, B_V, paral, err_paral, magV_A, magV_B; 
double year, WdsLastYear, WdsLastTheta, WdsLastRho, mag_A, mag_B;
int iline, status;
int found, ads_nber;
char b_in[256], ads_name[64], discov_name[64], wds_name[64], *pc;
char comp_name[64], spectral_type[64]; 

year = 0.;
discov_name[0] = '\0';
ads_name[0] = '\0';
wds_name[0] = '\0';

iline = 0;
while(!feof(fp_in))
{
/* Maximum length for a line will be 170 characters: */
  if(fgets(b_in,170,fp_in))
  {
  iline++;
  b_in[255] = '\0';
/* Remove ^M (Carriage Return) if present: */
  pc = b_in;
  while(*pc) {
  if(*pc == '\r') *pc = ' ';
  pc++;
  }
 
/* Check if input LateX line contains the name of the object:
*/
  status = astrom_decode_new_object_name(b_in, ads_name, discov_name, comp_name,
                                         wds_name, &year);
/* Remove the extra-blanks: */
  jlp_trim_string(ads_name, 64);
  jlp_trim_string(discov_name, 64);
  jlp_trim_string(comp_name, 64);
  jlp_trim_string(wds_name, 64);

  if(status == 0) {

#ifdef DEBUG
    printf("DEBUG: %s From decode: ads_name=%s< discov=%s< comp_name=%s< year=%.2f\n", 
            b_in, ads_name, discov_name, comp_name, year);
#endif
/* Remove the extra-blanks: */
    jlp_trim_string(ads_name, 64);
      if((ads_name[0] != '\0') && (discov_name[0] == '\0')) {
        get_discov_from_ads_wds_crossref(ADS_WDS_cross, ads_name, discov_name, 
                                         comp_name, wds_name);
        } else if((ads_name[0] == '\0') && (discov_name[0] != '\0')) {
        get_ads_from_ads_wds_crossref(ADS_WDS_cross, discov_name, ads_name);
        }

/* Remove the extra-blanks: */
  jlp_trim_string(ads_name, 64);
  jlp_trim_string(discov_name, 64);
  jlp_trim_string(comp_name, 64);
  jlp_trim_string(wds_name, 64);

/*
if(comp_name[0] != '\0') printf("DEBUG123: %s comp_name=%s -> ", discov_name, comp_name);
*/

// Change some companion names:
  if(!strcmp(comp_name, "ab")) strcpy(comp_name, "AB");
  if(!strcmp(comp_name, "bc")) strcpy(comp_name, "BC");
  if(!strcmp(comp_name, "cd")) strcpy(comp_name, "CD");
  if(!strcmp(comp_name, "ac")) strcpy(comp_name, "AC");
  if(!strcmp(comp_name, "ab-c")) strcpy(comp_name, "AB-C");
/*
if(comp_name[0] != '\0') printf("DEBUG1234/comp_name=%s\n ", comp_name);
*/

        get_data_from_WDS_catalog(WDS_catalog, discov_name, comp_name,
                                  wds_name,
                                  &WdsLastYear, &WdsLastRho, &WdsLastTheta, 
                                  &mag_A, &mag_B, spectral_type, 
                                  &found);
#ifdef DEBUG
printf("(ads_name=%s discov_name=%s< wds_name=>%s< stat=%d)\n", 
       ads_name, discov_name, wds_name, status);
#endif
         if(found) {
           sprintf(b_in, "%s = %s%s & %s & %d & & & & & & & WY=%d WT=%d WR=%3.1f \\\\\n",
                   wds_name, discov_name, comp_name, 
                   ads_name, (int)year, (int)WdsLastYear, 
                   (int)WdsLastTheta, WdsLastRho);
        } else {
           sprintf(b_in, "%s = %s%s & %s & %d & & & & & & & \\\\\n",
                   wds_name, discov_name, comp_name, ads_name, (int)year);
        } 
    } /* EOF status == 0 (new object name) */

// Copy current line to file:
        fputs(b_in, fp_out);
  } /* EOF if fgets() */
} /* EOF while loop */

return(0);
}
/*********************************************************************
* Check if current line is compatible with the syntax of a new object
* Example:
* 16564+6502 = STF 2118 AB-C & ADS 10279 & 2004. & & & & & & & \\
* Or:
*  & ADS 1079 CD & 2004. & & & & & & & \\
* Or:
  COU 2118 AC &  & 2004. & & & & & & & \\
* Or:
  & COU 2118 A-BC & 2004. & & & & & & & \\
*
* return 0 if syntax is OK, and that all the corresponding items in that case.
*
* INPUT:
*  b_in: input line to be decoded
*
* OUTPUT:
*  ads_name : ADS name (e.g., ADS 345)
*  discov_name : discoverers' designation (e.g. COU 234)
*  comp_name : companion name (null, AB, AC, A-BC, etc)
*  wds_name : eg 03421-4532
*  year: year of observation as indicated in the "astrom" file 
*********************************************************************/
int astrom_decode_new_object_name(char *b_in, char *ads_name, 
                                  char *discov_name, char *comp_name, 
                                  char *wds_name, double *year)
{
int ncol, icol, status = -1, iverbose = 0;
char *pc, *pc1, buffer[256], buff1[256];
double dval;

ads_name[0] = '\0';
discov_name[0] = '\0';
comp_name[0] = '\0';
wds_name[0] = '\0';
*year = 0.;

/* If comment, exit from here: */
if(b_in[0] == '%' || b_in[1] == '%' || b_in[0] == '\\') return(-1);

/* First check that there are 10 columns: */
pc = b_in;
ncol = 1;
while(*pc) {
  if(*pc == '&') ncol++;
  pc++;
 }
if(ncol != 10) return(-1);

// Read first column
icol = 1;
status = latex_read_svalue(b_in, buffer, icol);
if((status == 0) && (buffer[0] != '\0')) { 
  pc = strchr(buffer, '=');
// Case '01234+4532 = COU 184 AC'
  if(pc != NULL) {
// copy and truncate first part, before '=': 
// ads_name:
    strcpy(buff1, buffer);
    pc1 = strchr(buff1, '=');
    *pc1 = '\0';
    strcpy(wds_name, pc1);
// copy second part, after '=': 
// discov_name:
    pc++;
    if(*pc) strcpy(discov_name, pc);
    jlp_trim_string(discov_name, 32);
// Case 'COU 184 AC'
  } else {
    strcpy(discov_name, buffer);
  }
}

// Decode wds_name:
// Remove first and trailing blanks
jlp_trim_string(wds_name, 32);
if(wds_name[0] != '\0') {
 pc = wds_name;
 while(*pc) if(isdigit(*pc) || *pc == '-' || *pc == '+') pc++; 
 *pc = '\0';
 }

// Read second column
icol = 2;
// Case 'ADS 184 AC'
status = latex_read_svalue(b_in, buffer, icol);
jlp_trim_string(buffer, 32);
if((buffer[0] == 'A') && (buffer[1] == 'D')
   && (buffer[2] == 'S')) {
  strcpy(ads_name, buffer);
// Decode comp_name from ads_name:
// Remove first and trailing blanks (first time)
  jlp_trim_string(ads_name, 32);
  pc = ads_name;
// Case 'ADS 184 AC'
  while(*pc && isalpha(*pc)) pc++;
  while(*pc && (*pc == ' ')) pc++;
  while(*pc && isdigit(*pc)) pc++;
// Case ' AC'
  if(*pc != '\0') {
    strcpy(comp_name, pc);
// Remove first and trailing blanks 
    jlp_trim_string(comp_name, 32);
// Truncate ADS name before comppanion:
    *pc = '\0';
    }
// Case 'COU 184 AC'
  } else {
  strcpy(discov_name, buffer);
  jlp_trim_string(discov_name, 32);
  }

// Decode comp_name from discov name, and truncate discov_name before comp:
if(discov_name[0] != '\0') {
// Case 'COU 184 AC'
// Remove first and trailing blanks (first time)
  jlp_trim_string(discov_name, 32);
  pc = discov_name; 
// Skip alpha characters, ' ', and digits
  while(*pc && isalpha(*pc)) pc++; 
  while(*pc && *pc == ' ') pc++;
  while(*pc && isdigit(*pc)) pc++;
  if(*pc) {
      strcpy(comp_name, pc);
      jlp_trim_string(comp_name, 32);
// Truncate discov_name:
      *pc = '\0';
      jlp_trim_string(discov_name, 32);
      }
  }

/* Then check that there is a date in the 3rd column: */
icol = 3;
*year = 0.;
status = latex_read_dvalue(b_in, &dval, icol, iverbose);
if(status == 0){
  if(dval < 1900. || dval > 2100.){
    status = -2;
    } else {
    *year = dval;
    }
  }

return(status);
}
/**************************************************************************
* Get the discover name from the full ads name (ADS number and companion)
* (e.g., ADS 17180 AB-C will become 'BGH   1')
*
* INPUT:
*  ADS_WDS_cross: name of the ADS-WDS cross reference file
*  ads_name1 : full ads name (ADS number and companion, e.g. 'ADS 17180 AB-C') 
*
* OUTPUT:
*  discov_name1 : discoverer's name of ads_name1
*  comp_name1 : companion name
*  wds_name1 : WDS name corresponding to ads_name1
***************************************************************************/
int get_discov_from_ads_wds_crossref(char *ADS_WDS_cross, char *ads_name1,
                                     char *discov_name1, char *comp_name1,
                                     char *wds_name1)
{
char buffer[256], ads_name0[256], discov_name0[256], comp_name0[256], *pc;
char wds_name0[256], in_comp_name1[256];
int iline, ival, ads_nber0, ads_nber1, status = -1;
FILE *fp_in;

// Decode input ads name:

// ADS name:
 strcpy(buffer, &(ads_name1[3])); 
 ival = sscanf(buffer, "%d", &ads_nber1); 
 if(ival != 1) {
  fprintf(stderr, "read_ads_wds_crossref/error decoding %s input ads_name=%s ival=%d\n",
          buffer, ads_name1, ival); 
  return(-1);
  }
 if((ads_nber1 <= 0) || (ads_nber1 > 17180)) {
  fprintf(stderr, "read_ads_wds_crossref/error decoding %s input ads_name=%s ads_nber=%d\n",
          buffer, ads_name1, ads_nber1); 
  return(-1);
  }

// Companion name:
  strcpy(buffer, &(ads_name1[4]));
// Remove first and trailing blanks:
  jlp_trim_string(buffer, 256);
// First part with digital characters:
  pc = buffer;
  while(*pc) {
   if(isdigit(*pc)) pc++;
   else break;
   }
  strcpy(in_comp_name1, pc);
// Remove first and trailing blanks:
  jlp_trim_string(in_comp_name1, 256);
  if(!strcmp(in_comp_name1, "ab")) strcpy(in_comp_name1, "AB");

if((fp_in = fopen(ADS_WDS_cross, "r")) == NULL) {
  fprintf(stderr, "read_ads_wds_cross_reference/error opening %s\n",
          ADS_WDS_cross); 
  return(-1);
  }

iline = 0;
while(!feof(fp_in))
{
/* Maximum length for a line will be 256 characters: */
  if(fgets(buffer, 256, fp_in))
  {
  iline++;
  buffer[255] = '\0';
  if(buffer[0] != '%') {
// ads_name0 from '    1' to '17180'
    strcpy(ads_name0, &buffer[10]);
    ads_name0[5] ='\0';
    sscanf(ads_name0, "%d", &ads_nber0); 
// 7 characters: 'J   858' or  'HDS3301' for instance
    strcpy(discov_name0, &buffer[20]);
    discov_name0[7] ='\0';
// 4 characters maximum: A-BC, AB-C or AB-D...
    strcpy(comp_name0, &buffer[29]);
    comp_name0[4] ='\0';
// Remove first and trailing blanks:
    jlp_trim_string(comp_name0, 256);
    strcpy(wds_name0, &buffer[38]);
// -----     16949     STI1200                      {1}
    if(wds_name0[0] == ' ') wds_name0[0] = '\0';
    if(ads_nber1 == ads_nber0) {
      if(!strcmp(comp_name0, in_comp_name1)
         || (!strcmp(in_comp_name1, "AB") && (comp_name0[0] == '\0')) 
         || (!strcmp(comp_name0, "AB") && (in_comp_name1[0] == '\0'))) {
        status = 0;
        strcpy(discov_name1, discov_name0);
        strcpy(comp_name1, in_comp_name1);
        strcpy(wds_name1, wds_name0);
        break;
        }
      }
    } 
  } /* EOF if fgets() */
} /* EOF while loop */

fclose(fp_in);

return(status);
}
/**************************************************************************
* Get the ads name from the discover name
* (e.g., 'BGH   1' will become 'ADS 17180 AB-C')
*
* INPUT:
*  ADS_WDS_cross: name of the ADS-WDS cross reference file
*  discov_name1 : discoverer's name
*
* OUTPUT:
*  ads_name1 : ads name (ADS number, e.g. 'ADS 17180') 
***************************************************************************/
int get_ads_from_ads_wds_crossref(char *ADS_WDS_cross, char *discov_name1,
                                  char *ads_name1)
{
char buffer[256], discov_name0[256], in_discov_name1[256], *pc;
int iline, ival, ads_nber0, status = -1;
FILE *fp_in;

strcpy(in_discov_name1, discov_name1);
// Compact input name (i.e., 'BGH   1' becomes 'BGH1') 
// First remove first and trailing blanks:
jlp_trim_string(in_discov_name1, 256);
// Then remove remaining blanks in the middle:
pc = in_discov_name1;
while(*pc && isalpha(*pc)) pc++;
while(*pc && (*pc == ' ')) pc++;
while(*pc && isdigit(*pc)) pc++;
*pc = '\0';

/* DEBUG
printf("DEBUG: searching for discov_name=%s< (or >%s<) in ads_wds_crossref \n",
       discov_name1, in_discov_name1);
*/

if((fp_in = fopen(ADS_WDS_cross, "r")) == NULL) {
  fprintf(stderr, "read_ads_wds_cross_reference/error opening %s\n",
          ADS_WDS_cross); 
  return(-1);
  }

iline = 0;
while(!feof(fp_in))
{
/* Maximum length for a line will be 256 characters: */
  if(fgets(buffer, 256, fp_in))
  {
  iline++;
  buffer[255] = '\0';
  if(buffer[0] != '%') {
// ads_name0 from '    1' to '17180'
    ival = sscanf(&buffer[10], "%d", &ads_nber0);
    if(ival == 1) {
// 7 characters: 'J   858' or  'HDS3301' for instance
      strcpy(discov_name0, &buffer[20]);
      discov_name0[7] ='\0';
// Compact input name (i.e., 'BGH   1' becomes 'BGH1') 
// First remove first and trailing blanks:
      jlp_trim_string(discov_name0, 256);
// Then remove remaining blanks in the middle:
      pc = discov_name0;
      while(*pc && isalpha(*pc)) pc++;
      while(*pc && (*pc == ' ')) pc++;
      while(*pc && isdigit(*pc)) pc++;

// -----     16949     STI1200                      {1}
      if(!strcmp(in_discov_name1, discov_name0)) {
/* DEBUG
printf("DEBUG: discov_name=%d found in ads catalog\n", discov_name0);
*/
        status = 0;
        sprintf(ads_name1, "ADS %d", ads_nber0);
        break;
        }
      } // EOF if ival == 1
    } // EOF buffer[0] != '%'
  } /* EOF if fgets() */
} /* EOF while loop */

fclose(fp_in);

return(status);
}

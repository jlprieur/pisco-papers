/************************************************************************
* "update_PISCO_catalog.c"
* To add the WDS/HIP/CCDM names to all objects of the PISCO catalog containing 
* the list of targets (used by TAV1.EXE)
*
* JLP 
* Version 20/04/2011
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>                  /* exit() */
#include <string.h>
#include <ctype.h>                   /* isprint... */
#include <math.h>
#include <time.h>                    /* date */
#include "jlp_string.h"  /* (in jlplib/jlp_fits/ ) compact_string... */
#include "jlp_catalog_utils.h"       /* Routines used to read catalogs */
#include "WDS_catalog_utils.h"       /* Routines used to read WDS catalog */
#include "HIP_catalog_utils.h"       /* Routines used to read Hipparcos catalog */

#define DEGTORAD   (PI/180.00)

/*
#define DEBUG
#define DEBUG_1
*/

int update_PISCO_catalog_main(char* in_PISCO_catalog, char *WDS_catalog,
                              char *ADS_WDS_cross, char *HIP_HDS_WDS_cross,
                              char *HIC_catalog, char *HIP_catalog,
                              char *out_PISCO_catalog);
int scan_and_update_PISCO_catalog(FILE *fp_in_cat, char *WDS_catalog, 
                                  char *ADS_WDS_cross, char *HIP_HDS_WDS_cross,
                                  char *HIC_catalog, char *HIP_catalog,
                                  FILE *fp_out_cat);
static int update_second_line(char *in_line, FILE *fp_out_cat,
                              char *object_name, char *discov_name, 
                              char *comp_name, char *ADS_name, 
                              double alpha_Pcat, double delta_Pcat, 
                              double equinox_Pcat, 
                              char *WDS_catalog, char *ADS_WDS_cross, 
                              char *HIP_HDS_WDS_cross, char *HIC_catalog, 
                              char *HIP_catalog,
                              double D_tolerance, int *n_objects, int *n_HIP);

int main(int argc, char *argv[])
{
char in_PISCO_catalog[80], out_PISCO_catalog[80]; 
char ADS_WDS_cross[80], HIP_HDS_WDS_cross[80], WDS_catalog[80]; 
char HIC_catalog[80], HIP_catalog[80];
int status;

if(argc != 8) {
  printf("argc = %d\n", argc);
  printf("Syntax: update_PISCO_catalog out_PISCO_catalog in_PISCO_catalog WDS_catalog ADS_WDS_cross HIP_HDS_WDS_cross HIC_catalog HIP_main_catalog\n");
  printf("Example: update_PISCO_catalog test.cat zeiss_doppie.cat wdsweb_summ.txt bds2ads2wds.txt Hipparcos_HIC.dat Hipparcos_main.dat\n");
  return(-1);
}
strcpy(out_PISCO_catalog, argv[1]);
strcpy(in_PISCO_catalog, argv[2]);
strcpy(WDS_catalog, argv[3]);
strcpy(ADS_WDS_cross, argv[4]);
strcpy(HIP_HDS_WDS_cross, argv[5]);
strcpy(HIC_catalog, argv[6]);
strcpy(HIP_catalog, argv[7]);

printf("OK: in_cat=%s out_cat=%s\n", in_PISCO_catalog, out_PISCO_catalog);
printf("OK: WDS=%s ADS_WDS=%s HIP_HDS_WDS=%s HIC=%s HIP=%s\n", 
        WDS_catalog, ADS_WDS_cross, HIP_HDS_WDS_cross, HIC_catalog, HIP_catalog);

/* Call update_PISCO_catalog_main that does the main job: */
status = update_PISCO_catalog_main(in_PISCO_catalog, WDS_catalog, 
                                   ADS_WDS_cross, HIP_HDS_WDS_cross, 
                                   HIC_catalog, HIP_catalog, 
                                   out_PISCO_catalog); 

return(status);
}
/************************************************************************
* update_PISCO_catalog_main
* main routine
*
* INPUT:
* in_PISCO_catalog: name of the input Latex catalog containing the list 
*             of PISCO targets 
* WDS_catalog: name of the WDS catalog
* ADS_WDS_cross: name of the ADS/WDS cross-reference file 
* HIP_HDS_WDS_cross: name of the HIP/WDS cross-reference file
* HIC_catalog: name of the Hipparcos Input Catalog
* HIP_catalog: name of the Hipparcos-Tycho main catalog
* out_PISCO_catalog: name of modified catalog with the list of PISCO targets
*              and the WDS names 
*
*************************************************************************/
int update_PISCO_catalog_main(char* in_PISCO_catalog, char *WDS_catalog,
                              char *ADS_WDS_cross, char *HIP_HDS_WDS_cross,
                              char *HIC_catalog, char *HIP_catalog,
                              char *out_PISCO_catalog)
{
FILE *fp_in_cat, *fp_out_cat; 
time_t tti = time(NULL);

/* Open input file containing the original Latex PISCO catalog */
if((fp_in_cat = fopen(in_PISCO_catalog, "r")) == NULL) {
   fprintf(stderr, "update_PISCO_catalog_main/Fatal error opening input PISCO catalog: %s\n",
           in_PISCO_catalog);
   return(-1);
  }

/* Open output file containing the modified Latex PISCO catalog */
if((fp_out_cat = fopen(out_PISCO_catalog, "w")) == NULL) {
   fprintf(stderr, "update_PISCO_catalog_main/Fatal error opening output PISCO catalog: %s\n",
           out_PISCO_catalog);
   fclose(fp_in_cat);
   return(-1);
  }

/* Header of the output Latex PISCO catalog: */
fprintf(fp_out_cat, "%% From %s, modified by update_PISCO_catalog.c (%s)\n", 
        in_PISCO_catalog, ctime(&tti));

/* Scan the input file and compute the update_PISCO_catalog 
*/
scan_and_update_PISCO_catalog(fp_in_cat, WDS_catalog, ADS_WDS_cross, 
                              HIP_HDS_WDS_cross, HIC_catalog, HIP_catalog,
                              fp_out_cat);

/* Close opened files:
*/
fclose(fp_in_cat);
fclose(fp_out_cat);
return(0);
}
/***************************************************************************
* Scan the input PISCO catalog and insert the WDS name
* Example:
* ADS 326AB (3.40, 340) &  0 23 59.8 &-03 28 31 & 2000.0
* & A: 7.4, B: 9.3 A2 & BU   488,   & \cr
* will become:
* ADS 326AB (3.40, 340) &  0 23 59.8 &-03 28 31 & 2000.0
* & A: 7.4, B: 9.3 A2 & BU   488AB, & WDS0000+000 HIP234 PI_HIP=0.34(+/-0.02) V_HIP=7.6 B-V_HIP=0.34 \cr
*
* INPUT:
* fp_in_cat: pointer to the input Latex PISCO catalog
* WDS_catalog: name of the WDS catalog
* ADS_WDS_cross: name of the ADS/WDS cross-reference file 
* HIP_HDS_WDS_cross: name of the HIP/WDS cross-reference file
* HIC_catalog: name of the Hipparcos Input catalog
* HIP_catalog: name of the Hipparcos-Tycho main catalog
* fp_out_cat: pointer to the output Latex PISCO catalog
***************************************************************************/
int scan_and_update_PISCO_catalog(FILE *fp_in_cat, char *WDS_catalog, 
                                  char *ADS_WDS_cross, char *HIP_HDS_WDS_cross,
                                  char *HIC_catalog, char *HIP_catalog,
                                  FILE *fp_out_cat)
{
char in_line[128], spectral_type[40];
char object_name[40], comp_name[20], ADS_name[20], discov_name[20];
double alpha_Pcat, delta_Pcat, equinox_Pcat, D_tolerance;
double magV_A, magV_B;
int status, first_line_is_read, iline, n_objects, n_HIP, n_focusing;
int has_an_orbit;

/* Tolerance of the coordinates in degrees 
* used for searching for Hipparcos names.
* 0.1 arcminute is a good value: */
D_tolerance = 0.1/60.;

/* Scan the input Latex PISCO catalog: */
first_line_is_read = 0;
iline = 0;
n_objects = 0;
n_HIP = 0;
/* Stars used to control Risley prisms and focusing: */
n_focusing = 0;
/* DEBUG: */
#ifdef DEBUG
printf("DEBUG/WARNING: limitation to the first 100 lines! \n");
while(!feof(fp_in_cat) && iline < 100) {
#else
while(!feof(fp_in_cat)) {
#endif
  if(fgets(in_line, 120, fp_in_cat)) {
    iline++;
/* Check if it is a first line or a second line */
    if(in_line[0] != '%') {
/* Copy the first line */
      if(in_line[0] != '&') {
         status = PISCO_catalog_read_line1(in_line, object_name, comp_name, 
                                           ADS_name, &alpha_Pcat, &delta_Pcat, 
                                           &equinox_Pcat);
         if(status) {
           fprintf(stderr, "Fatal error/Bad object syntax in line #%d\n", iline);
           exit(-1);
           }
#ifdef DEBUG_1
if(comp_name[0] != '\0') printf("iline=%d object name=>%s< comp_name=>%s<\n", 
                iline, object_name, comp_name);
#endif
         first_line_is_read = 1;
         cleanup_string(in_line, 128);
         fprintf(fp_out_cat, "%s\n", in_line);
/* Decode and modify the second line */
      } else {
         if(!first_line_is_read) {
            fprintf(stderr, "Fatal error/Bad syntax in line #%d\n", iline);
            exit(-1);
            }
         first_line_is_read = 0;
         PISCO_catalog_read_line2_discov(in_line, discov_name);
         status = PISCO_catalog_read_line2(in_line, &magV_A, &magV_B,
                                           spectral_type, &has_an_orbit);
         if (status) {
            fprintf(stderr, "Fatal error/Bad syntax in line #%d\n", iline);
            exit(-1);
            } 
         if(discov_name[0] == '\0') {
            n_focusing++;
#ifdef DEBUG_1
            printf("Is %s a star used fo focusing the telescope ?\n", 
                    object_name);
#endif
/* Simply copy the input to the output when not a discovered binary: */
            cleanup_string(in_line, 128);
            fprintf(fp_out_cat, "%s\n", in_line);
            } else if (status == -1) {
            fprintf(stderr, "Fatal error/Bad syntax in line #%d\n", iline);
            exit(-1);
            } else {
/* Update the second line with extra parameters: */
            update_second_line(in_line, fp_out_cat, object_name, discov_name, 
                               comp_name, ADS_name, alpha_Pcat, delta_Pcat, 
                               equinox_Pcat,
                               WDS_catalog, ADS_WDS_cross, HIP_HDS_WDS_cross,
                               HIC_catalog, HIP_catalog, D_tolerance, 
                               &n_objects, &n_HIP);
            }
#ifdef DEBUG_1
if(comp_name[0] != '\0') printf("iline=%d discov name=>%s<\n", iline, 
                                 discov_name);
#endif
/* & A: 7.4, B: 9.3 A2 & BU   488AB, & WDS0000+000 HIP234 PI_HIP=0.34(+/-0.02) V_HIP=7.6 B-V_HIP=0.34 \cr
*/
      }
/* If commented line, simply copy this line to the ouput file: */
    } else {
      cleanup_string(in_line, 128);
      fprintf(fp_out_cat, "%s\n", in_line);
    }
  }
}

printf("End of process: %d lines read (%d objects modified, %d focusing stars)\n", 
        iline, n_objects, n_focusing);

printf("%d Hipparcos objects identified among the known binaries (coordinate-tolerance = %.3f arcmin)\n", 
        n_HIP, D_tolerance * 60.);

return(0);
}
/*************************************************************************
* PISCO catalog (used by TAV1.EXE)
* Modify the second line relative to an object
* Example:
* & A: 7.0, B: 8.2 G5 & BU   862, orbit & \cr
* will become:
* & A: 7.0, B: 8.2 G5 & BU   862, orbit & WDS0000+000 HIP00 PI_HIP=3.5(+/-2)mas \cr
*
* INPUT:
* in_line: current (second) line of the input PISCO catalog to be processed
* fp_out_cat: pointer to the ouput PISCO catalog 
* object_name, discov_name, comp_name: names of the object, of the discoverer's
*                                      and of the companion as read in the
*                                      first line of the input PISCO catalog.
* alpha_Pcat, delta_Pcat, equinox_Pcat : coordinates of the object from the
*                         first line of the input PISCO catalog.
* WDS_catalog: name of the input WDS catalog
* ADS_WDS_cross: name of the cross-catalog between ADS and WDS names
* HIP_HDS_WDS_cross: name of the cross-catalog between HIP and WDS names
*
* INPUT/OUTPUT
* n_objects: counter of the objects successfuly processed
* n_HIP: counter of the Hipparcos objects
*
*************************************************************************/
static int update_second_line(char *in_line, FILE *fp_out_cat,
                              char *object_name, char *discov_name, 
                              char *comp_name, char *ADS_name, 
                              double alpha_Pcat, double delta_Pcat, 
                              double equinox_Pcat, 
                              char *WDS_catalog, char *ADS_WDS_cross, 
                              char *HIP_HDS_WDS_cross, char *HIC_catalog, 
                              char *HIP_catalog,
                              double D_tolerance, int *n_objects, int *n_HIP)
{
char WDS_name[40], HIP_name_cross[40], HIP_name[40], CCDM_name[40]; 
char *pc, buffer[128];
double paral, err_paral, V_mag, B_V_index;
int found, found_in_cross, found_in_Hip_cat, is_OK;

/* Search for WDS number in WDS catalog using discov_name */
search_discov_name_in_WDS_catalog(WDS_catalog, discov_name, WDS_name, &found);

/* If not found, copy the line without changes */
if(!found) {
  fprintf(stderr, "update_second_line/WARNING: %s (=%s) not found in WDS catalog\n",
          object_name, discov_name);
  cleanup_string(in_line, 128);
/* Copy input line to output catalog without any modification: */
  fprintf(fp_out_cat, "%s\n", in_line);
  return(0);
  } 

/* Check the consistency between the WDS name and the coordinates contained
in PISCO catalog:
*/
  check_consistency_coord_WDSname(alpha_Pcat, delta_Pcat, equinox_Pcat, 
                                  object_name, WDS_name, &is_OK);
/* Check the consistency between the ADS and the WDS names 
*/
  check_consistency_ADSname(ADS_name, WDS_name, discov_name, ADS_WDS_cross, 
                            &is_OK);

/* Search for WDS number in HIPHIP//HDS/WDS cross-file 
* (only for debugging, since there are very few objects in this list)
*/
  HIP_name_from_HIP_HDS_WDS_cross(WDS_name, HIP_HDS_WDS_cross, HIP_name_cross, 
                                  &found_in_cross);
  if(found_in_cross) {
    printf("%s=%s=%s was found in HIP/HDS/WDS cross-file (=HIP %s)\n", 
            object_name, discov_name, WDS_name, HIP_name_cross);
   } 

/* DEBUG: search in Hipparcos catalog */
  search_object_in_HIC_catalog(HIC_catalog, alpha_Pcat, delta_Pcat, 
                               equinox_Pcat, HIP_name, CCDM_name,
                               &V_mag, &B_V_index, D_tolerance,
                               &found_in_Hip_cat);
  if(found_in_Hip_cat) {
    printf("%s=%s=%s was found in Hipparcos catalog (=HIP%s=CCDM%s)\n", 
            object_name, discov_name, WDS_name, HIP_name, CCDM_name);
    printf("V=%.3f B-V=%.3f\n", V_mag, B_V_index);
    (*n_HIP)++;
    read_data_in_HIP_catalog(HIP_catalog, HIP_name, &paral, &err_paral, 
                             &found); 
    if(found) {
       printf("Paral=%.2f+/-%.2f\n", paral, err_paral);
      } else {
       fprintf(stderr, "Fatal error: HIP=%s not in Hipparcos main catalog!\n", HIP_name);
       exit(-1);
      }
   }  

if(found_in_cross && strcmp(HIP_name_cross, HIP_name)) {
  fprintf(stderr, "Fatal error: inconsistency: object=%s = HIP %s in cross-reference file\n", 
          object_name, HIP_name_cross);
  fprintf(stderr, "and coordinates correspond to HIP %s in Hipparcos catalog!\n", 
          HIP_name);
  exit(-1);
  }

/* Copy the modified input line to ouput catalog: */
  strncpy(buffer, in_line, 128);
  pc = buffer;
  while(*pc && strncmp(pc, "\\cr", 3)) pc++;
  *pc = '\0';
  fprintf(fp_out_cat, "%sWDS%s ", buffer, WDS_name);
  if(*HIP_name) {
    fprintf(fp_out_cat, "HIP%s ", HIP_name);
    compact_string(CCDM_name, 20);
    if(*CCDM_name) fprintf(fp_out_cat, "CCDM%s ", CCDM_name);
    if(V_mag != 100.) fprintf(fp_out_cat, "V=%.2f ", V_mag);
    if(B_V_index != 100.) fprintf(fp_out_cat, "B-V=%.2f ", B_V_index);
/* Parallax in mas: */
    if(paral > 0.) {
      if(err_paral > 0.) 
        fprintf(fp_out_cat, "PI=%.2f(%.2f) ", paral, err_paral);
      else
        fprintf(fp_out_cat, "PI=%d ", (int)paral);
      }
   }
  fprintf(fp_out_cat, "\\cr\n");
  (*n_objects)++;

return(0);
}

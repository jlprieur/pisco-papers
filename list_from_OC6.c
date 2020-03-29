/************************************************************************
* "list_from_OC6.c"
*
# Extract the orbital parameters from OC6 for the objects contained
# in the LaTeX calibrated table of measurements
* (in order to compute the residuals with residuals_1.c)
*
* JLP 
* Version 16/11/2009
*************************************************************************/
#include "jlp_catalog_utils.h"
#include "residuals_utils.h"

/*
#define DEBUG
#define DEBUG_1
*/

static int extract_orbits_from_OC6(FILE *fp_out, FILE *fp_calib, 
                                   char *OC6_fname, int is_master_file,
                                   int *n_required_orbits,
                                   int *n_found_orbits, int norbits_per_object,
                                   int search_for_all);

int main(int argc, char *argv[])
{
char calib_fname[80], OC6_fname[80], out_orbit_list[80];
FILE *fp_out, *fp_calib;
int n_required_orbits, n_found_orbits, norbits_per_object, search_for_all;
/* is_master_file: flag set to 1 if master file ("orb6.master", 
*                           to 0 if OC6 file ("orb6orbits.txt")
*/
int is_master_file;


if(argc == 7) {
  if(*argv[6]) argc = 7;
  else if(*argv[5]) argc = 6;
  else if(*argv[4]) argc = 5;
  else if(*argv[3]) argc = 4;
  else if(*argv[2]) argc = 3;
  else if(*argv[1]) argc = 2;
  else argc = 1;
}
if(argc != 6 && argc != 7) {
  printf("Syntax: list_from_OC6 calibrated_latex_table OC6_catalog out_orbit_list is_master_file search_for_all_objects [norbits_per_object]\n");
  printf("if search_for_all_objects = 1 will search for an orbit\n");
  printf("       even when orbit=0 in last column of tab_calib.tex\n");
  printf("norbits_per_object = nber of orbits per object (default=0=all)\n");
  printf("(0=all orbits, 1=last orbit, 2=last two orbits, etc)\n");
/* calibrated Latex file, OC6, output list, flag: number of orbits 
* (0=all 1=last 2=last two orbits)
*/
  return(-1);
}
/* Calibrated latex table: */
strcpy(calib_fname, argv[1]);
/* File name of the Sixth Orbit Catalog: */
strcpy(OC6_fname, argv[2]);
/* Output orbit list (selection of all the orbits necessary to compute
* the residuals): */
strcpy(out_orbit_list, argv[3]);
sscanf(argv[4], "%d", &is_master_file); 
sscanf(argv[5], "%d", &search_for_all); 
if(argc == 7) {
 sscanf(argv[6], "%d", &norbits_per_object); 
 } else { 
 norbits_per_object = 0;
 }

/* is_master_file: flag set to 1 if master file ("orb6.master", 
*                           to 0 if OC6 file ("orb6orbits.txt")
*/

#ifdef DEBUG
printf("OK: calib_fname=%s out_orbit_list=%s norbits_per_object=%d \n",
       calib_fname, out_orbit_list, norbits_per_object);
printf("search_for_all=%d\n", search_for_all);
printf("OK: OC6=%s (is_master_file=%d)\n", OC6_fname, is_master_file);
#endif

fp_out = NULL;
fp_calib = NULL;

/* Open output orbit list: */
if((fp_out = fopen(out_orbit_list, "w")) == NULL) {
   fprintf(stderr, "list_from_OC6/Fatal error opening output file: %s\n",
           out_orbit_list);
    return(-1);
  }

/* Open input latex calibrated file: */
if((fp_calib = fopen(calib_fname, "r")) == NULL) {
   fprintf(stderr, "list_from_OC6/Fatal error opening input LaTeX file: %s\n",
           calib_fname);
    return(-1);
  }

/* Scan the input LaTeX file and retrieve the orbits from OC6 
*/
extract_orbits_from_OC6(fp_out, fp_calib, OC6_fname, is_master_file,
                        &n_required_orbits, &n_found_orbits, norbits_per_object,
                        search_for_all);

/* Epilog for orbit file: */
fprintf(fp_out, "# %d series of orbits were required in %s and %d were found in OC6\n", 
                 n_required_orbits, calib_fname, n_found_orbits);
printf(" %d series of orbits were required in %s and %d were found in OC6...\n", 
                 n_required_orbits, calib_fname, n_found_orbits);

/* Close opened files:
*/
if(fp_calib != NULL) fclose(fp_calib);
if(fp_out != NULL) fclose(fp_out);
return(0);
}

/*************************************************************************
* Scan the calibrated table and when "orb" is found,
* look for the corresponding orbital parameters in the OC6 catalog 
*
* INPUT:
* fp_out: pointer to the output orbit list file
* fp_calib: pointer to the input LaTeX calibrated file 
* OC6_fname: name of the input Sixth Orbit Catalog
* is_master_file: flag set to 1 if master file ("orb6.master", 
*                          to 0 if OC6 file ("orb6orbits.txt")
*
* OUTPUT:
* n_required_orbits: number of orbits required in Latex calibrated file
* n_found_orbits: number of corresponding orbits actually found in OC6 
**************************************************************************/
static int extract_orbits_from_OC6(FILE *fp_out, FILE *fp_calib, 
                                   char *OC6_fname, int is_master_file,
                                   int *n_required_orbits,
                                   int *n_found_orbits, int norbits_per_object,
                                   int search_for_all)
{
char in_line[256], ads_name[40], discov_name[40], comp_name[40], wds_name[40]; 
char buffer[128], previous_ads_name[40], previous_comp_name[40];
char previous_discov_name[40];
int iline, status, orbit, verbose_if_error = 1, found, candidate_found;

*n_required_orbits = 0;
*n_found_orbits = 0;

previous_ads_name[0] = '\0';
previous_discov_name[0] = '\0';
previous_comp_name[0] = '\0';
iline = 0;
/* Scan the input LaTeX calibrated file */
while(!feof(fp_calib)) {
  if(fgets(in_line, 256, fp_calib)) {
    iline++;
/* Good lines start with a digit (WDS names...) */
    if(in_line[0] != '%' && (isdigit(in_line[0]) 
        || !strncmp(in_line, "\\idem", 5))) {
/* Remove the end of line '\n' from input line: */
      cleanup_string(in_line, 256);

/* Read orbit flag in column 11: */
      status = latex_get_column_item(in_line, buffer, 11, verbose_if_error);
      if(status || (sscanf(buffer,"%d", &orbit) != 1)) {
         fprintf(stderr,"Fatal error reading orbit: %s (line=%d)\n", 
                 in_line, iline);
         exit(-1);
         }

/* If orbital binary : */
     if(orbit || search_for_all) {
/* Retrieve the residuals for this object and epoch in the residual table: */
/* Read object name and companion name from columns 2 and 3: */
      if(isdigit(in_line[0])) {
         read_object_name_from_CALIB_line(in_line, wds_name, discov_name, 
                                          comp_name, ads_name);
         compact_string(comp_name,40);
         (*n_required_orbits)++;
#ifdef DEBUG
         printf("Orbit needed for ads=>%s< disc=>%s< comp=>%s<\n", 
               ads_name, discov_name, comp_name);
#endif
/* If some orbits are found, copy the full lines to the output list file: */
       line_extraction_from_OC6_catalog(OC6_fname, is_master_file,
                                        ads_name, discov_name, comp_name,
                                        fp_out, &found, &candidate_found,
                                        norbits_per_object);
       if(found) {
        (*n_found_orbits)++;
       } else if (candidate_found) {
         printf("Candidate orbit found for ads=>%s< disc=>%s< comp=>%s<\n", 
               ads_name, discov_name, comp_name);
       } else if(!search_for_all) { 
         printf("*** WARNING: orbit not found for a=>%s< d=>%s< c=>%s<!\n", 
               ads_name, discov_name, comp_name);
       }
      } else {
         strcpy(ads_name, previous_ads_name);
         strcpy(discov_name, previous_discov_name);
         strcpy(comp_name, previous_comp_name);
      }
     } /* EOF orbit case */

    }/* EOF if !isdigit ... */
  } /* EOF if fgets */ 
/* Load object name to handle the case of "idem" (i.e. multiple measurements
* of the same object, without repeating the object name in cols. 1 2 3)*/
 strcpy(previous_ads_name, ads_name);
 strcpy(previous_discov_name, discov_name);
 strcpy(previous_comp_name, comp_name);
 } /* EOF while ... */
printf("extract_orbits_from_OC6: %d lines sucessfully processed in LaTeX file\n", 
        iline);
return(0);
}

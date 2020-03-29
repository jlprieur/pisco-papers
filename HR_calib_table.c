/************************************************************************
* "HR_calib_table.c"
*
* To create a file with photometric and miscelaneous data
* in order to prepare various plots (HR diagram or statistics)
* Can create a Hertzprung-Russel diagram of the objects contained 
* in the calibrated table with the list of targets zeiss_doppie_new.cat
* (that is used by TAV1.EXE)
*
* JLP 
* Version 03/05/2009
*************************************************************************/
#include "jlp_catalog_utils.h"

#define DEBUG
#define DEBUG_1
/*
*/

static int create_HR_of_calib_table(FILE *fp_calib, char *PISCO_catalog_fname, 
                                    FILE *fp_out, double paral_rel_error0);

int main(int argc, char *argv[])
{
char calib_fname[80], PISCO_catalog_fname[80], out_fname[80];
FILE *fp_calib, *fp_out;
time_t ttime = time(NULL);
double paral_rel_error0;

if(argc != 4 && argc != 5) {
  printf("Syntax: HR_calib_table calibrated_table PISCO_catalog out_table [relative_error_of_parallax]\n");
  printf("Enter a negative relative error to output photometric data of all the objects\n");
  printf("Enter a positive relative error to output absolute magnitudes of objects measured by Hipparcos\n");
  return(-1);
}
strcpy(calib_fname, argv[1]);
strcpy(PISCO_catalog_fname, argv[2]);
strcpy(out_fname, argv[3]);
if(argc == 5) {
 sscanf(argv[4], "%lf", &paral_rel_error0);
 } else {
 paral_rel_error0 = 0.8;
 }

printf("OK: calib=%s catalog=%s output=%s \n", calib_fname, 
        PISCO_catalog_fname, out_fname); 

/* Open input calibrated table: */
if((fp_calib = fopen(calib_fname, "r")) == NULL) {
   fprintf(stderr, "HR_calib_table/Fatal error opening calib. table %s\n",
           calib_fname);
    return(-1);
  }

/* Open output file: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "HR_calib_table/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the output photometric file: */
if(paral_rel_error0 < 0) {
  fprintf(fp_out, "%% Photometric data of %s using %s\n%% Created on %s",
        calib_fname, PISCO_catalog_fname, ctime(&ttime));
  fprintf(fp_out, "%% magV magV_A B-V Delta_mag discov_name\n");
/* Header of the output HR file: */
  } else {
  fprintf(fp_out, "%% HR diagram of %s using %s (paral_rel_error0=%.3f) \n%% Created on %s",
        calib_fname, PISCO_catalog_fname, paral_rel_error0, ctime(&ttime));
  fprintf(fp_out, "%% magV magV_A B-V Delta_mag MagV MagV_A discov_name\n");
  }

/* Scan the input calibrated table and generate output file 
* (compute M_V, B-V, and build a HR diagram if paral_rel_error0 > 0)
*/
create_HR_of_calib_table(fp_calib, PISCO_catalog_fname, fp_out, 
                         paral_rel_error0);

/* Close opened files:
*/
fclose(fp_calib);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the calibrated table 
*
* INPUT:
* fp_calib: pointer to the input file containing the calibrated table 
* PISCO_catalog_name: name of the PISCO catalog ("zeiss_doppie.cat") used
*               for retrieving photometric data 
* fp_out: pointer to the output Latex file with the full table 
*
*************************************************************************/
static int create_HR_of_calib_table(FILE *fp_calib, char *PISCO_catalog_fname, 
                                    FILE *fp_out, double paral_rel_error0) 
{
double magV, B_V, paral, err_paral, magV_A, magV_B, magV_abs, magV_A_abs;
int iline, n_objects, status, has_an_orbit;
char in_line[256], wds_name[40]; 
char ads_name[40], object_name[40], discov_name[40], comp_name[40]; 
char spectral_type[40], ads_name2[40], wds_name2[40], discov_name2[40];

iline = 0;
n_objects = 0;
while(!feof(fp_calib)) {
  if(fgets(in_line, 256, fp_calib)) {
    iline++;
/* Good lines start with a digit (WDS names...) */
    if(in_line[0] != '%' && (isdigit(in_line[0]) 
        || !strncmp(in_line, "\\idem", 5))) {

/* Select only one entry for each object (not "\idem", etc) */
      if(isdigit(in_line[0])) {
         read_object_name_from_CALIB_line(in_line, wds_name, discov_name,
                                          comp_name, ads_name);
/* Generate the object name: */
        if(*ads_name)
          strcpy(object_name, ads_name);
        else
          strcpy(object_name, discov_name);
 
/* Look for object_name in file PISCO_catalog_name ("zeiss_doppie.cat"), 
* and determine values of: alpha, delat, coord_equinox,
*                          has_an_orbit, discov_name0, comp_name0.
*/
      status = get_data_from_PISCO_catalog(PISCO_catalog_fname, 
                       object_name, &magV, &B_V, &paral, &err_paral,
                       &magV_A, &magV_B, spectral_type, discov_name2, 
                       ads_name2, wds_name2, &has_an_orbit);

/* JLP2010: possibility of retrieving ADS name from PISCO catalog: */
        if(ads_name[0] == '\0') strcpy(ads_name, ads_name2);

         if(magV != 100.  && magV_A != 100.
            && magV_B != 100. && B_V != 100.) {
/* Photometric data: */
            if(paral_rel_error0 < 0) {
              compact_string(discov_name, 20);
              fprintf(fp_out, "%.3f %.3f %.3f %.3f %s %s\n",
                      magV, magV_A, B_V, ABS(magV_A - magV_B),
                      discov_name, spectral_type);
              n_objects++;
/* HR diagram: */
            } else { 
              if(paral > 0.) {
                if(err_paral/paral < paral_rel_error0) {
/* MV = mV - 2.5 log10 (d/d_10)^2
*  MV = mV - 5 log10 (d) + 5 log10(d_10)
*  MV = mV - 5 log10(d) + 5
* paral = 1/d
* MV = mV + 5 log10(paral) + 5
/* Parallax in mas: */
                magV_abs = magV + 5. * log10(paral * 0.001) + 5.;
                magV_A_abs = magV_A + 5. * log10(paral * 0.001) + 5.;
                compact_string(discov_name, 20);
                fprintf(fp_out, "%.3f %.3f %.3f %.3f %.3f %.3f %s %s\n",
                        magV, magV_A, B_V, ABS(magV_A - magV_B), magV_abs,
                        magV_A_abs, discov_name, spectral_type);
                n_objects++;
                }
              }
            }
           } /* EOF magv != 100 */
      }
    }/* EOF if !isdigit ... */
  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("HR_curve_of_calib_table: %d lines sucessfully read (%d objects)\n", 
        iline, n_objects);
fprintf(fp_out, "%% HR_curve_of_calib_table: %d lines sucessfully read (%d objects)\n", 
        iline, n_objects);
return(0);
}

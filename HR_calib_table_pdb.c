/************************************************************************
* "HR_calib_table_pdb.c"
*
* To create a file with photometric and miscelaneous data
* in order to prepare various plots (HR diagram or statistics)
* Can create a Hertzprung-Russel diagram of the objects contained 
* in the calibrated table with the WDS catalog 
*
* JLP 
* Version 03/05/2009
*************************************************************************/
#include "jlp_catalog_utils.h"

#define DEBUG
#define DEBUG_1
/*
*/

static int create_HR_of_calib_table_pdb(FILE *fp_calib, char *WDS_catalog, 
                                        char *HIC_catalog,
                                        char *HIP_catalog, FILE *fp_out, 
                                        double paral_rel_error0);

int main(int argc, char *argv[])
{
char calib_fname[80], out_fname[80];
char WDS_catalog[80], HIC_catalog[80], HIP_catalog[80];
FILE *fp_calib, *fp_out;
time_t ttime = time(NULL);
double paral_rel_error0;
int status;

if(argc != 6 && argc != 7) {
  printf("Syntax: HR_calib_table_pdb calib_table WDS_cat HIC_catalog HIP_main_catalog out_table [relative_error_of_parallax]\n");
  printf("Enter a negative relative error to output photometric data of all the objects\n");
  printf("Enter a positive relative error to output absolute magnitudes of objects measured by Hipparcos\n");
  printf("Example: HR_calib_table_pdb tab_calib.tex wdsweb_summ.txt wds2hds2hip.txt Hipparcos_HIC.dat Hipparcos_main.dat Photom2a.dat -1\n");
  return(-1);
}
strcpy(calib_fname, argv[1]);
strcpy(WDS_catalog, argv[2]);
strcpy(HIC_catalog, argv[3]);
strcpy(HIP_catalog, argv[4]);
strcpy(out_fname, argv[5]);
if(argc == 7) {
 sscanf(argv[6], "%lf", &paral_rel_error0);
 } else {
 paral_rel_error0 = 0.8;
 }

printf("OK: calib=%s output=%s \n", calib_fname, out_fname); 
printf("OK: WDS=%s HIC=%s HIP=%s\n", WDS_catalog, HIC_catalog, HIP_catalog);

/* Open input calibrated table: */
if((fp_calib = fopen(calib_fname, "r")) == NULL) {
   fprintf(stderr, "HR_calib_table_pdb/Fatal error opening calib. table %s\n",
           calib_fname);
    return(-1);
  }

/* Open output file: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "HR_calib_table_pdb/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the output photometric file: */
if(paral_rel_error0 < 0) {
  fprintf(fp_out, "%% Photometric data of %s\n%% Created on %s",
          calib_fname, ctime(&ttime));
  fprintf(fp_out, "%% magV magV_A B-V Delta_mag rho discov_name spectral_type\n");
/* Header of the output HR file: */
  } else {
  fprintf(fp_out, "%% HR diagram of %s (paral_rel_error0=%.3f) \n%% Created on %s",
        calib_fname, paral_rel_error0, ctime(&ttime));
  fprintf(fp_out, "%% magV magV_A B-V Delta_mag MagV MagV_A rho discov_name spectral_type\n");
  }

/* Scan the input calibrated table and generate output file 
* (compute M_V, B-V, and build a HR diagram if paral_rel_error0 > 0)
*/
status = create_HR_of_calib_table_pdb(fp_calib, WDS_catalog,
                                      HIC_catalog, HIP_catalog, fp_out, 
                                      paral_rel_error0);

/* Close opened files:
*/
fclose(fp_calib);
fclose(fp_out);
return(status);
}
/************************************************************************
* Scan the calibrated table 
*
* INPUT:
* fp_calib: pointer to the input file containing the calibrated table 
* fp_out: pointer to the output Latex file with the full table 
*
*************************************************************************/
static int create_HR_of_calib_table_pdb(FILE *fp_calib, char *WDS_catalog, 
                                        char *HIC_catalog, char *HIP_catalog, 
                                        FILE *fp_out, double paral_rel_error0) 
{
double magV, B_V_index, paral, err_paral, magV_A, magV_B, magV_abs, magV_A_abs;
double flux_A, flux_B;
double epoch0, rho0, err_rho0, theta0, err_theta0;
int iline, status, has_an_orbit, found_in_WDS, found_in_Hip_cat;
int n_wds_objects, n_hip_objects, n_hip_good_objects;
char in_line[256], ads_name[64], wds_name[64]; 
char object_name[64], discov_name[64], comp_name[64], spectral_type[64];
char old_discov_name[64];

iline = 0;
n_wds_objects = 0;
n_hip_objects = 0;
n_hip_good_objects = 0;
old_discov_name[0] = '\0';

while(!feof(fp_calib)) {
  if(fgets(in_line, 256, fp_calib)) {
    iline++;
/* Good lines start with a digit (WDS names...) */
    if(in_line[0] != '%' && (isdigit(in_line[0]) 
        || !strncmp(in_line, "\\idem", 5))) {

/* Select only one entry for each object (not "\idem", etc) */
// NB: in nw format (after 2015,...) the name of the object is always
// duplicated for all measurements (no \idem ...)
      if(isdigit(in_line[0])) {

      read_object_name_from_CALIB_line(in_line, wds_name, discov_name,
                                       comp_name, ads_name);
/* Remove the extra-blanks: */
      trim_string(wds_name, 64);
      trim_string(ads_name, 64);
      trim_string(discov_name, 64);
      trim_string(comp_name, 64);

// New: after 2018, now also read the measurements
// to be able to plot Delta mag versus rho measurements...
      status = read_measures_from_CALIB_line(in_line, &epoch0, &rho0, 
                                             &err_rho0, &theta0, &err_theta0);

// Check if it is a new object and that it was resolved:
      if(strcmp(old_discov_name, discov_name) && (status == 0)) {


       strcpy(old_discov_name, discov_name);

// Get magV, B_V_index, parall, err_paral from Hipparcos
// Get magV_A, magV_B, spectral_type from WDS 
       get_data_from_WDS_and_HIP_catalogs(WDS_catalog, HIC_catalog, HIP_catalog,
                                          discov_name, wds_name, &magV, 
                                          &B_V_index, &paral, &err_paral, 
                                          &magV_A, &magV_B, spectral_type, 
                                          &found_in_WDS, &found_in_Hip_cat);

       if(found_in_WDS) n_wds_objects++;
       if(found_in_Hip_cat) n_hip_objects++;

/* Photometric data: */
         if(found_in_WDS && (paral_rel_error0 < 0)) {
              compact_string(discov_name, 20);
            if(found_in_Hip_cat) {
              fprintf(fp_out, "%.3f %.3f %.3f %.3f %.3f %s %s\n",
                      magV, magV_A, B_V_index, ABS(magV_A - magV_B),
                      rho0, discov_name, spectral_type);
// If not found in Hipparcos, use WDS data only:
            } else {
              B_V_index = 100.;
              flux_A = pow(10., -0.4 * magV_A);
              flux_B = pow(10., -0.4 * magV_B);
              magV = -2.5 * log10( flux_A + flux_B);
              fprintf(fp_out, "%.3f %.3f %.3f %.3f %.3f %s %s\n",
                      magV, magV_A, B_V_index, ABS(magV_A - magV_B),
                      rho0, discov_name, spectral_type);
            }
/* HR diagram: */
         } else if(found_in_Hip_cat && (paral_rel_error0 > 0)) {
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
                        magV, magV_A, B_V_index, ABS(magV_A - magV_B), 
                        magV_abs, magV_A_abs, discov_name, spectral_type);
                n_hip_good_objects++;
                }
              } // If paral > 0
           } /* EOF found in Hip_cat and paral_rel_error > 0 */
      } // If discov_name != old_discov_name
     }/* EOF if !isdigit ... */
    }// EOF if inline ...
  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("HR_curve_of_calib_table_pdb: %d lines sucessfully read (%d wds_objects with %d in Hipparcos catalog and %d with parallax)\n", 
        iline, n_wds_objects, n_hip_objects, n_hip_good_objects);
fprintf(fp_out, "%% HR_curve_of_calib_table_pdb: %d lines sucessfully read (%d wds_objects with %d in Hipparcos catalog and %d with parallax)\n", 
        iline, n_wds_objects, n_hip_objects, n_hip_good_objects);
return(0);
}

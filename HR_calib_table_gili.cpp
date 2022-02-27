/************************************************************************
* "HR_calib_table_gili.c"
*
* To create a file with photometric and miscelaneous data
* in order to prepare various plots (HR diagram or statistics)
* Can create a Hertzprung-Russel diagram of the objects contained 
* in the calibrated table with the WDS catalog 
*
* For Gili's papers (without ADS name)
* 
* 2021: now ouput two other files: hr_without_dwarfs.txt hr_dwarfs_only.txt
*
* JLP 
* Version 03/05/2020
*************************************************************************/
#include <string.h>  // strstr() 
#include "jlp_catalog_utils.h"
#include "WDS_catalog_utils.h"
#include "jlp_string.h"  // jlp_compact_string

/*
#define DEBUG
#define DEBUG_1
*/

static int create_HR_of_calib_table_gili(FILE *fp_calib, char *WDS_catalog, 
                                        char *HIC_catalog,
                                        char *HIP_catalog, FILE *fp_out, 
                                        FILE *fp_out_dwarfs_tex,
                                        double paral_rel_error0, int ioption);
static int check_if_dwarf(double B_V_index, double magV_Hip_abs, 
                          double magV_A_WDS_abs, char *spectral_type, 
                          int *is_dwarf);

int main(int argc, char *argv[])
{
char calib_fname[80], out_fname[80], out_tex_fname[80];
char WDS_catalog[80], HIC_catalog[80], HIP_catalog[80];
int ioption;
FILE *fp_calib, *fp_out, *fp_out_dwarfs_tex;
time_t ttime = time(NULL);
double paral_rel_error0;
int status;

if(argc != 6 && argc != 7) {
  printf("Syntax: HR_calib_table_gili calib_table WDS_cat HIC_catalog HIP_main_catalog out_table [relative_error_of_parallax]\n");
  printf("Enter a negative relative error to output photometric data of all the objects\n");
  printf("Enter a positive relative error to output absolute magnitudes of objects measured by Hipparcos\n");
  printf("Example: HR_calib_table_gili tab_calib.tex wdsweb_summ.txt wds2hds2hip.txt Hipparcos_HIC.dat Hipparcos_main.dat Photom2a.dat -1\n");
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
/*
* ioption: 
*          1 photometric curve
*          2 HR curve
*/
if(paral_rel_error0 < 0) {
  ioption = 1;
  printf("OK: will output HR curve (not the photometric curve !)\n");
  printf("\n Now ouput two other files: hr_without_dwarfs.txt hr_dwarfs_only.txt \n");
} else {
  ioption = 2;
  printf("OK: will output photometric curve (not the HR curve !)\n");
}


printf("paral_rel_error=%f ioption=%d\n", paral_rel_error0, ioption);
printf("OK: calib=%s output=%s \n", calib_fname, out_fname); 
printf("OK: WDS=%s HIC=%s HIP=%s\n", WDS_catalog, HIC_catalog, HIP_catalog);

/* Open input calibrated table: */
if((fp_calib = fopen(calib_fname, "r")) == NULL) {
   fprintf(stderr, "HR_calib_table_gili/Fatal error opening calib. table %s\n",
           calib_fname);
    return(-1);
  }

/* Open output file: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "HR_calib_table_gili/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

if(ioption == 2) {
  strcpy(out_tex_fname, "HR_dwarfs.tex");
  if((fp_out_dwarfs_tex = fopen(out_tex_fname, "w")) == NULL) {
    fprintf(stderr, "HR_calib_table_gili/Fatal error opening output file: %s\n",
           out_tex_fname);
    return(-1);
   }
  fprintf(fp_out_dwarfs_tex, "\\begin{table} \n");
  fprintf(fp_out_dwarfs_tex, "\\centerline{ \n");
  fprintf(fp_out_dwarfs_tex, "\\begin{tabular}{lllcccc} \n");
  fprintf(fp_out_dwarfs_tex, "\\hline \n");
  fprintf(fp_out_dwarfs_tex, "wds_name & discov_name & spectral_type & magV_A_WDS & magV_Hip & B-V & Delta_mag & MagV \\\\ \n");
  fprintf(fp_out_dwarfs_tex, " & & & & & & \\\\ \n");
  fprintf(fp_out_dwarfs_tex, "\\hline \n");
 }

/* Header of the output photometric file: */
if(ioption == 1) {
  fprintf(fp_out, "%% Photometric data of %s\n%% Created on %s",
          calib_fname, ctime(&ttime));
  fprintf(fp_out, "%% magV_Hip magV_A_WDS B-V Delta_mag rho discov_name spectral_type\n");
  } else {
/* Header of the output HR file: */
  fprintf(fp_out, "%% HR diagram of %s (paral_rel_error0=%.3f) \n%% Created on %s",
        calib_fname, paral_rel_error0, ctime(&ttime));
  fprintf(fp_out, "%% magV_Hip magV_A_WDS B-V Delta_mag MagV MagV_A rho discov_name spectral_type\n");
  }

/* Scan the input calibrated table and generate output file 
* (compute M_V, B-V, and build a HR diagram if paral_rel_error0 > 0)
*/
status = create_HR_of_calib_table_gili(fp_calib, WDS_catalog,
                                      HIC_catalog, HIP_catalog, fp_out, 
                                      fp_out_dwarfs_tex,
                                      paral_rel_error0, ioption);

/* Close opened files:
*/
fclose(fp_calib);
fclose(fp_out);
if(ioption == 2) { 
  fprintf(fp_out_dwarfs_tex, "\\hline \n");
  fprintf(fp_out_dwarfs_tex, "\\end{tabular} \n");
  fprintf(fp_out_dwarfs_tex, "} \n");
  fprintf(fp_out_dwarfs_tex, "\\end{table} \n");
  fclose(fp_out_dwarfs_tex);
  }
return(status);
}
/************************************************************************
* Scan the calibrated table 
*
* INPUT:
* fp_calib: pointer to the input file containing the calibrated table 
* fp_out: pointer to the output Latex file with the full table 
* ioption: 
*          1 if photometric curve
*          2 if HR curve
*************************************************************************/
static int create_HR_of_calib_table_gili(FILE *fp_calib, char *WDS_catalog, 
                                        char *HIC_catalog, char *HIP_catalog, 
                                        FILE *fp_out, FILE *fp_out_dwarfs_tex,
                                        double paral_rel_error0, int ioption) 
{
double magV_Hip, magV_Hip_abs, magV_WDS, B_V_index, paral, err_paral;
double magV_A_WDS, magV_B_WDS, magV_A_WDS_abs, flux_A, flux_B;
double epoch0, rho0, err_rho0, theta0, err_theta0;
int iline, status, has_an_orbit, found_in_WDS, found_in_Hip_cat;
int n_wds_objects, n_hip_objects, n_hip_good_objects;
int is_dwarf, n_dwarfs, n_non_dwarfs, gili_format;
char in_line[256], wds_name[64], wds_name_calib[64]; 
char object_name[64], discov_name[64], comp_name[64], spectral_type[64];
char old_discov_name[64], out_fname[64];
FILE *fp_out_dwarfs, *fp_out_non_dwarfs;

/* Open output HR dwarfs and non dwarfs files: */
printf("ioption=%d\n", ioption);
if(ioption == 2) {
  n_dwarfs = 0;
  n_non_dwarfs = 0;
  strcpy(out_fname, "HR_dwarfs_only.txt");
  if((fp_out_dwarfs = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "HR_calib_table_gili/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }
  strcpy(out_fname, "HR_non_dwarfs.txt");
  if((fp_out_non_dwarfs = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "HR_calib_table_gili/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }
  }

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

      read_object_name_from_CALIB_line_gili(in_line, wds_name_calib, 
                                            discov_name, comp_name);
/* Remove the extra-blanks: */
      jlp_trim_string(wds_name_calib, 64);
      jlp_trim_string(discov_name, 64);
      jlp_trim_string(comp_name, 64);

// New: after 2018, now also read the measurements
// to be able to plot Delta mag versus rho measurements...
      gili_format = 1;
      status = read_measures_from_CALIB_line_gili(in_line, &epoch0, &rho0, 
                                             &err_rho0, &theta0, &err_theta0,
                                             gili_format);

// Check if it is a new object and that it was resolved:
      if(strcmp(old_discov_name, discov_name) && (status == 0)) {


       strcpy(old_discov_name, discov_name);

// Get magV_Hip, B_V_index, parall, err_paral from Hipparcos
// Get magV_A_WDS, magV_B_WDS, spectral_type from WDS 
       found_in_WDS = -1;
       found_in_Hip_cat = -1;
       get_data_from_WDS_and_HIP_catalogs(WDS_catalog, HIC_catalog, HIP_catalog,
                                          discov_name, comp_name,
                                          wds_name, &magV_Hip, 
                                          &B_V_index, &paral, &err_paral, 
                                          &magV_A_WDS, &magV_B_WDS, 
                                          spectral_type, 
                                          &found_in_WDS, &found_in_Hip_cat);
/* DEBUG
printf("wds_name=%s (wds_calib=%s) discov_name=%s magV_Hip=%f found_in_WDS=%d found_in_Hip_cat=%d\n",
      wds_name, wds_name_calib, discov_name, magV_Hip, found_in_WDS, found_in_Hip_cat);
*/

       if(found_in_WDS) n_wds_objects++;
       if(found_in_Hip_cat) n_hip_objects++;

/*
* ioption: 
*          1 photometric curve
*          2 HR curve
*/
/* Photometric data: */
         if(found_in_WDS && (ioption == 1)) {
              jlp_compact_string(discov_name, 20);
            if(found_in_Hip_cat) {
              fprintf(fp_out, "%.3f %.3f %.3f %.3f %.3f %s %s (magV_Hip and B_V from Hip)\n",
                      magV_Hip, magV_A_WDS, B_V_index, ABS(magV_A_WDS - magV_B_WDS),
                      rho0, discov_name, spectral_type);
// If not found in Hipparcos, use WDS data only:
// So has to generate magV, and  B_V_index 
            } else {
              B_V_index = 100.;
              flux_A = pow(10., -0.4 * magV_A_WDS);
              flux_B = pow(10., -0.4 * magV_B_WDS);
              magV_WDS = -2.5 * log10( flux_A + flux_B);
              fprintf(fp_out, "%.3f %.3f %.3f %.3f %.3f %s %s (magV and B_V from WDS)\n",
                      magV_WDS, magV_A_WDS, B_V_index, 
                      ABS(magV_A_WDS - magV_B_WDS),
                      rho0, discov_name, spectral_type);
            }
/* HR diagram: */
// magV_Hip, B_V_index, parall, err_paral from Hipparcos
// magV_A_WDS, magV_B_WDS, spectral_type from WDS 
         } else if(found_in_Hip_cat && (ioption == 2)) {
              if(paral > 0.) {
                if(err_paral/paral < paral_rel_error0) {
/* MV = mV - 2.5 log10 (d/d_10)^2
*  MV = mV - 5 log10 (d) + 5 log10(d_10)
*  MV = mV - 5 log10(d) + 5
* paral = 1/d
* MV = mV + 5 log10(paral) + 5
/* Parallax in mas: */
                magV_Hip_abs = magV_Hip + 5. * log10(paral * 0.001) + 5.;
                magV_A_WDS_abs = magV_A_WDS + 5. * log10(paral * 0.001) + 5.;
                jlp_compact_string(discov_name, 20);
                fprintf(fp_out, "%.3f %.3f %.3f %.3f %.3f %.3f %s %s\n",
                        magV_Hip, magV_A_WDS, B_V_index, 
                        ABS(magV_A_WDS - magV_B_WDS), 
                        magV_Hip_abs, magV_A_WDS_abs, discov_name, 
                        spectral_type);
                check_if_dwarf(B_V_index, magV_Hip_abs, magV_A_WDS_abs, 
                               spectral_type, &is_dwarf);
                if(is_dwarf == 1) {
                   n_dwarfs++;
                   fprintf(fp_out_dwarfs, "%.3f %.3f %.3f %.3f %.3f %.3f %s %s\n",
                        magV_Hip, magV_A_WDS, B_V_index, 
                        ABS(magV_A_WDS - magV_B_WDS), 
                        magV_Hip_abs, magV_A_WDS_abs, discov_name, 
                        spectral_type);
// wds_name & discov_name & spectral_type & magV_WDS_A & magV_Hip & B-V & Delta_mag & MagV 
                   fprintf(fp_out_dwarfs_tex, "%s & %s & %s & %.2f & %.2f & %.2f & %.2f \\\\ \n",
                        wds_name, discov_name, spectral_type, magV_A_WDS,
                        magV_Hip, B_V_index, ABS(magV_A_WDS - magV_B_WDS), 
                        magV_Hip_abs); 
                   } else {
                   fprintf(fp_out_non_dwarfs, "%.3f %.3f %.3f %.3f %.3f %.3f %s %s\n",
                        magV_Hip, magV_A_WDS, B_V_index, 
                        ABS(magV_A_WDS - magV_B_WDS), 
                        magV_Hip_abs, magV_A_WDS_abs, discov_name, 
                        spectral_type);
                   n_non_dwarfs++;
                   }
                n_hip_good_objects++;
                }
              } // If paral > 0
           } /* EOF found in Hip_cat and paral_rel_error > 0 */
      } // If discov_name != old_discov_name
     }/* EOF if !isdigit ... */
    }// EOF if inline ...
  } /* EOF if fgets */ 
 } /* EOF while ... */

printf("HR_curve_of_calib_table_gili: %d lines sucessfully read (%d wds_objects with %d in Hipparcos catalog and %d with parallax)\n", 
        iline, n_wds_objects, n_hip_objects, n_hip_good_objects);
fprintf(fp_out, "%% HR_curve_of_calib_table_gili: %d lines sucessfully read (%d wds_objects with %d in Hipparcos catalog and %d with parallax)\n", 
        iline, n_wds_objects, n_hip_objects, n_hip_good_objects);

if(ioption == 2) {
 fprintf(fp_out_dwarfs_tex, "%% HR_curve_of_calib_table_gili: %d wds_objects with %d in Hipparcos catalog and %d with parallax, n_dwarfs=%d\n", 
        n_wds_objects, n_hip_objects, n_hip_good_objects, n_dwarfs);
 fprintf(fp_out_dwarfs, "%% HR_curve_of_calib_table_gili: %d wds_objects with %d in Hipparcos catalog and %d with parallax, n_dwarfs=%d\n", 
        n_wds_objects, n_hip_objects, n_hip_good_objects, n_dwarfs);
 fprintf(fp_out_non_dwarfs, "%% HR_curve_of_calib_table_gili: %d wds_objects with %d in Hipparcos catalog and %d with parallax, n_non_dwarfs=%d\n", 
        n_wds_objects, n_hip_objects, n_hip_good_objects, n_non_dwarfs);

 fclose(fp_out_dwarfs);
 fclose(fp_out_non_dwarfs);
 }
return(0);
}
/****************************************************************************
*
****************************************************************************/
static int check_if_dwarf(double B_V_index, double magV_Hip_abs, 
                          double magV_A_WDS_abs, char *spectral_type, int *is_dwarf)
{
char M_word[4], K_word[4];
strcpy(M_word, "M");
strcpy(K_word, "K");

*is_dwarf = 0;
#ifdef DEBUG
printf("spectral type: >%s< ", spectral_type);
#endif

// strstr returns a pointer to the start of the word in "spectral_type" if
// the word is found
if(((strstr(spectral_type, "M") != NULL) || 
   (strstr(spectral_type, "K") != NULL)) && 
// I add the need for V:
   (strstr(spectral_type, "V") != NULL)) {
    *is_dwarf = 1;
  }
// Look then for super-giants, giants and sub-giants and discard them:
if((strstr(spectral_type, "III") != NULL) || 
   (strstr(spectral_type, "II") != NULL) ||
   (strstr(spectral_type, "IV") != NULL)) {
    *is_dwarf = 0;
  }
#ifdef DEBUG
printf("is_dwarf=%d\n", *is_dwarf);
#endif
return(0);
}

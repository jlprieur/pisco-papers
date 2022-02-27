/************************************************************************
* "residuals_gili_2.c"
*
* To compute the residuals of measurements of binary stars
* from known orbits.
*
* From o-c.for (version of 2008)
* written by Marco SCARDIA - Osservatorio Astronomico di Brera-Merate
*
* New features:
* - It is no longer necessary to indicate the number of objects
*   in the first line.
* - A line starting with % is interpreted as comments, and not
*   processed by the program.
* - Automatic precession correction of the input measurements
*   when the equinox of the orbit is very old (> 10 years)
* - It can retrieve PISCO measurements from the output LaTeX table 
*   created by latex_calib or any other LaTeX table in this format
* - It can read orbits from a subset of the OC6 catalog.
*
*
* OUTPUT:
*    *_curve.dat : ASCII file with data used for computing curves
*    *.txt : ASCII file with O-C
*    *.tex : LaTeX ASCII table with O-C
*    *_ref1.tex : Latex ASCII file with compacted references
*    *_ref2.tex : Latex ASCII file with full references
*
* JLP 
* Version 13/05/2021
*************************************************************************/
#include "jlp_catalog_utils.h"
#include "residuals_utils.h"
#include "OC6_catalog_utils.h"
#include "jlp_string.h"  // jlp_compact_string

/*
#define DEBUG
#define DEBUG_1
*/

static int residuals_gili_2_main(char *input_orbit_list, char *output_ext, 
                                 char *calib_fname, char *OC6_references_fname,
                                 int orbit_grade_max, int gili_format);
static int compute_residuals_from_calib_gili(FILE *fp_out_txt, 
                               FILE *fp_out_latex, FILE *fp_out_curve, 
                               FILE *fp_out_ref1, FILE *fp_out_ref2, 
                               char *calib_fname, char *input_orbit_list, 
                               char *OC6_references_fname, int orbit_grade_max,
                               int gili_format);
static int get_orbit_from_OC6catalog(char *input_orbit_list, 
                                     char *OC6_references_fname,
                                     char *object_name1, char *comp_name1, 
                                     int comp1_is_AB, char *object_name2,
                                     char *discov_name2, char *comp_name2,
                                     char *WDS_name2,
                                     double *Omega_node, double *omega_peri, 
                                     double *i_incl, double *e_eccent, 
                                     double *T_periastron, 
                                     double *orbit_equinox,
                                     double *mean_motion, double *a_smaxis, 
                                     double *Period, int *orbit_grade,
                                     char *author2, char *refer0, char *refer1);
static int process_measurement_gili(FILE *fp_out_txt, FILE *fp_out_latex,
              FILE *fp_out_curve, char *OC6_references_fname,
              char *object_name, char *discov_name, char *comp_name, 
              char *author, int orbit_grade,
              double Omega_node, double omega_peri, 
              double i_incl, double e_eccent,
              double T_periastron, double Period, double a_smaxis, 
              double mean_motion, double orbit_equinox, double epoch_o, 
              double rho_o, double theta_o, double err_rho_o, 
              double err_theta_o);

int main(int argc, char *argv[])
{
char input_orbit_list[80], output_ext[40], calib_fname[80];
char OC6_references_fname[128];
int orbit_grade_max, gili_format;

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
  printf("Syntax: residuals_gili_2 input_orbit_list output_ext calibrated_latex_table grade_max gili_format [reference_list] \n");
  printf(" Assume OC6 format (P, a, i, Omega, T, e, omep, [equinox] )\n");
  printf(" Test that is used: orbit_grade less or equal to orbit_grade_max\n");
  return(-1);
}
strcpy(input_orbit_list, argv[1]);
strcpy(output_ext, argv[2]);
/* Calibrated latex table */
strcpy(calib_fname, argv[3]);
/* orbit_grade_max */
if(sscanf(argv[4], "%d", &orbit_grade_max) != 1) {
  printf("Fatal error reading orbit_grade_max ! (argv[4]=%s)\n", argv[4]);
  return(-1);
  }
sscanf(argv[5], "%d", &gili_format);

/* File with full references (not necessary) */
 if(argc == 7) 
   strcpy(OC6_references_fname, argv[6]);
 else 
   OC6_references_fname[0] = '\0';

#ifdef DEBUG
printf("OK: input_orbit_list=%s output_ext=%s calib_fname=%s\n", 
       input_orbit_list, output_ext, calib_fname);
printf("OK: orbit_grade_max=%d\n", orbit_grade_max);
printf("OK: OC6_references_fname=>%s< gili_format=%d\n", 
        OC6_references_fname, gili_format);
#endif

/* Call residuals1_main that does the main job: */
residuals_gili_2_main(input_orbit_list, output_ext, calib_fname, 
                      OC6_references_fname, orbit_grade_max, gili_format);

return(0);
}
/************************************************************************
* residuals_gili_2_main
* main routine of "residuals_2.c"
* Open the input/output files needed to compute the residuals 
* and call all other routines.
*
* INPUT:
* input_orbit_list: name of the file containing the orbital elements
* output_ext: extension of the output files
* calib_fname: name of the file containing the Latex calibrated table
*              (final version, ready for publication)
* OC6_references_fname: name of the file containing the OC6 biblio. references
*
*************************************************************************/
static int residuals_gili_2_main(char* input_orbit_list, char *output_ext, 
                                 char *calib_fname, char *OC6_references_fname,
                                 int orbit_grade_max, int gili_format)
{
char out_filename[100];
FILE *fp_out_txt, *fp_out_latex, *fp_out_curve; 
FILE *fp_out_ref1, *fp_out_ref2;
time_t t = time(NULL);

fp_out_ref1 = NULL;
fp_out_ref2 = NULL;
fp_out_curve = NULL;

/* Open output curve for O-C plot: */
  sprintf(out_filename, "%s_curve.dat", output_ext);
  if((fp_out_curve = fopen(out_filename, "w")) == NULL) {
    fprintf(stderr, "residuals2_main/Fatal error opening output calibrated latex file: %s\n",
           out_filename);
    return(-1);
   }

if(*OC6_references_fname) {
/* Open output curve for compacted references: */
  sprintf(out_filename, "%s_ref1.tex", output_ext);
  if((fp_out_ref1 = fopen(out_filename, "w")) == NULL) {
    fprintf(stderr, "residuals2_main/Fatal error opening output latex file: %s\n",
           out_filename);
    return(-1);
   }
fprintf(fp_out_ref1, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out_ref1, "%% Residuals from: %s and %s, (orbit_grade_max=%d), computed on %s", 
        calib_fname, input_orbit_list, orbit_grade_max, ctime(&t));
fprintf(fp_out_ref1, "%% Created by residuals_2.c -- JLP version of 23/05/2021 --\n%% \n");
fprintf(fp_out_ref1, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");


/* Open output curve for full references: */
  sprintf(out_filename, "%s_ref2.tex", output_ext);
  if((fp_out_ref2 = fopen(out_filename, "w")) == NULL) {
    fprintf(stderr, "residuals2_main/Fatal error opening output latex file: %s\n",
           out_filename);
    return(-1);
   }
fprintf(fp_out_ref2, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out_ref2, "%% Residuals from: %s and %s, (orbit_grade_max=%d), computed on %s", 
        calib_fname, input_orbit_list, orbit_grade_max, ctime(&t));
fprintf(fp_out_ref2, "%% Created by residuals_2.c -- JLP version of 23/05/2021 --\n%% \n");
fprintf(fp_out_ref2, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
}  /* EOF case *OC6_references_fname != 0 */

/* Open output ASCII text file in plain format: */
sprintf(out_filename, "%s.txt", output_ext);
if((fp_out_txt = fopen(out_filename, "w")) == NULL) {
   fprintf(stderr, "residuals1_main/Fatal error opening output text file: %s\n",
           out_filename);
    return(-1);
  }
fprintf(fp_out_txt, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out_txt, "%% Residuals from: %s and %s, (orbit_grade_max=%d), computed on %s", 
        calib_fname, input_orbit_list, orbit_grade_max, ctime(&t));
fprintf(fp_out_txt, "%% Created by residuals_2.c -- JLP version of 23/05/2021 --\n%% \n");
fprintf(fp_out_txt, "%% Name  Epoch  rho_O  rho_C  Drho_O-C  theta_O  theta_C  Dtheta_O-C  Author\n");
fprintf(fp_out_txt, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

/* Open output Latex file: */
sprintf(out_filename, "%s.tex", output_ext);
if((fp_out_latex = fopen(out_filename, "w")) == NULL) {
   fprintf(stderr, "Compute_residuals/Fatal error opening output Latex file: %s\n",
           out_filename);
    return(-1);
  }

/* Header of the output Latex table: */
fprintf(fp_out_latex, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out_latex, "%% Residuals from: %s and %s, (orbit_grade_max=%d), computed on %s", 
        calib_fname, input_orbit_list, orbit_grade_max, ctime(&t));
fprintf(fp_out_latex, "%% Created by residuals_gili_2.c -- JLP version of 23/05/2021 --\n");
fprintf(fp_out_latex, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out_latex, "\\begin{table} \n\\begin{center} \n\
\\caption{Residuals of the measurements of Table 1 with published orbits.} \n\
\\begin{tabular}{llcccrrl} \n \\hline \n\
 Name  &   Orbit   & Epoch  & $\\rho$(O) $\\theta$(O) & \n\
& $\\Delta \\rho$(O-C) & $\\Delta \\theta$(O-C) & Grade\\\\ \n\
& & & & (\\arcsec) & ($^\\circ$) \\\\ \n\
\\hline \n & & & & & & & \\\\ \n");

/* Header of the output curve: */
fprintf(fp_out_curve, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out_curve, "%% Residuals from: %s and %s, (orbit_grade_max=%d), computed on %s", 
        calib_fname, input_orbit_list, orbit_grade_max, ctime(&t));
fprintf(fp_out_curve, "%% Created by residuals_2.c -- JLP version of 13/04/2010 --\n");
fprintf(fp_out_curve, "%% Drho(O-C) Dtheta(O-C) err_rho_O err_theta_O\n");
fprintf(fp_out_curve, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

/* Scan the input orbit file and compute the residuals 
*/
compute_residuals_from_calib_gili(fp_out_txt, fp_out_latex, 
                                  fp_out_curve, fp_out_ref1, fp_out_ref2, 
                                  calib_fname, input_orbit_list, 
                                  OC6_references_fname, orbit_grade_max,
                                  gili_format);

/* Epilog for Latex file: */
fprintf(fp_out_latex, " & & & & & & & \\\\ \n \\hline \n \\end{tabular} \n\
\\end{center} \n \\end{table*} \n");

/* Close opened files:
*/
fclose(fp_out_txt);
fclose(fp_out_latex);
if(fp_out_ref1) fclose(fp_out_ref1);
if(fp_out_ref2) fclose(fp_out_ref2);
if(fp_out_curve) fclose(fp_out_curve);
return(0);
}
/************************************************************************
* compute_residuals_from_calib_gili
* Scan the input list from the orbit file,
* look for the measures in the tab_calib.tex 
* and compute the residuals 
*
* INPUT:
* fp_out_txt: pointer to the output file with the residuals and various data 
*             in plain ASCII format 
* fp_out_latex: pointer to the output Latex file with the residuals
* fp_out_curve: pointer to the file containing the O-C curve
* fp_out_ref1: pointer to the file with compacted references
* fp_out_ref2: pointer to the file with full references
*
*************************************************************************/
static int compute_residuals_from_calib_gili(FILE *fp_out_txt, 
                              FILE *fp_out_latex, FILE *fp_out_curve, 
                              FILE *fp_out_ref1, FILE *fp_out_ref2, 
                              char *calib_fname, char *input_orbit_list, 
                              char *OC6_references_fname, int orbit_grade_max,
                              int gili_format)
{
#define NMAX 1024
char author[NMAX*60], refer0[NMAX*130], refer1[NMAX*130], object_name[NMAX*60]; 
char author2[60], refer20[130], refer21[130], object_name2[60];
char discov_name2[64], comp_name2[64], WDS_name2[64];
double Omega_node, omega_peri, i_incl, e_eccent, T_periastron, orbit_equinox;
double mean_motion, a_smaxis, Period; 
double epoch_o, rho_o, theta_o, err_rho_o, err_theta_o;
char in_line[300], object_name1[40], comp_name1[40];
int iline, comp1_is_AB, kk, orbit_grade = -1, status, status1;
FILE *fp_in_latex;

/* Open input file containing the input Latex table: */
if((fp_in_latex = fopen(calib_fname, "r")) == NULL) {
  fprintf(stderr, "compute_residuals_from_calib_gili/Error opening %s\n",
          calib_fname);
  return(-1);
 }

/* Scan all the LaTeX file looking for the object name */
/* Example:
11000$-$0328  &  STF1500  & 2010.387 & 2 & 1.370 & 0.015 & 299.9  & 0.3 & 1.20 &   \\
*/
iline = 0;
kk = 0;
while(!feof(fp_in_latex)){
  if(fgets(in_line,300,fp_in_latex)) {
  iline++;
/* Process only the meaningful lines (skipping the header...)*/
  if(isdigit(in_line[0])) {
    read_full_name_from_CALIB_line_gili(in_line, object_name1, comp_name1,
                                        &comp1_is_AB);
printf("From read_full_name_from_CALIB_line_gili: object=%s comp=%s\n", object_name1, comp_name1);
// JLP2020: reduce the name if comp_is_AB is true:
    if(comp1_is_AB == 1) remove_AB_from_object_name(object_name1);

    status = read_measures_from_CALIB_line_gili(in_line, &epoch_o,
                                       &rho_o, &err_rho_o,
                                       &theta_o, &err_theta_o, gili_format);
printf("ZZZZ: rho_o=%f \n", rho_o);
    if(status > 0) {
       fprintf(stderr, "Error retrieving data for object_name=%s %s (Unresolved ?)\n",
            object_name1, comp_name1);
       fprintf(stderr, ">%s<\n", in_line);
     } else if(status < 0) {
       fprintf(stderr, "read_measures_from_CALIB_table_gili/Fatal error reading calibrated file \n");
       fprintf(stderr, ">%s<\n", in_line);
       exit(-1);
     } else {
       jlp_compact_string(object_name1, 40);
       jlp_compact_string(comp_name1, 40);
       status1 = get_orbit_from_OC6catalog(input_orbit_list, 
                                          OC6_references_fname,
                                          object_name1, comp_name1, comp1_is_AB,
                                          object_name2, discov_name2,
                                          comp_name2, WDS_name2, &Omega_node, 
                                          &omega_peri, &i_incl, &e_eccent, 
                                          &T_periastron, &orbit_equinox,
                                          &mean_motion, &a_smaxis, &Period,
                                          &orbit_grade, author2, refer20, 
                                          refer21); 
      if(status1 != 0) {
        fprintf(stderr, 
          "compute_residuals_gili/Error in get_orbit_from_OC6catalog object=%s comp=%s status1%d\n",
               object_name1, comp_name1, status1); 
// Case when an orbit has been found with a good grade for that object: 
       } else if((orbit_grade > 0) 
                && (orbit_grade <= orbit_grade_max)) {
         strcpy(&object_name[kk * 60], object_name2); 
         strcpy(&author[kk * 60], author2); 
         strcpy(&refer0[kk * 130], refer20); 
         strcpy(&refer1[kk * 130], refer21); 
/**********************
/* Process this measurement
* correct the measurement from precession, compute the ephemerides,
* and derive the O-C residuals
*/
         status = process_measurement_gili(fp_out_txt, fp_out_latex, 
                                     fp_out_curve, OC6_references_fname, 
                                     object_name2, discov_name2, comp_name2, 
                                     author2, orbit_grade,  
                                     Omega_node, omega_peri, i_incl, e_eccent, 
                                     T_periastron, Period, a_smaxis, 
                                     mean_motion, orbit_equinox, epoch_o, 
                                     rho_o, theta_o, err_rho_o, 
                                     err_theta_o);
         if(status) {
           fprintf(stderr, 
               "compute_residuals_gili/Error processing measurement in line #%d\n",
                iline); 
          }

/* Update kk the orbit counter: */
         kk++;
         }
      } // EOF status == 0 ...
    } // EOF isdigit...
  } // EOF if gets ...
} // EOF while...

fclose(fp_in_latex);
/* Sort the references by alphabetic order and save them to output
* Latex files: */
if(fp_out_ref1 != NULL && fp_out_ref2 != NULL) 
            sort_references(fp_out_ref1, fp_out_ref2, object_name, author, 
                            refer0, refer1, kk);

printf("Compute_residuals: %d lines sucessfully read (n_orbits=%d)\n", 
        iline, kk);

return(0);
}
/************************************************************************
* Process one measurement
* correct the measurement from precession, compute the ephemerides,
* and derive the O-C residuals
*
* INPUT:
* fp_out_txt: pointer to the output file with the residuals and various data 
*             in plain ASCII format 
* fp_out_latex: pointer to the output Latex file with the residuals
* fp_out_curve: pointer to the output file containing the O-C curve 
*
* OUTPUT:
* O-C residuals in "fp_out_txt" and "fp_out_latex" files
************************************************************************/
static int process_measurement_gili(FILE *fp_out_txt, FILE *fp_out_latex,
              FILE *fp_out_curve, 
              char *OC6_references_fname, char *object_name, char *discov_name,
              char *comp_name, char *author, int orbit_grade, 
              double Omega_node, 
              double omega_peri, double i_incl, double e_eccent,
              double T_periastron, double Period, double a_smaxis, 
              double mean_motion, double orbit_equinox, double epoch_o, 
              double rho_o, double theta_o, double err_rho_o, 
              double err_theta_o)
{
char my_name[60], quadrant_discrep[20];
double c_tolerance, rho_c, theta_c, Drho, Dtheta;
int status;
register int i;

/*  c_tolerance = smallest increment allowed in the iterative process
*                used for solving Kepler's equation
*/
c_tolerance = ABS(1.5E-5 * cos(i_incl) 
                   / sqrt((1.0 + e_eccent)/(1.0 - e_eccent)));

 printf("process_measurement/Using orbit for object=%s\n ", object_name);

/* Conversion to radians: */
 theta_o *= DEGTORAD;

/* Compute the ephemerids corresponding to the observation epoch: */
compute_ephemerid(Omega_node, omega_peri, i_incl, e_eccent, T_periastron, 
                  Period, a_smaxis, mean_motion, epoch_o, c_tolerance, 
                  &theta_c, &rho_c);

/* Conversion to degrees: */
 theta_o /= DEGTORAD;
 Dtheta = theta_o - theta_c;
 if(Dtheta < -300.0) Dtheta += 360.0;
 if(Dtheta > 300.0) Dtheta -= 360.0;

/* Special handling of Dtheta when close to 180 or -180 degrees: */
 strcpy(quadrant_discrep, "");
 if(ABS(Dtheta - 180.) < 60.) {
   Dtheta -= 180.; 
   strcpy(quadrant_discrep, "$^Q$");
   } if(ABS(Dtheta + 180.) < 60.) {
   Dtheta += 180.; 
   strcpy(quadrant_discrep, "$^Q$");
   }

 Drho = rho_o - rho_c;
/* Trick to have a constant width */
 sprintf(my_name, "%s %s", object_name, comp_name);
/* Left justified text is obtained with a minus sign in the format:*/
 fprintf(fp_out_txt, "%-18.18s %9.3f %9.3f %9.3f %8.2f %8.2f %8.2f %7.1f %s %d\n",
         my_name, epoch_o, rho_o, rho_c, Drho, theta_o, theta_c, 
         Dtheta, author, orbit_grade);
 fprintf(fp_out_latex, "%s %s & %s & %9.3f & %9.3f & %6.1f & %8.2f & %8.2f%s & %d \\\\\n",
         object_name, comp_name, author, epoch_o, rho_o, theta_o, Drho, Dtheta,
         quadrant_discrep, orbit_grade);
 fprintf(fp_out_curve, "%8.3f %7.2f %8.3f %7.2f %9.3f %-18.18s %s %d\n",
         Drho, Dtheta, err_rho_o, err_theta_o, epoch_o, my_name, author,
         orbit_grade);

return(0);
}
/********************************************************************
*
********************************************************************/
static int get_orbit_from_OC6catalog(char *input_orbit_list, 
                                     char *OC6_references_fname,
                                     char *object_name1, char *comp_name1, 
                                     int comp1_is_AB, char *object_name2,
                                     char *discov_name2, char *comp_name2,
                                     char *WDS_name2,
                                     double *Omega_node, double *omega_peri, 
                                     double *i_incl, double *e_eccent, 
                                     double *T_periastron, 
                                     double *orbit_equinox,
                                     double *mean_motion, double *a_smaxis, 
                                     double *Period, int *orbit_grade,
                                     char *author2, char *refer0, char *refer1)
{
double Omega_node2, omega_peri2, i_incl2, e_eccent2, T_periastron2; 
double orbit_equinox2, mean_motion2, a_smaxis2, Period2; 
char comp1_really_compacted[32], comp2_really_compacted[32]; 
FILE *fp_in;
/* Maximum line seems to be 265 for OC6 catalog... */
char in_line1[300];
int iline, is_master_file, line_length, status, comp2_is_AB;

jlp_really_compact_companion(comp_name1, comp1_really_compacted, 40);

/* Open input file containing the orbital parameters: */
if((fp_in = fopen(input_orbit_list, "r")) == NULL) {
   fprintf(stderr, "get_orbit_from_OC6catalog/Fatal error opening file: %s\n",
           input_orbit_list);
    return(-1);
  }

// fp_in: pointer to the input file containing the orbital elements
iline = 0;
while(!feof(fp_in)){
/* Read input line and retrieve the orbital parameters 
*/
  if(fgets(in_line1,300,fp_in)) {
    line_length = (int)strlen(in_line1);
/* strlen(line1) = 279 if master file
*  strlen(line1) = 265 if OC6 file
*/
    is_master_file = (line_length > 270) ? 1 : 0;
    iline++;
/* Commented lines can start with % or # : */
    if(in_line1[0] != '%' && in_line1[0] != '#') {
/* Check if very short line (such as the one with n_orbits, used by Marco's
* program
*/
    if(line_length < 10) {
      printf("WARNING line #%d is very short (length=%d): >%s<\n", 
              iline, line_length, in_line1);
    } else {

/* OC6 format: */
     status = get_orbit_from_OC6_list_gili(in_line1, iline, 
                                 is_master_file, WDS_name2, discov_name2, 
                                 comp_name2, object_name2, author2, 
                                 &Omega_node2, &omega_peri2, &i_incl2, 
                                 &e_eccent2, &T_periastron2, &Period2, 
                                 &a_smaxis2, &mean_motion2, &orbit_equinox2,
                                 orbit_grade);
     if(status) {
      fprintf(stderr, "get_orbit_from_OC6catalog/WARNING: error reading orbital parameters in line #%d (status=%d)\n", iline, status); 
      } else {
      jlp_compact_string(object_name2, 60);
/* Check if AB companion (default) or something else */
      comp2_is_AB = 0;
      if((comp_name2[0] == '\0') || !strcmp(comp_name2,"AB")
         || !strncmp(comp_name1,"Aa-B",4)) comp1_is_AB = 1;
// JLP2020: reduce the name if comp_is_AB is true:
      if(comp2_is_AB == 1) remove_AB_from_object_name(object_name2);
      jlp_really_compact_companion(comp_name2, comp2_really_compacted, 40);
      if(!strcmp(object_name1, object_name2)
         && !strcmp(comp1_really_compacted, comp2_really_compacted)) {
        *Omega_node = Omega_node2;
        *omega_peri = omega_peri2; 
        *i_incl = i_incl2;
        *e_eccent = e_eccent2; 
        *T_periastron = T_periastron2;
        *orbit_equinox = orbit_equinox2; 
        *mean_motion = mean_motion2;
        *a_smaxis = a_smaxis2; 
        *Period = Period2;
        fclose(fp_in);
        if(*OC6_references_fname) { 
          get_OC6_full_reference(object_name2, author2, OC6_references_fname, 
                                 refer0, refer1); 
/* DEBUG:
printf("object=%s author=%s refer0=%s refer1=%s\n", object_name2, 
         author2, refer0, refer1);
*/
        } else {
          strcpy(author2, "");
          strcpy(refer0, "");
          strcpy(refer1, "");
        }
// Return 0 if object was found 
        return(0);
       } // EOF !strncmp(object_name1, object_name2);
      } // EOF if status == 0
     } // EOF if long line 
   } // EOF not commented line
  } // EOF if fgets
} // EOF while

fclose(fp_in);
// Return -2 if object was not found 
return(-2);
}

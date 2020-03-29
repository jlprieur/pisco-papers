/************************************************************************
* "residuals_gili_1.c"
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
* Version 13/04/2009
*************************************************************************/
#include "jlp_catalog_utils.h"
#include "residuals_utils.h"
#include "OC6_catalog_utils.h"

/*
#define DEBUG
#define DEBUG_1
*/

static int residuals_gili_1_main(char* input_filename, char *output_ext, 
                                 char *calib_fname, 
                                 char *OC6_references_fname, int iformat);
static int compute_residuals_gili(FILE *fp_in, FILE *fp_out_txt, 
                                  FILE *fp_out_latex, FILE *fp_out_curve, 
                                  FILE *fp_out_ref1, FILE *fp_out_ref2, 
                                  char *calib_fname, 
                                  char *OC6_references_fname, int iformat);
static int get_orbit_from_Marco_list(char *in_line1, char *in_line2, 
              int iline, 
              char *object_name, char *WDS_name, char *ADS_name, 
              char *discov_name, char *comp_name, char *author,
              double *Omega_node, double *omega_peri, double *i_incl, 
              double *e_eccent, double *T_periastron, double *Period, 
              double *a_smaxis, double *mean_motion, double *orbit_equinox);
static int get_measures_from_file(FILE *fp_in, int *iline, 
                                  double *epoch_o, double *rho_o, 
                                  double *theta_o, double *err_rho_o, 
                                  double *err_theta_o, int *nmeas);
static int process_measurements_gili(FILE *fp_out_txt, FILE *fp_out_latex,
              FILE *fp_out_curve, char *OC6_references_fname,
              char *object_name, char *discov_name, char *comp_name, 
              char *author,
              double Omega_node, double omega_peri, double i_incl, double e_eccent,
              double T_periastron, double Period, double a_smaxis, 
              double mean_motion, double orbit_equinox, double *epoch_o, 
              double *rho_o, double *theta_o, double *err_rho_o, 
              double *err_theta_o, int nmeas);
static int read_object_name1(char *in_line, char *object_name, 
                             char *WDS_name, char *ADS_name, char *discov_name,
                             char *comp_name, char *author, int iline);

int main(int argc, char *argv[])
{
char input_filename[80], output_ext[40], calib_fname[80];
char OC6_references_fname[128];
int iformat;

if(argc == 7) {
  if(*argv[6]) argc = 7;
  else if(*argv[5]) argc = 6;
  else if(*argv[4]) argc = 5;
  else if(*argv[3]) argc = 4;
  else if(*argv[2]) argc = 3;
  else if(*argv[1]) argc = 2;
  else argc = 1;
}
if(argc != 4 && argc != 5 && argc != 6) {
  printf("Syntax: residuals_1 input_list input_format output_ext [calibrated_latex_table] [reference_list] \n");
  printf("Format: -1 if Marco's format (Omega=node, omep=longitude of periastron, i, e, T, P, a, [equinox]) with measures\n");
  printf("        1 if Marco's format without measures\n");
  printf("        2 if OC6 format (P, a, i, Omega, T, e, omep, [equinox] )\n");
  return(-1);
}
strcpy(input_filename, argv[1]);
sscanf(argv[2], "%d", &iformat);
strcpy(output_ext, argv[3]);
/* Calibrated latex table (not necessary if iformat < 0) */
if(argc >= 6) {
strcpy(calib_fname, argv[4]);
} else {
calib_fname[0] = '\0';
}
/* Error messages: */
if((iformat > 0) && (calib_fname[0] == '\0')) {
  fprintf(stderr, "Fatal error\n");
  fprintf(stderr, "iformat=%d: Calibrated table is needed: please provide filename\n",
          iformat);
  return(-1);
  }
/* File with full references (not necessary if iformat != 2) */
if(iformat == 2) {
 if(argc == 7) {
   strcpy(OC6_references_fname, argv[6]);
 } else {
 OC6_references_fname[0] = '\0';
 }
}

#ifdef DEBUG
printf("OK: input=%s iformat=%d output_ext=%s calib_fname=%s\n", 
       input_filename, iformat, output_ext, calib_fname);
printf("OK: OC6_references_fname=>%s<\n", OC6_references_fname);
#endif

/* Call residuals1_main that does the main job: */
residuals_gili_1_main(input_filename, output_ext, calib_fname, 
                      OC6_references_fname, iformat);

return(0);
}
/************************************************************************
* residuals_gili_1_main
* main routine of "residuals_1.c"
* Open the input/output files needed to compute the residuals 
* and call all other routines.
*
* INPUT:
* input_filename: name of the file containing the measurements and the
*                 orbital elements
* output_ext: extension of the output files
* calib_fname: name of the file containing the Latex calibrated table
*              (final version, ready for publication)
* OC6_references_fname: name of the file containing the OC6 biblio. references
* iformat: format of the input file (-1=Marco with measures, 
*          1=Marco without measures, 2=OC6 without measures)
*
*************************************************************************/
static int residuals_gili_1_main(char* input_filename, char *output_ext, 
                                 char *calib_fname, 
                                 char *OC6_references_fname, int iformat)
{
char out_filename[100];
FILE *fp_in, *fp_out_txt, *fp_out_latex, *fp_out_curve; 
FILE *fp_out_ref1, *fp_out_ref2;
time_t t = time(NULL);

fp_out_ref1 = NULL;
fp_out_ref2 = NULL;
fp_out_curve = NULL;

/* Open input file containing the measurements and the orbital parameters: */
if((fp_in = fopen(input_filename, "r")) == NULL) {
   fprintf(stderr, "residuals1_main/Fatal error opening input file: %s\n",
           input_filename);
    return(-1);
  }

/* Open output curve for O-C plot: */
  sprintf(out_filename, "%s_curve.dat", output_ext);
  if((fp_out_curve = fopen(out_filename, "w")) == NULL) {
    fprintf(stderr, "residuals1_main/Fatal error opening output calibrated latex file: %s\n",
           out_filename);
    return(-1);
   }

if(*OC6_references_fname) {
/* Open output curve for compacted references: */
  sprintf(out_filename, "%s_ref1.tex", output_ext);
  if((fp_out_ref1 = fopen(out_filename, "w")) == NULL) {
    fprintf(stderr, "residuals1_main/Fatal error opening output latex file: %s\n",
           out_filename);
    return(-1);
   }
fprintf(fp_out_ref1, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out_ref1, "%% Residuals from: %s, computed on %s", 
        input_filename, ctime(&t));
fprintf(fp_out_ref1, "%% Created by residuals_1.c -- JLP version of 13/04/2010 --\n%% \n");
fprintf(fp_out_ref1, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");


/* Open output curve for full references: */
  sprintf(out_filename, "%s_ref2.tex", output_ext);
  if((fp_out_ref2 = fopen(out_filename, "w")) == NULL) {
    fprintf(stderr, "residuals1_main/Fatal error opening output latex file: %s\n",
           out_filename);
    return(-1);
   }
fprintf(fp_out_ref2, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out_ref2, "%% Residuals from: %s, computed on %s", 
        input_filename, ctime(&t));
fprintf(fp_out_ref2, "%% Created by residuals_1.c -- JLP version of 13/04/2010 --\n%% \n");
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
fprintf(fp_out_txt, "%% Residuals from: %s, computed on %s", 
        input_filename, ctime(&t));
fprintf(fp_out_txt, "%% Created by residuals_1.c -- JLP version of 13/04/2010 --\n%% \n");
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
fprintf(fp_out_latex, "%% Residuals from: %s, computed on %s", 
        input_filename, ctime(&t));
fprintf(fp_out_latex, "%% Created by residuals_gili_1.c -- JLP version of 13/02/2020 --\n");
fprintf(fp_out_latex, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out_latex, "\\begin{table} \n\\begin{center} \n\
\\caption{Residuals of the measurements of Table 1 with published orbits.} \n\
\\begin{tabular}{llccrr} \n \\hline \n\
 Name  &   Orbit   & Epoch  & $\\rho$(O) \n\
& $\\Delta \\rho$(O-C) & $\\Delta \\theta$(O-C) \\\\ \n\
& & & & (\\arcsec) & ($^\\circ$) \\\\ \n\
\\hline \n & & & & & \\\\ \n");

/* Header of the output curve: */
fprintf(fp_out_curve, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out_curve, "%% Residuals from: %s, computed on %s", 
        input_filename, ctime(&t));
fprintf(fp_out_curve, "%% Created by residuals_1.c -- JLP version of 13/04/2010 --\n");
fprintf(fp_out_curve, "%% Drho(O-C) Dtheta(O-C) err_rho_O err_theta_O\n");
fprintf(fp_out_curve, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

/* Scan the input orbit file and compute the residuals 
*/
compute_residuals_gili(fp_in, fp_out_txt, fp_out_latex, fp_out_curve, 
                       fp_out_ref1, fp_out_ref2, calib_fname, 
                       OC6_references_fname, iformat);

/* Epilog for Latex file: */
fprintf(fp_out_latex, " & & & & & \\\\ \n \\hline \n \\end{tabular} \n\
\\end{center} \n \\end{table*} \n");

/* Close opened files:
*/
fclose(fp_in);
fclose(fp_out_txt);
fclose(fp_out_latex);
if(fp_out_ref1) fclose(fp_out_ref1);
if(fp_out_ref2) fclose(fp_out_ref2);
if(fp_out_curve) fclose(fp_out_curve);
return(0);
}

/***************************************************************************
* get_measures_from_file
* Retrieve the measurements from the input file (when iformat < 0)
*
* Examples:
*
* 2007.969 0.656 155.3 1
* 2007.969 0.653 154.4 0
*
* or:
*
* 2007.9693 0.825 18.0 0
*
* INPUT:
* fp_in: pointer to the input file containing the measurements and the
*        orbital elements
* 
* OUTPUT:
* object designation, orbital parameters and arrays with the measurements
*
***************************************************************************/
static int get_measures_from_file(FILE *fp_in, int *iline, 
                                  double *epoch_o, double *rho_o, 
                                  double *theta_o, double *err_rho_o, 
                                  double *err_theta_o, int *nmeas)
{
char in_line[300];
double epoch, rho, theta;
int kk, nval;

/* Read the measurements from the input file if iformat < 0: */
   *nmeas = 0;
   do {
   in_line[0] = '%';
   while(!feof(fp_in) && in_line[0] == '%') {
    (*iline)++;
    if(!fgets(in_line,300,fp_in)) {
    fprintf(stderr, "Process_new_orbit/Error reading input file at line #%d \n %s \n", 
            *iline, in_line);
    return(-2);
    }
   } 

   nval = sscanf(in_line, "%lf %lf %lf %d", &epoch, &rho, &theta, &kk);
   if(nval != 4) {
    fprintf(stderr, "Process_new_orbit/Error reading input file at line #%d \n %s \n", 
            *iline, in_line);
    return(-2);
    }
    epoch_o[*nmeas] = epoch;
    rho_o[*nmeas] = rho;
    theta_o[*nmeas] = theta;
    err_rho_o[*nmeas] = 0.;
    err_theta_o[*nmeas] = 0.;
    (*nmeas)++;

#ifdef DEBUG_1
   printf("iline=%d/measurements: %f %f %f %d\n", 
           *iline, epoch, rho, theta, kk);
#endif
  } while(kk == 1); /* EOF while kk == 1 */

return(0);
}
/***************************************************************************
* get_orbit_from_Marco_list
* Read input line and retrieve the orbital parameters and the measurements
* in Marco's format (i.e., iformat = 1 or -1)
*
* Input lines contain
*  - name of object 
*  - the orbital parameters in Marco's format (if iformat = 1 or -1) 
*  - the measurements (if iformat = -1)
*    otherwise, (if format = 1) retrieve the measurements 
*      from the calibrated table.
*
* Examples:
*
* With iformat = -1:
*
* ADS 293 - Ole2001
* 149.0 182.7 128.3 0.9424 1928.76 335.385 0.3935 
* 2007.969 0.656 155.3 1
* 2007.969 0.653 154.4 0
*
* or:
*
* 00093+7943 STF   2 = ADS 102 - Hei1997 - Heintz (1997)
* 171.2 333.7 110.1 0.715 1887.5 540.0 0.995
* 2007.9693 0.825 18.0 0
*
* With iformat = 1:
*
* ADS 293 - Ole2001
* 149.0 182.7 128.3 0.9424 1928.76 335.385 0.3935 
*
* INPUT:
* in_line1: line from the the input file, containing the object name
* in_line2: line from the the input file, containing the orbital elements
* iline: number of the line corresponding to in_line2 in input file
* iformat: format of the input file (-1=Marco with measures, 
*          1=Marco without measures, 2=OC6 without measures)
* 
* OUTPUT:
* object designation, orbital parameters and arrays with the measurements
*
***************************************************************************/
static int get_orbit_from_Marco_list(char *in_line1, char *in_line2, 
              int iline, char *object_name, char *WDS_name, char *ADS_name, 
              char *discov_name, char *comp_name, char *author,
              double *Omega_node, double *omega_peri, double *i_incl, 
              double *e_eccent, double *T_periastron, double *Period, 
              double *a_smaxis, double *mean_motion, double *orbit_equinox)
{
int nval, status;

/* Decode the first line containing the object designation
* and the author of the orbit
*/
status = read_object_name1(in_line1, object_name, WDS_name, ADS_name, 
                          discov_name, comp_name, author, iline - 1);

/* Read orbital elements from input file: 
* (in Marco's format) */
 nval = sscanf(in_line2, "%lf %lf %lf %lf %lf %lf %lf %lf", 
                Omega_node, omega_peri, i_incl, e_eccent, T_periastron, 
                Period, a_smaxis, orbit_equinox);
 if(nval != 7 && nval != 8) {
  fprintf(stderr, "Process_new_orbit/Error reading input file at line #%d (nval=%d) \n %s \n", 
          iline, nval, in_line2);
  return(-2);
  }
 if(nval == 7) *orbit_equinox = 2000.0;

#ifdef DEBUG_1
 printf("nval=%d Omega_node=%f omega_peri=%f incl=%f e=%f T=%f P=%f a=%f Equinox=%f\n", 
        nval, *Omega_node, *omega_peri, *i_incl, *e_eccent, *T_periastron, 
        *Period, *a_smaxis, *orbit_equinox);
#endif

/* Conversion to radians: */
 *Omega_node *= DEGTORAD;
 *omega_peri *= DEGTORAD;
 *i_incl *= DEGTORAD;
 *mean_motion = (360.0 / *Period) * DEGTORAD;

return(0);
}
/************************************************************************
* compute_residuals_gili
* Scan the input list from the orbit file,
* look for the measures in the tab_calib.tex 
* and compute the residuals 
*
* INPUT:
* fp_in: pointer to the input file containing the measurements and the
*        orbital elements
* fp_out_txt: pointer to the output file with the residuals and various data 
*             in plain ASCII format 
* fp_out_latex: pointer to the output Latex file with the residuals
* fp_out_curve: pointer to the file containing the O-C curve
* fp_out_ref1: pointer to the file with compacted references
* fp_out_ref2: pointer to the file with full references
* iformat: format of the input file (-1=Marco with measures, 
*          1=Marco without measures, 2=OC6 without measures)
*
*************************************************************************/
static int compute_residuals_gili(FILE *fp_in, FILE *fp_out_txt, 
                              FILE *fp_out_latex, FILE *fp_out_curve, 
                              FILE *fp_out_ref1, FILE *fp_out_ref2, 
                              char *calib_fname, 
                              char *OC6_references_fname, int iformat)
{
#define NMAX 1024
double Omega_node, omega_peri, i_incl, e_eccent, T_periastron, orbit_equinox;
double mean_motion, a_smaxis, Period; 
double epoch_o[50], rho_o[50], theta_o[50], err_rho_o[50], err_theta_o[50];
int nmeas, iline, status, is_master_file, line_length, n_names, kk;
char object_name[NMAX*60], discov_name[40], comp_name[10], WDS_name[40]; 
char ADS_name[40], author[NMAX*60], refer0[NMAX*130], refer1[NMAX*130];
/* Maximum line seems to be 265 for OC6 catalog... */
char in_line1[300], in_line2[300];

kk = 0;
iline = 0;

while(!feof(fp_in)) {
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
/* Read input line and retrieve the orbital parameters and the
*  measurements
*/
/* Marco Scardia's format: */
     if(ABS(iformat) == 1) { 
/* Read next line: */
        while(!feof(fp_in)) {
          if(fgets(in_line2,300,fp_in)) {
           iline++;
          if(in_line2[0] != '%') break; 
          }
        }

     status = get_orbit_from_Marco_list(in_line1, in_line2, 
                               iline, &object_name[kk * 60], 
                               WDS_name, ADS_name, discov_name, comp_name, 
                               &author[kk * 60], &Omega_node, 
                               &omega_peri, &i_incl, 
                               &e_eccent, &T_periastron, &Period, &a_smaxis, 
                               &mean_motion, &orbit_equinox);
     object_name[(kk+1)*60 -1] = '\0';
     author[(kk+1)*60 -1] = '\0';
/* OC6 format: */
     } else {
     status = get_orbit_from_OC6_list_gili(in_line1, iline, 
                                 is_master_file, WDS_name, 
                                 discov_name, comp_name, &object_name[kk * 60],
                                 &author[kk * 60], &Omega_node, 
                                 &omega_peri, &i_incl, 
                                 &e_eccent, &T_periastron, &Period, &a_smaxis, 
                                 &mean_motion, &orbit_equinox);
     object_name[(kk+1)*60 -1] = '\0';
     author[(kk+1)*60 -1] = '\0';
     if(!status && *OC6_references_fname) 
          get_OC6_full_reference(&object_name[kk * 60], &author[kk * 60], 
                                 OC6_references_fname, &refer0[kk * 130], 
                                 &refer1[kk * 130]); 
        refer0[(kk+1)*130 -1] = '\0';
        refer1[(kk+1)*130 -1] = '\0';
/* DEBUG:
printf("object=%s author=%s refer0=%s refer1=%s\n", &object_name[kk * 60], 
         &author[kk*60], &refer0[kk * 130], &refer1[kk * 130]);
*/

     }
     if(status) {
     fprintf(stderr, "compute_residuals_gili/WARNING: error reading orbital parameters in line #%d (status=%d)\n", iline, status); 
     } else {
/* Retrieve the measurements from the input files
 and compute the corresponding residuals */
    
/* Retrieve the measurements from the input file (if iformat < 0)
* or retrieve the measurements from the calibrated table.
*/
     if(iformat > 0) {
      printf("UUUTR:From OC6 catalog: iformat=%d object=%s discov=%s comp=>%s<\n", 
              iformat, &object_name[kk * 60], discov_name, comp_name);
/* Read the measurements from the calibrated file if iformat > 0: */
      status = get_measures_from_CALIB_table_gili(calib_fname, 
                       &object_name[kk * 60], comp_name, epoch_o, rho_o, 
                       theta_o, err_rho_o, err_theta_o, &nmeas); 
      printf("UUUTR:From calib table: status=%d epoch_o=%f rho_o=%f theta_o=%f (nmeas=%d)\n", 
             status, epoch_o[0], rho_o[0], theta_o[0], nmeas);
/* Otherwise read the measurements from the input file fp_in: */
     } else {
       status = get_measures_from_file(fp_in, &iline, epoch_o, 
                       rho_o, theta_o, err_rho_o, err_theta_o, &nmeas);
     }
/* Process all measurements */
     if(status < 0) {
     fprintf(stderr, "\ncompute_residuals_gili/Error processing line #%d (%s)\n", 
             iline, in_line1); 
     } else if (status > 0 || nmeas == 0) {
     fprintf(stderr, "compute_residuals_gili/Warning missing measurements in line #%d\n", 
             iline); 
     } else {
/* Process all nmeas measurements relative to a given object:
* correct the measurements from precession, compute the ephemerides,
* and derive the O-C residuals
*/
       status = process_measurements_gili(fp_out_txt, fp_out_latex, 
                                     fp_out_curve, OC6_references_fname, 
                                     &object_name[kk * 60], 
                                     discov_name, comp_name, &author[kk * 60], 
                                     Omega_node, omega_peri, i_incl, e_eccent, 
                                     T_periastron, Period, a_smaxis, 
                                     mean_motion, orbit_equinox, epoch_o, 
                                     rho_o, theta_o, err_rho_o, 
                                     err_theta_o, nmeas);
       if(status) {
       fprintf(stderr, 
               "compute_residuals_gili/Error processing measurements in line #%d\n",
               iline); 
       } else {
/* Update kk the orbit counter: */
       kk++;
       }
     } /* EOF !status (from get_measures_from... */
     } /* EOF !status (succes reading orbital parameters of the current line */
    } /* EOF line is not short (i.e., i > 10) */
    } /* EOF if in_line1 != % */
  } /* EOF if fgets */ 
 }

n_names = kk;
printf("Compute_residuals: %d lines sucessfully read (n_names=%d)\n", 
        iline, n_names);

/* Sort the references by alphabetic order and save them to output
* Latex files: */
if(fp_out_ref1 != NULL && fp_out_ref2 != NULL) 
            sort_references(fp_out_ref1, fp_out_ref2, object_name, author, 
                            refer0, refer1, n_names);
return(0);
}
/************************************************************************
* Process all measurements relative to a given object:
* correct the measurements from precession, compute the ephemerides,
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
static int process_measurements_gili(FILE *fp_out_txt, FILE *fp_out_latex,
              FILE *fp_out_curve, 
              char *OC6_references_fname, char *object_name, char *discov_name,
              char *comp_name, char *author, double Omega_node, 
              double omega_peri, double i_incl, double e_eccent,
              double T_periastron, double Period, double a_smaxis, 
              double mean_motion, double orbit_equinox, double *epoch_o, 
              double *rho_o, double *theta_o, double *err_rho_o, 
              double *err_theta_o, int nmeas)
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

 printf("UUTR/process_measurements/ object=%s nmeas=%d\n ", object_name, nmeas);

/* Main loop on all the measures 
*/
for(i = 0; i < nmeas; i++) {

/* Conversion to radians: */
 theta_o[i] *= DEGTORAD;

/* Compute the ephemerids corresponding to the observation epoch: */
compute_ephemerid(Omega_node, omega_peri, i_incl, e_eccent, T_periastron, 
                  Period, a_smaxis, mean_motion, epoch_o[i], c_tolerance, 
                  &theta_c, &rho_c);

/* Conversion to degrees: */
 theta_o[i] /= DEGTORAD;
 Dtheta = theta_o[i] - theta_c;
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

 Drho = rho_o[i] - rho_c;
/* Trick to have a constant width */
 sprintf(my_name, "%s %s", object_name, comp_name);
/* Left justified text is obtained with a minus sign in the format:*/
 fprintf(fp_out_txt, "%-18.18s %9.3f %9.3f %9.3f %8.2f %8.2f %8.2f %7.1f %s\n",
         my_name, epoch_o[i], rho_o[i], rho_c, Drho, theta_o[i], theta_c, 
         Dtheta, author);
 fprintf(fp_out_latex, "%s %s & %s & %9.3f & %9.3f & %8.2f & %8.2f%s \\\\\n",
         object_name, comp_name, author, epoch_o[i], rho_o[i], Drho, Dtheta,
         quadrant_discrep);
 fprintf(fp_out_curve, "%8.3f %7.2f %8.3f %7.2f %9.3f %-18.18s %s \n",
         Drho, Dtheta, err_rho_o[i], err_theta_o[i], epoch_o[i], my_name, author);
}

return(0);
}
/***************************************************************************
* Read object name in Marco or Luigi's format
*
* Example of input lines:
* 15348+1032 STF1954AB = ADS 9701 - WSI2004a - Mason et al. (2004a)
* 15278+2906 JEF   1 - Tok1984 - Tokovinin (1984)
* or
* ADS 504 - Nov2008a
* JEF   1 - Tok1984
*
* INPUT:
* in_line: line corresponding to the object
* iline: number of "in_line" in the input file
*
* OUTPUT:
* object_name, WDS_name, ADS_name, discov_name, comp_name, author
****************************************************************************/
static int read_object_name1(char *in_line, char *object_name, char *WDS_name,
                             char *ADS_name, char *discov_name, char *comp_name,
                             char *author, int iline)
{
char *pc1, buffer[80];
int i, j, idash, iwds, nval;

/* ADS 11635Cc-D - Doc1984b
*/
 strncpy(object_name, in_line, 40);
 idash = 0;
 for(i = 0; i < 39 && object_name[i]; i++) {
   if(object_name[i] == ' ' && object_name[i+1] == '-') {
     idash = i+1;
     break;
     }
   }
object_name[i] = '\0';

if(idash <= 0) {
  fprintf(stderr, "Process_new_orbit/Bad syntax of object name in line #%d\n (%s)\n", 
          iline, in_line);
  return(-1);
  }

/* Remove the extra-blanks: */
trim_string(object_name, 40);

/* WDS_name (e.g., 15278+2906) */
WDS_name[0] = '\0';
iwds = 0;
if(object_name[5] == '+' || object_name[5] == '-') {
 nval = sscanf(object_name, "%05d%05d", &i, &j);
 if(nval == 2) { 
   strncpy(WDS_name, object_name, 10); 
   WDS_name[10] = '\0';
   iwds = 10;
   }
 }

/* Discover name:  fixed format with 8 characters (e.g. JEF  2 AB) */
strcpy(buffer, &object_name[iwds]);
pc1 = buffer;
i = 0;
while(*pc1  && strncmp(pc1, "ADS", 3) && *pc1 != '=') 
      discov_name[i++] = *(pc1++);
discov_name[i] = '\0';
/* Remove the extra-blanks: */
trim_string(discov_name, 40);

/* Compagnon name (AB, Aa, BC, etc) from discover name 
* (after the numbers) */
comp_name[0] = '\0';
if(*discov_name) {
  pc1 = &discov_name[0];
  while(*pc1  && (isalpha(*pc1) || *pc1 == ' ')) pc1++; 
  while(*pc1  && (isdigit(*pc1) || *pc1 == ' ')) pc1++; 
  strcpy(comp_name, pc1);
  *pc1 = '\0';
}

/* ADS name: */
strcpy(buffer, object_name);
pc1 = buffer;
while(*pc1  && strncmp(pc1, "ADS", 3)) pc1++; 
if(!strncmp(pc1, "ADS", 3))
  strcpy(ADS_name, pc1);
else
  ADS_name[0] = '\0';

/* Remove the extra-blanks: */
trim_string(ADS_name, 40);

/* Compagnon name (AB, Aa, BC, etc) from ADS name 
* (after the numbers) */
if((comp_name[0] == '\0') && (ADS_name[0] != '\0')) {
  pc1 = &ADS_name[0];
  while(*pc1  && (isalpha(*pc1) || *pc1 == ' ')) pc1++; 
  while(*pc1  && (isdigit(*pc1) || *pc1 == ' ')) pc1++; 
  strcpy(comp_name, pc1);
  *pc1 = '\0';
}

/* Copy the author (removing the Carriage Return/EOF if present): */
i = 0;
pc1 = &in_line[idash + 1];
while(*pc1 && i < 59) {
   if(isprint(*pc1)) author[i++] = *pc1;
   pc1++;
   }
author[i] = '\0';

/* Remove the extra-blanks: */
trim_string(author, 60);

#ifdef DEBUG
printf("iline=%d\n object: >%s< author: >%s<\n", iline, object_name, author);
printf("WDS=%s, discov=%s, (comp=%s) ADS=%s. \n", 
        WDS_name, discov_name, comp_name, ADS_name);
#endif

/* Restriction of the object name to ADS name or discoverer name */
if(ADS_name[0] != '\0') strcpy(object_name, ADS_name);
else if (discov_name[0] != '\0') strcpy(object_name, discov_name);

return(0);
}

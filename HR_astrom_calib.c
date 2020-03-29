/*************************************************************************
* Program HR_astrom_calib
* - to generate a table with the results of all the observations
* made in Merate
* - to compute statistics and HR diagram 
* 
* Read all the files such as: astrom05a_PaperIII.tex, astrom05b.tex, etc...
*
* Format of input files:
\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}
& Name  & Epoch & Filter & Eyep. & $\rho$ & $\sigma_\rho$ & $\theta$ &
$\sigma_\theta$ & comments \\
 &       &       &        & (mm)
& (pixels) & (pixels) & (deg.) & (deg.) & \\
\hline % -----------------------------------------------------------------------
%%
16564+6502 = STF 2118 AB & ADS 10279 & 2004. & & & & & & & orb \\
%%
%% 090904_ads10279_Vd_a
%% Thresholds: -1,8
%% rho_center=14.62 (+/-1.87) theta_center=156.71 (+/-7.32)or theta=-23.29
%% rho_center=14.72 theta_center=156.70  or theta = -23.30
%% rho_center=14.37 (+/-1.89) theta_center=-23.64 (+/-7.54)
%% rho_center=14.34 theta_center=-23.73
%% mean: rho = 14.50 \pm 0.17  theta = -23.61 \pm 0.27 (too small -> 0.3)(n=8)
& 090904\_ads10279\_Vd\ & 09/09/2004 & V & 20 & 14.50 & 0.17 & -23.61 & 0.3 & \\
....
*
Keywords in the notes:
EP=2004.133 (epoch)
WR=127      (Rho, WDS)
WT=127      (Theta, WDS)
WY=2002     (Year of WDS observation) 
Q=2         (Quadrant with restricted triple-correplation)
LQ=3        (Quadrant with long integration)
*
* JLP
* Version 26/06/2009
*************************************************************************/
#include "jlp_trim.h"  /* (in jlplib/jlp_fits/ ) trim_string... */
#include "astrom_utils1.h"
#include "astrom_utils2.h"

/*
#define DEBUG
*/

static int HR_calib(FILE *fp_list, FILE *fp_out, double scale_10, 
                    double scale_20, double scale_32, double theta0, 
                    double sign, int i_filename, int i_date, int i_filter,
                    int i_eyepiece, int i_rho, int i_drho, int i_theta, 
                    int i_dtheta, int i_notes, int comments_wanted, 
                    char *file_with_list);

int main(int argc, char *argv[])
{
double scale_10, scale_20, scale_32, theta0, sign;
char file_with_list[80], fileout[80];
int i_filename, i_date, i_filter, i_eyepiece; 
int i_rho, i_drho, i_theta, i_dtheta; 
int i_notes, comments_wanted;
FILE *fp_list, *fp_out;

/* If command line with "runs" */
if(argc == 7){
 if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
 else argc = 1;
 }


if(argc != 3)
  {
  printf("Syntax: HR_astrom_calib file_with_list out_file \n");
  printf("Example: runs HR_astrom_calib astrom_list.txt HR_calib_full.tex \n\n");
  printf("Example of contents of astrom_list.txt: \n\n");
  printf("astrom05a_PaperIII.tex\n");
  printf("astrom05b_PaperIV.tex\n");
  printf("astrom06a_PaperV.tex\n");
  exit(-1);
  }
else
  {
  strcpy(file_with_list, argv[1]);
  strcpy(fileout, argv[2]);
  }

printf(" OK: file_with_list=%s fileout=%s \n", file_with_list, fileout);

/* out_rho = in_rho * scale_10  (or in_rho * scale_20 if 20 mm eyepiece)
   out_theta = in_theta * sign + theta0 
*/
/* Calibration values obtained with the slit mask in 2005-2006: 
  scale_10 = 0.0320;
  scale_20 = 0.0752;
  scale_32 = 0.158;
  sign = 1.;
  theta0 = 89.77;
*/
/* Calibration values obtained with the slit mask in 2011: */
  scale_10 = 0.0320;
  scale_20 = 0.0754;
  scale_32 = 0.158;
  sign = 1.;
  theta0 = 89.94;
  printf(" out_rho = in_rho * %g (for 10 mm eyepiece)\n", scale_10);
  printf(" out_rho = in_rho * %g (for 20 mm eyepiece)\n", scale_20);
  printf(" out_rho = in_rho * %g (for 32 mm eyepiece)\n", scale_32);
  printf(" out_theta = in_theta * %g + %g \n", sign, theta0);


if((fp_list = fopen(file_with_list,"r")) == NULL) {
  printf(" Fatal error opening input list file %s \n", file_with_list);
  exit(-1);
  }


if((fp_out = fopen(fileout,"w")) == NULL) {
  printf(" Fatal error opening output file %s \n", fileout);
  fclose(fp_list);
  exit(-1);
  }

fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(fp_out,"%%%% Calibration with HR_astrom of files contained in %s\n", file_with_list);
fprintf(fp_out,"%%%% JLP / Version of 05/07/2012 \n");
fprintf(fp_out,"%%%% out_rho = in_rho * %g (for 10 mm eyepiece)\n", scale_10);
fprintf(fp_out,"%%%% out_rho = in_rho * %g (for 20 mm eyepiece)\n", scale_20);
fprintf(fp_out,"%%%% out_rho = in_rho * %g (for 32 mm eyepiece)\n", scale_32);
fprintf(fp_out,"%%%% out_theta = in_theta * %g + %g \n", sign, theta0);
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

fprintf(fp_out,"%%%% File automatically generated with HR_astrom (Version 26/06/2009)\n");
fprintf(fp_out,"%%%% Convention for Papers II, III & IV (Merate): \n");
fprintf(fp_out,"%%%% Publication mode: only direct measurement is kept when rho > 0.3\" \n");
fprintf(fp_out,"%%%% Weighted mean between direct and recorded measurements when rho <= 0.3\" \n");
fprintf(fp_out,"%%%% Minimum error for rho: 0.1 pixel  or 0.5%% \n");
fprintf(fp_out,"%%%% Minimum error for theta: 0.3 degree \n");

/* Scan the file and make the conversion: */
/* date in column 3 
*  eyepiece in column 5
*  rho in column 6
*  drho in column 7
*  theta in column 8
*  dtheta in column 9
*/
  i_filename = 2;
  i_date = 3;
  i_filter = 4;
  i_eyepiece = 5;
  i_rho = 6;
  i_drho = 7;
  i_theta = 8;
  i_dtheta = 9;
  i_notes = 10;
/* To generate a laTeX array formatted for publication:
* sort out the objects, compute mean values, perform calibration, etc.
*/ 
  HR_calib(fp_list,fp_out, scale_10, scale_20, scale_32, theta0, sign, 
           i_filename, i_date, i_filter, i_eyepiece, i_rho, i_drho, 
           i_theta, i_dtheta, i_notes, comments_wanted, file_with_list);

/* Version "publi", to generate a laTeX array formatted for publication:
* sort out the objects, compute mean values, perform calibration, etc.
  astrom_latex_calib_publi(fp_in, fp_out, scale_10, scale_20, scale_32, theta0,
                           sign, i_filename, i_date, i_filter, i_eyepiece, 
                           i_rho, i_drho, i_theta, i_dtheta, i_notes, 
                           comments_wanted, filein);
*/

fclose(fp_list);
fclose(fp_out);
return(0);
}
/*************************************************************************
* Calibration of astrom files
*
* INPUT:
* i_filename: column nber of the filename used for this measurement
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
* comments_wanted: 1 if comments are wanted, 0 otherwise
*************************************************************************/
static int HR_calib(FILE *fp_list, FILE *fp_out, double scale_10, 
                    double scale_20, double scale_32, double theta0, 
                    double sign, int i_filename, int i_date, int i_filter,
                    int i_eyepiece, int i_rho, int i_drho, int i_theta, 
                    int i_dtheta, int i_notes, int comments_wanted, 
                    char *file_with_list)
{
OBJECT *obj;
int *index_obj, nobj, tabular_only; 
char infile_name[80];
FILE *fp_in;

if((obj = (OBJECT *)malloc(NOBJ_MAX * sizeof(OBJECT))) == NULL) {
  printf("Fatal error allocating memory space for OBJECT: nobj_max=%d\n", 
          NOBJ_MAX);
  exit(-1);
  }
if((index_obj = (int *)malloc((NOBJ_MAX) * sizeof(int))) == NULL) {
  printf("Fatal error allocating memory space for index_obj: nobj_max=%d\n", 
          NOBJ_MAX);
  exit(-1);
  }

nobj = 0;

while(!feof(fp_list))
{
/* Reading the name of the next "astrom" file to be read: */
  if(fgets(infile_name, 80, fp_list)){
/* Possibility of comments, starting with % or #: */
    if(infile_name[0] != '#' && infile_name[0] != '%') {
/* Removing blanks, if needed: */
    trim_string(infile_name, 80);
/* Opening this "astrom" file: */
    if((fp_in = fopen(infile_name,"r")) == NULL) {
      printf(" Fatal error opening input astrom file >%s< \n", infile_name);
      exit(-1);
      }

    astrom_read_measures(fp_in, comments_wanted, obj, &nobj, i_filename, 
                         i_date, i_filter, i_eyepiece, i_rho, i_drho, i_theta, 
                         i_dtheta, i_notes);
printf("After %s : total value of nobj = %d\n", infile_name, nobj);
    fclose(fp_in);
    } /* EOF infile[0] != % */
    } /* EOF fgets ... */
} /* EOF while !feof */

#ifdef DEBUG
printf("Returned by astrom_read_measures: nobj = %d\n", nobj);
#endif

astrom_calibrate_measures(obj, nobj, scale_10, scale_20, scale_32,
                          theta0, sign);

/* Sort the objects according to their Right Ascension: */
astrom_ra_sort_objects(obj, index_obj, nobj);

/* Compute mean values for the objects with rho < 0.3" 
*  and discard all measurements on recorded data if rho >= 0.3" (Merate-Paper II) */
/*
astrom_mean_for_paper2(obj, nobj);
*/
astrom_mean_for_full_table(obj, nobj);

/* For big tables, should set this parameter to 1: */
tabular_only = 1; 
astrom_write_publi_table(fp_out, comments_wanted, obj, index_obj, nobj,
                         tabular_only);

astrom_compute_statistics(fp_out, obj, nobj, infile_name);

free(index_obj);
free(obj);
return(0);
}

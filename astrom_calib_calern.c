/*************************************************************************
* Program astrom_calib_calern
* (formely called latex_calib)
* To convert a Latex array with measurements in pixels to array in arcseconds
* (and conversion of angles from XY axis to Noth-East referential)
*
* derived from astromp_calib, with more possibilities for the eyepieces
* with a parameter file
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
* Version 02/07/2018
*************************************************************************/
#include "astrom_utils1.h" 
#include "astrom_utils2.h" 

static int read_calib_file(char *filecalib, double *calib_scale1, 
                           int *calib_eyepiece1, int *n_eyepieces1,
                           double *sign, double *theta0);

int main(int argc, char *argv[])
{
double calib_scale1[3], theta0, sign;
int i, calib_eyepiece1[3], n_eyepieces1 = 3; 
char filein[128], fileout[128], filecalib[128];
int i_filename, i_date, i_filter, i_eyepiece, i_rho, i_drho, i_theta, i_dtheta; 
int i_notes, comments_wanted, publi_mode, input_with_header, status;
int gili_format;
FILE *fp_in, *fp_out;

/* If command line with "runs" */
if(argc == 7){
 if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
 else argc = 1;
 }


if(argc != 5)
  {
  printf(" Syntax: astrom_calib calibration_file in_file out_file comments_wanted,publi_mode,input_with_header \n");
  printf(" out_rho = in_rho * scale_10  (or in_rho * scale_20 if 20 mm eyepiece)\n");
  printf(" out_theta = in_theta * sign + theta0 \n");
/*
* publi mode : 1,1,0 with WDS enriched astrom file
* copy mode: 1,0,0 with raw astrom file
*/
  printf(" publi_mode = 0 (calibrated, with no other modifications)\n");
  printf(" publi_mode = 1 (calibrated and ready for publication)\n");
  printf(" input_with_header = 0 (input LaTeX file without header)\n");
  printf(" input_with_header = 1 (input LaTeX file with header)\n");
  printf(" Example: runs astrom_calib 0.03202,0.07517,1.,89.77 astrom05a.tex tab_calib.tex 0,1,0  (without comments, publication mode, input file without header)\n");
  printf(" Example: runs astrom_calib 0.03202,0.07517,1.,89.77 astrom05a.tex tab_calib.tex 1,0,1  (with comments, no publication, input file with header)\n");
  exit(-1);
  }
else
  {
  strcpy(filecalib,argv[1]);
  strcpy(filein,argv[2]);
  strcpy(fileout,argv[3]);
  sscanf(argv[4],"%d,%d,%d", &comments_wanted, &publi_mode, &input_with_header);
  }

// Read calibration file:
 status = read_calib_file(filecalib, calib_scale1, calib_eyepiece1, 
                          &n_eyepieces1, &sign, &theta0);
 if(status != 0) return(-1);
for(i = 0; i < n_eyepieces1; i++) {
  printf(" OK: out_rho = in_rho * %g (for %d mm eyepiece)\n", 
         calib_scale1[i], calib_eyepiece1[i]);
  }
printf(" OK: out_theta = in_theta * %g + %g \n", sign, theta0);
printf(" OK: filein=%s fileout=%s \n", filein, fileout);
printf(" OK: comments_wanted=%d publi_mode=%d input_with_header=%d\n",
       comments_wanted, publi_mode, input_with_header);

if((fp_in = fopen(filein,"r")) == NULL)
{
printf(" Fatal error opening input file %s \n",filein);
exit(-1);
}

if((fp_out = fopen(fileout,"w")) == NULL)
{
printf(" Fatal error opening output file %s \n",fileout);
fclose(fp_in);
exit(-1);
}
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(fp_out,"%%%% Automatic calibration of %s with astrom_calib\n", filein);
fprintf(fp_out,"%%%% JLP / Version of 23/09/2008 \n");
for(i = 0; i < n_eyepieces1; i++) {
  fprintf(fp_out,"%%%% out_rho = in_rho * %g (for %d mm eyepiece)\n", 
         calib_scale1[i], calib_eyepiece1[i]);
  }
fprintf(fp_out,"%%%% out_theta = in_theta * %g + %g \n", sign, theta0);
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

if(publi_mode == 1) {
  fprintf(fp_out,"%%%% File automatically generated with astrom_calib (Version 12/11/2009)\n");
  fprintf(fp_out,"%%%% Convention for Papers II, III & IV (Merate): \n");
  fprintf(fp_out,"%%%% Publication mode: only direct measurement is kept when rho > 0.3\" \n");
  fprintf(fp_out,"%%%% Weighted mean between direct and recorded measurements when rho <= 0.3\" \n");
  fprintf(fp_out,"%%%% Minimum error for rho: 0.1 pixel  or 0.5%% \n");
  fprintf(fp_out,"%%%% Minimum error for theta: 0.3 degree \n");
  } else {
    if(input_with_header == 0) {
    fprintf(fp_out,"%%%% Header was automatically added since it was missing in original file\n");
    fprintf(fp_out,"%%%% (No other modifications since publi_mode = %d)\n", publi_mode);
    }
  }

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
/* Version "publi", to generate a laTeX array formatted for publication:
* sort out the objects, compute mean values, perform calibration, etc.
*
* publi mode : 1,1,0 with WDS enriched astrom file
* copy mode: 1,0,0 with raw astrom file
*/ 
  if(publi_mode == 1) {
  gili_format = 0;
  astrom_calib_publi(fp_in,fp_out, calib_scale1, calib_eyepiece1, n_eyepieces1,
                     theta0, sign, 
                     i_filename, i_date, i_filter, i_eyepiece, i_rho, i_drho, 
                     i_theta, i_dtheta, i_notes, comments_wanted, filein,
                     gili_format);
/* Version "copy", to copy all the contents, and simply convert the 
* measurements to arcseconds and degrees relative to North 
*/ 
  } else { 
  astrom_calib_copy(fp_in,fp_out, calib_scale1, calib_eyepiece1, n_eyepieces1,
                     theta0, sign, 
                    i_date, i_eyepiece, i_rho, i_drho, i_theta, i_dtheta,
                    comments_wanted, input_with_header);
  }

fclose(fp_in);
fclose(fp_out);
return(0);
}
/**************************************************************************
* Read calibration file:
*
* Example:

# Calern November 2016 (Andor):
# 10mm: Mean scale = 0.0180 +/- 0.0014
# 20mm: Mean scale = 0.0424 +/- 0.0033
# 32mm: Mean scale = 0.0899 +/- 0.0066
# scale_10, scale_20, scale_32, sign, theta0
10mm = 0.0180
20mm = 0.0424
32mm = 0.0899
sign = 1.
theta0 = 89.0

***************************************************************************/
static int read_calib_file(char *filecalib, double *calib_scale1, 
                           int *calib_eyepiece1, int *n_eyepieces1,
                           double *sign, double *theta0)
{
FILE *fp_calib;
char buffer[256];
double scale0;
int status = -1, ival, jval;

 if((fp_calib = fopen(filecalib, "r")) == NULL) {
   fprintf(stderr, "read_calib_file/Error opening file: >%s<\n",
                filecalib); 
   return(status);
   }

*sign = 0.;
*theta0 = 0.;

ival = 0;
jval = 0;
 while(!feof(fp_calib)) {
   if(fgets(buffer, 256, fp_calib)) {
// Comments start with % or #
     if((buffer[0] != '%') && (buffer[0] != '#')) {
       if(!strncmp(buffer, "10mm =", 6)) {
         sscanf(&buffer[6], "%lf", &scale0);
         calib_scale1[ival] = scale0;
         calib_eyepiece1[ival] = 10;
         ival++;
       } else if(!strncmp(buffer, "20mm =", 6)) {
         sscanf(&buffer[6], "%lf", &scale0);
         calib_scale1[ival] = scale0;
         calib_eyepiece1[ival] = 20;
         ival++;
       } else if(!strncmp(buffer, "32mm =", 6)) {
         sscanf(&buffer[6], "%lf", &scale0);
         calib_scale1[ival] = scale0;
         calib_eyepiece1[ival] = 32;
         ival++;
       } else if(!strncmp(buffer, "sign =", 6)) {
         sscanf(&buffer[6], "%lf", sign);
         jval++;
       } else if(!strncmp(buffer, "theta0 =", 8)) {
         sscanf(&buffer[8], "%lf", theta0);
         jval++;
       } else {
         fprintf(stderr, "read_calib_file/Bad syntax in line=>%s<\n",
                buffer);
         }
     if(jval == 2) {
       status = 0;
       break;
       }
     } // EOF if buffer...
   } // EOF if fgets...
 } // EOF while(!feof...)

*n_eyepieces1 = ival;

 fclose(fp_calib);

return(status);
}

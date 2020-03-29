/*************************************************************************
* Program astrom_calib
* (formely called latex_calib)
* To convert a Latex array with measurements in pixels to array in arcseconds
* (and conversion of angles from XY axis to Noth-East referential)
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
* Version 02/07/2012
*************************************************************************/
#include "astrom_utils1.h" 
#include "astrom_utils2.h" 

int main(int argc, char *argv[])
{
double calib_scale1[3], scale_10, scale_20, scale_32, theta0, sign;
int calib_eyepiece1[3];
char filein[60], fileout[60];
int i_filename, i_date, i_filter, i_eyepiece, i_rho, i_drho, i_theta, i_dtheta; 
int i_notes, comments_wanted, publi_mode, input_with_header;
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
  printf(" Syntax: astrom_calib scale_10,scale_20,scale_32,sign,theta0 in_file out_file comments_wanted,publi_mode,input_with_header \n");
  printf(" out_rho = in_rho * scale_10  (or in_rho * scale_20 if 20 mm eyepiece)\n");
  printf(" out_theta = in_theta * sign + theta0 \n");
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
  sscanf(argv[1],"%lf,%lf,%lf,%lf,%lf", &scale_10, &scale_20, &scale_32, &sign, &theta0);
  strcpy(filein,argv[2]);
  strcpy(fileout,argv[3]);
  sscanf(argv[4],"%d,%d,%d", &comments_wanted, &publi_mode, &input_with_header);
  }
printf(" OK: out_rho = in_rho * %g (for 10 mm eyepiece)\n", scale_10);
printf(" OK: out_rho = in_rho * %g (for 20 mm eyepiece)\n", scale_20);
printf(" OK: out_rho = in_rho * %g (for 32 mm eyepiece)\n", scale_32);
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
fprintf(fp_out,"%%%% out_rho = in_rho * %g (for 10 mm eyepiece)\n", scale_10);
fprintf(fp_out,"%%%% out_rho = in_rho * %g (for 20 mm eyepiece)\n", scale_20);
fprintf(fp_out,"%%%% out_rho = in_rho * %g (for 32 mm eyepiece)\n", scale_32);
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
*/ 
  calib_eyepiece1[0] = 10;
  calib_eyepiece1[1] = 20;
  calib_eyepiece1[2] = 32;
  calib_scale1[0] = scale_10;
  calib_scale1[1] = scale_20;
  calib_scale1[2] = scale_32;
  if(publi_mode == 1) {
// n_eyepieces=3 here
  gili_format = 0;
  astrom_calib_publi(fp_in,fp_out, calib_scale1, calib_eyepiece1, 3,
                     theta0, sign, 
                     i_filename, i_date, i_filter, i_eyepiece, i_rho, i_drho, 
                     i_theta, i_dtheta, i_notes, comments_wanted, filein,
                     gili_format);
/* Version "copy", to copy all the contents, and simply convert the 
* measurements to arcseconds and degrees relative to North 
*/ 
  } else { 
// n_eyepieces=3 here
  astrom_calib_copy(fp_in,fp_out, calib_scale1, calib_eyepiece1, 3,
                     theta0, sign, 
                    i_date, i_eyepiece, i_rho, i_drho, i_theta, i_dtheta,
                    comments_wanted, input_with_header);
  }

fclose(fp_in);
fclose(fp_out);
return(0);
}

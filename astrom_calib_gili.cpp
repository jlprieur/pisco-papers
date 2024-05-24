/*************************************************************************
* Program astrom_calib_gili
* (formely called latex_calib)
* To convert a Latex array with measurements in pixels to array in arcseconds
* (and conversion of angles from XY axis to Noth-East referential)
*
* derived from astromp_calib, with more possibilities for the eyepieces
* with a parameter file
*
* Gili's format of input file, from Gdpisco, old version (in_astrom_fmt=1) :
\hline % -------------------------------------------------------------
14595+1753 = COU188 &  & 2013 & & & & & & &  WY=2010 WT=227 WR=0.3 \\
%%
%% Input file: cou188_a.fits
%% Back to original image
....
\hline % -------------------------------------------------------------
15186+2356 = COU307 &  & 2013 & & & & & & & orb WY=2010 WT=3 WR=0.4 \\
*
* Calern format of input file, from Gdpisco, new version (in_astrom_fmt = 2):
\hline % -------------------------------------------------------------
00226+5417 = A907 & 2017 & & & & & & & WY=2016 WT=211 WR=0.8 \\
%%
%% Input file: 090117_a907_R_a.fits
%% rho=19.18 theta=-56.77 xc=117.49 yc=144.04 (117.6,144.2,7.5,3) Bary.
...
& 090117\_a907\_R\_a & 09/01/2017 & R & 20 & 19.16 & 0.10 & -56.88 & 0.3 & EP=2017.0264  Dm=0.00+/-0.02\\
*
* The only difference of the astrom input files is an empty (ADS) column 
* for Gili's files and the orb information on the last column for Calern files
*
* in_astrom_fmt = 1 (Gili, with WDS=name, Empty ADS col, Epoch, ...)
* in_astrom_fmt = 2 (Calern, with WDS=name, Epoch, ...)               
* out_calib_fmt = 1 (Gili, with WDS, Name, Epoch, Eyepiece, rho, err_rho, ...)
* out_calib_fmt = 2 (Calern, with WDS, Name, Epoch, Filter, Eyepiece, rho, err_r
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
* Version 15/09/2021
*************************************************************************/
#include "astrom_utils1.h" 
#include "astrom_utils2.h" 
#include "jlp_fitsio.h"  // BesselToJulian

static int read_calib_file(char *filecalib, double *calib_scale1, 
                           int *calib_eyepiece1, int*n_eyepieces1,
                           int *sign1, double *theta01, double *calib_date1,
                           int *calib_dd1, int *calib_mm1, int *calib_year1,
                           int *ncalib1, int ndim);

/*************************************************************************
*
*************************************************************************/
#define NDIM 16
int main(int argc, char *argv[])
{
double calib_scale1[NDIM * 8], theta01[NDIM], calib_date1[NDIM];
int sign1[NDIM]; 
int i, icalib, calib_eyepiece1[NDIM * 8], n_eyepieces1[NDIM], ncalib1; 
int calib_dd1[NDIM], calib_mm1[NDIM], calib_year1[NDIM];
int ndim = NDIM, nval;
char filein[128], fileout[128], filecalib[128];
int i_filename, i_date, i_filter, i_eyepiece, i_rho, i_drho, i_theta, i_dtheta; 
int i_notes, comments_wanted, publi_mode, input_with_header, status;
int in_astrom_fmt = 2, out_calib_fmt = 2;
FILE *fp_in, *fp_out;

/* If command line with "runs" */
if(argc == 7){
 if(*argv[5]) argc = 6;
 else if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
 else argc = 1;
 }

if(argc != 6)
  {
  printf(" Syntax: astrom_calib calibration_file in_file out_file comments_wanted,publi_mode,input_with_header \n");
  printf(" out_rho = in_rho * scale_bin1  (or in_rho * 2 * scale_bin1 if bin=2)\n");
  printf(" out_theta = in_theta * sign1 + theta01 \n");
/*
* publi mode : 1,1,0 with WDS enriched astrom file
* copy mode: 1,0,0 with raw astrom file
*/
  printf(" publi_mode = 0 (calibrated, with no other modifications)\n");
  printf(" publi_mode = 1 (calibrated and ready for publication)\n");
  printf(" input_with_header = 0 (input LaTeX file without header)\n");
  printf(" input_with_header = 1 (input LaTeX file with header)\n");
  printf(" Example: runs astrom_calib gilicalib.txt astrom05a.tex tab_calib.tex 0,1,0  (without comments, publication mode, input file without header) in_astrom_fmt,out_calib_fmt \n");
  printf(" Example: runs astrom_calib gilicalib.txt astrom05a.tex tab_calib.tex 1,0,1  (with comments, no publication, input file with header) 1 \n");
  printf(" in_astrom_fmt: 1 if (Gili) empty ADS column, 2 if (Calern) no ADS column\n");
  printf(" out_calib_fmt: 1 if (Gili) no filt column, 2 if (Calern) filter columnbut\n");
  exit(-1);
  }
else
  {
  strcpy(filecalib,argv[1]);
  strcpy(filein,argv[2]);
  strcpy(fileout,argv[3]);
  nval = sscanf(argv[4],"%d,%d,%d", &comments_wanted, &publi_mode, &input_with_header);
  if(nval != 3) {
     fprintf(stderr, "Bad syntax: nval=%d (!= 3)\n", nval);
     exit(-1);
     }
  nval = sscanf(argv[5],"%d,%d", &in_astrom_fmt, &out_calib_fmt);
  if((nval != 2) || ((in_astrom_fmt != 1) && (in_astrom_fmt != 2))
       || ((out_calib_fmt != 1) && (out_calib_fmt != 2)) ) {
     fprintf(stderr, "Bad syntax: nval=%d (!= 2 ?) in_astrom_fmt=%d out_calib_fmt=%d\n", 
         nval, in_astrom_fmt, out_calib_fmt);
     exit(-1);
     }
  }

printf(" OK: filein=%s fileout=%s \n", filein, fileout);
printf(" OK: comments_wanted=%d publi_mode=%d input_with_header=%d in_astrom_fmt=%d out_calib_fmt=%d\n",
       comments_wanted, publi_mode, input_with_header, in_astrom_fmt, out_calib_fmt);

// Read calibration file:
status = read_calib_file(filecalib, calib_scale1, calib_eyepiece1, 
                          n_eyepieces1, sign1, theta01, calib_date1, 
                          calib_dd1, calib_mm1, calib_year1,
                          &ncalib1, ndim);
if(status != 0) return(-1);
for(icalib = 0; icalib < ncalib1; icalib++) {
  printf(" calib #%d\n calib_date1=%.4f (%02d/%02d/%d)\n",
          icalib, calib_date1[icalib], calib_dd1[icalib], 
          calib_mm1[icalib], calib_year1[icalib]);
  for(i = 0; i < n_eyepieces1[icalib]; i++) {
  printf(" out_rho = in_rho * %g (for %d mm eyepiece)\n",
         calib_scale1[ndim * icalib + i], calib_eyepiece1[ndim * icalib + i]);
  }
  printf(" out_theta = in_theta * %d + %g \n", sign1[icalib], theta01[icalib]);
}

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
fprintf(fp_out,"%%%% JLP / Version of 23/09/2020 \n");
for(icalib = 0; icalib < ncalib1; icalib++) {
fprintf(fp_out,"%%%% calib #%d calib_date1=%.4f (%02d/%02d/%d)\n",
        icalib, calib_date1[icalib], calib_dd1[icalib], 
        calib_mm1[icalib], calib_year1[icalib]);
for(i = 0; i < n_eyepieces1[icalib]; i++) {
  fprintf(fp_out,"%%%% out_rho = in_rho * %g (for eyepiece/binning=%d)\n", 
         calib_scale1[icalib * ndim + i], calib_eyepiece1[icalib * ndim + i]);
  }
fprintf(fp_out,"%%%% out_theta = in_theta * %d + %g \n", 
        sign1[icalib], theta01[icalib]);
}
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

if(publi_mode == 1) {
  fprintf(fp_out,"%%%% File automatically generated with astrom_calib_gili (Version 16/09/2020)\n");
  fprintf(fp_out,"%%%% Convention for Papers of R. Gili (Nice): \n");
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
// Example of input line of measures (output line created by Gdpisco.exe):
// mean: rho = 27.78 \pm 0.00 (too small->0.1) theta = 22.85 \pm 0.00 (too small->0.3) (n=1)
// mean: Dmag = 4.66 \pm 0.00 (too small->0.02) (n=1)
// & rst5109ly\_c10w & 17/08/2013 &  & 1 & 27.78 & 0.10 & 22.85 & 0.3 & EP=2013.6290  Dm=4.66+/-0.02\\\
  i_filename = 2; // ex: "a153\aw"
  i_date = 3; // ex: "28/08/2013"
  i_filter = 4; // ex: ""
  i_eyepiece = 5; // bin factor, ex: "1"
  i_rho = 6; // ex: "27.78"
  i_drho = 7; // ex: "0.10"
  i_theta = 8; // ex: "22.85"
  i_dtheta = 9; // ex: "0.3"
  i_notes = 10; // ex: EP=2013.5567 Dm=4.66 XYBIN=2,2 Q=3
/* Version "publi", to generate a laTeX array formatted for publication:
* sort out the objects, compute mean values, perform calibration, etc.
*
* publi mode : 1,1,0 with WDS enriched astrom file
* copy mode: 1,0,0 with raw astrom file
*/ 
  if(publi_mode == 1) {
/*
* in_astrom_fmt : =1 if (Gili) empty ADS column, no orb info
* 22388+4419 = HO 295 AB &  & 2004. & & & & & & & \\
*                 =2 if (Calern) no ADS column, no orb info
* 22388+4419 = HO 295 AB & 2004. & & & & & & & \\
* in_astrom_fmt = 1 (Gili, with WDS=name, Empty ADS col, Epoch, ...)
* in_astrom_fmt = 2 (Calern, with WDS=name, Epoch, ...)               
* out_calib_fmt = 1 (Gili, with WDS, Name, Epoch, Eyepiece, rho, ..., Dmag)
* out_calib_fmt = 2 (Calern, with WDS, Name, Epoch, Filter, Eyepiece, rho, err_r
*/
// gili_format: in_astrom_fmt=1
// calern_format: in_astrom_fmt=2
  printf("Will run astrom_calib_publi with in_astrom_fmt=%d\n", in_astrom_fmt);
  astrom_calib_publi(fp_in, fp_out, 
                     calib_date1, calib_dd1, calib_mm1, calib_year1,
                     calib_scale1, calib_eyepiece1, n_eyepieces1,
                     theta01, sign1, ncalib1, ndim, 
                     i_filename, i_date, i_filter, i_eyepiece, i_rho, i_drho, 
                     i_theta, i_dtheta, i_notes, comments_wanted, filein,
                     in_astrom_fmt, out_calib_fmt);
/* Version "copy", to copy all the contents, and simply convert the 
* measurements to arcseconds and degrees relative to North 
*/ 
  } else { 
  astrom_calib_copy(fp_in, fp_out, 
                    calib_date1, calib_dd1, calib_mm1, calib_year1,
                    calib_scale1, calib_eyepiece1, n_eyepieces1,
                    theta01, sign1, ncalib1, ndim, 
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
# Nice November 2011 (Andor):
bin1 = 0.0738
bin2 = 0.1476
bin4 = 0.2952
sign = 1
theta0 = 89.0

***************************************************************************/
static int read_calib_file(char *filecalib, double *calib_scale1, 
                           int *calib_eyepiece1, int*n_eyepieces1,
                           int *sign1, double *theta01, double *calib_date1,
                           int *calib_dd1, int *calib_mm1, int *calib_year1,
                           int *ncalib1, int ndim)
{
FILE *fp_calib;
char buffer[256];
int dd0, mm0, yy0, sign;
double scale0, theta0, year0, time0, JulianDate0;
int status = -1, i, ival, jval, icalib;

 if((fp_calib = fopen(filecalib, "r")) == NULL) {
   fprintf(stderr, "read_calib_file/Error opening file: >%s<\n",
                filecalib); 
   return(status);
   }

icalib = -1;
ival = 0;
jval = 0;
 while(!feof(fp_calib)) {
   if(fgets(buffer, 256, fp_calib)) {
// Comments start with % or #
     if((buffer[0] != '%') && (buffer[0] != '#')) {
//     printf("ZZZ %s", buffer);
       if(buffer[0] == '=') {
         if(icalib >= 0) { 
           n_eyepieces1[icalib] = ival;
           if(jval < 2) {
            fprintf(stderr, "read_calib_file/Fatal error loading icalib=%d jval=%d\n", 
                   icalib, jval);
            exit(-1);
           }
         }
         icalib++;
         ival = 0;
         jval = 0;
         sscanf(&buffer[1], "%d/%d/%d", &dd0, &mm0, &yy0);
         year0 = (double)yy0;
         time0 = 0.;
// int JLP_julian_epoch(double aa, int mm, int idd, double time, double *j_date);
         JLP_julian_epoch(year0, mm0, dd0, time0, &JulianDate0);
         calib_dd1[icalib] = dd0;
         calib_mm1[icalib] = mm0;
         calib_year1[icalib] = year0;
         calib_date1[icalib] = JulianDate0;
         printf("ZZZ icalib=%d calib_date1=%f %02d/%02d/%d\n", 
                 icalib, calib_date1[icalib], calib_dd1[icalib], 
                 calib_mm1[icalib], calib_year1[icalib]);
         }
       if(icalib >= 0) {
       if(!strncmp(buffer, "10mm =", 6)) {
         sscanf(&buffer[6], "%lf", &scale0);
         calib_scale1[icalib * ndim + ival] = scale0;
         calib_eyepiece1[icalib * ndim + ival] = 10;
         ival++;
       } else if(!strncmp(buffer, "15mm =", 6)) {
         sscanf(&buffer[6], "%lf", &scale0);
         calib_scale1[icalib * ndim + ival] = scale0;
         calib_eyepiece1[icalib * ndim + ival] = 15;
         ival++;
       } else if(!strncmp(buffer, "20mm =", 6)) {
         sscanf(&buffer[6], "%lf", &scale0);
         calib_scale1[icalib * ndim + ival] = scale0;
         calib_eyepiece1[icalib * ndim + ival] = 20;
         ival++;
       } else if(!strncmp(buffer, "24mm =", 6)) {
         sscanf(&buffer[6], "%lf", &scale0);
         calib_scale1[icalib * ndim + ival] = scale0;
         calib_eyepiece1[icalib * ndim + ival] = 24;
         ival++;
       } else if(!strncmp(buffer, "32mm =", 6)) {
         sscanf(&buffer[6], "%lf", &scale0);
         calib_scale1[icalib * ndim + ival] = scale0;
         calib_eyepiece1[icalib * ndim + ival] = 32;
         ival++;
       } else if(!strncmp(buffer, "bin1 =", 6)) {
         sscanf(&buffer[6], "%lf", &scale0);
         calib_scale1[icalib * ndim + ival] = scale0;
         calib_eyepiece1[icalib * ndim + ival] = 1;
         ival++;
       } else if(!strncmp(buffer, "bin2 =", 6)) {
         sscanf(&buffer[6], "%lf", &scale0);
         calib_scale1[icalib * ndim + ival] = scale0;
         calib_eyepiece1[icalib * ndim + ival] = 2;
         ival++;
       } else if(!strncmp(buffer, "bin4 =", 6)) {
         sscanf(&buffer[6], "%lf", &scale0);
         calib_scale1[icalib * ndim + ival] = scale0;
         calib_eyepiece1[icalib * ndim + ival] = 4;
         ival++;
       } else if(!strncmp(buffer, "sign =", 6)) {
         sscanf(&buffer[6], "%d", &sign);
         sign1[icalib] = sign;
         jval++;
       } else if(!strncmp(buffer, "theta0 =", 8)) {
         sscanf(&buffer[8], "%lf", &theta0);
         theta01[icalib] = theta0;
         jval++;
       } else if(buffer[0] != '=') {
         fprintf(stderr, "read_calib_file/Bad syntax in line=>%s<\n",
                buffer);
         }
       } // if icalib >=0
     } // EOF if buffer...
   } // EOF if fgets...
 } // EOF while(!feof...)

*ncalib1 = icalib;
for(icalib = 0; icalib < *ncalib1; icalib++) {
   printf("read_calib_file_gili/calibration #%d\ncalib_date1=%.4f (%02d/%02d/%d)\n",
           icalib, calib_date1[icalib], calib_dd1[icalib], 
           calib_mm1[icalib], calib_year1[icalib]);
   for(i = 0; i < n_eyepieces1[icalib]; i++) {
     printf("calib_eyepiece[%d]=%d calib_scale1=%.4f\n", 
            i, calib_eyepiece1[icalib * ndim +i],
            calib_scale1[icalib * ndim + i]); 
     }
   printf("read_calib_file_gili/sign=%d\n", sign1[icalib]);
   printf("read_calib_file_gili/theta0=%.1f\n", theta01[icalib]);
 }

 status = 0;
 fclose(fp_calib);

return(status);
}

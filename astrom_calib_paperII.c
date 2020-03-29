/*************************************************************************
* Program latex_calib
* To convert a Latex array with measurements in pixels to array in arcseconds
* (and conversion of angles from XY axis to Noth-East referential)
*
* Format of input files:
\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}
& Name  & Epoch & Filter & Eyep. & $\rho$ & $\sigma_\rho$ & $\theta$ &
$\sigma_\theta$ & comments \\
 &       &       &        & (mm)
& (pixels) & (pixels) & (deg.) & (deg.) & \\
\hline % ------------------------------------------------------------------------$
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
* JLP
* Version 16/06/2005
*************************************************************************/
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

/* Maximum length for a line will be 170 characters: */
#define NMAX 360

#define NO_DATA -10000.
#define SQUARE(x) ((x)*(x))
#define MINI(x,y) ((x) < (y)) ? (x) : (y)
#define MAXI(x,y) ((x) < (y)) ? (y) : (x)
 
#define DEBUG
/*
*/

/* Maximum number of measurements per object: */
#define NMEAS 30
/* Maximum number of objects: */
#define NOBJ_MAX 500

/* Structure to define a measurement: */
typedef struct {
char comments[30][80]; 	/* Comments (lines starting with %) */
char data[180];         /* Data: line with filename, filter, measurements */
char filename[40];      /* Name of the FITS file used for this measurement */
char date[12];          /* Date, e.g., 29/12/2004 */
double epoch;            /* Julian day as a fraction of year, e;g., 2004.23 */
char filter[10];        /* Filter name: V, R, Sf */
int eyepiece;           /* Focal length of the magnifying eyepiece */
int quadrant;           /* Quadrant value */
int dquadrant;          /* Quadrant uncertainty (0 or 1) */
double rho;              /* Angular separation of the binary (arcseconds) */
double drho;             /* Error on rho (arcseconds) */
double theta;            /* Position angle of the companion (relative to North)*/
double dtheta; 		/* Error on theta (degrees) */
char notes[40];         /* Notes in the last column of the data line */
int flagged_out;        /* Flag to cancel output (for publication) */
} MEASURE; 

/* Structure to define an object */
typedef struct {
char wds[40];		/* WDS name */
char name[40];		/* Official binary name */
char ads[40];		/* ADS name */
MEASURE measure[NMEAS];	/* Measurements concerning this object */
char notes[40];		/* Notes which are common to all measurements */
double ra; 		/* Right ascension */
int dec; 		/* Declination */
int orbit;              /* Flag, set to one if orbit */
int nmeas;              /* Nber of measurements for this object*/
} OBJECT;

static int latex_calib_publi(FILE *fp_in, FILE *fp_out, double scale_10,
                             double scale_20, double theta0, double sign,
                             int i_filename, int i_date, int i_filter,
                             int i_eyepiece, int i_rho, int i_drho,
                             int i_theta, int i_dtheta, int i_notes,
                             int comments_wanted, int work_mode);
static int latex_calib_copy(FILE *fp_in, FILE *fp_out, double scale_10, 
                            double scale_20, double theta0, double sign, 
                            int i_date, int i_eyepiece, int i_rho, int i_drho, 
                            int i_theta, int i_dtheta, int comments_wanted);
static int calibrate_measures(OBJECT *obj, int nobj, double scale_10, 
                              double scale_20, double theta0, double sign);
static int calib_data_copy(char *b_data, char *b_out, double scale_10, 
                           double scale_20, double theta0, double sign, 
                           int i_date, int i_eyepiece, int i_rho, int i_drho, 
                           int i_theta, int i_dtheta);
static int decode_data(char *b_data, char *date, double *epoch, double *rho, 
                       double *drho, 
                       double *theta, double *dtheta, int *eyepiece, int i_date, 
                       int i_eyepiece, int i_rho, int i_drho, int i_theta, 
                       int i_dtheta);
static int read_dvalue(char *b_data, int *value, int icol);
static int read_fvalue(char *b_data, double *value, int icol); 
static int read_svalue(char *b_data, char *value, int icol); 
static int read_object_name(char *b_data, char *wds_name, char *official_name,
                            char *ads_name, int *orbit);
static int read_quadrant(char *notes, int *quadrant, int *dquadrant);
static int write_fvalue(char *b_data, char *b_out, double value, int icol,
                         int nber_of_decimals);
static int read_epoch_value(char *b_data, char *date, double *epoch , int icol); 
int julian(double aa, int mm, int idd, double time, double *djul);

static int add_new_measure(char *b_data, OBJECT *obj, int nobj, int i_filename, 
                           int i_date, int i_filter, int i_eyepiece, int i_rho,
                           int i_drho, int i_theta, int i_dtheta, int i_notes);
static int add_new_object(char *b_data, OBJECT *obj, int *nobj);
static int check_measure(char *b_data, int i_eyepiece, int i_rho, int i_drho, 
                         int i_theta, int i_dtheta);
static int read_measures(FILE *fp_in, int comments_wanted, OBJECT *obj, 
                         int *nobj, int i_filename, int i_date, 
                         int i_filter, int i_eyepiece, int i_rho, 
                         int i_drho, int i_theta, int i_dtheta, int i_notes);
static int write_publi_table(FILE *fp_out, int comments_wanted, OBJECT *obj, 
                             int *index_obj, int nobj);
static int compute_statistics(FILE *fp_out, OBJECT *obj, int nobj);
static int ra_sort_objects(OBJECT *obj, int *index_obj, int nobj);
static int JLP_QSORT_INDX(double *array, int *index, int *nn);
static void qs2(double *array, int *index, int left, int right);
static int mean_for_paper2(OBJECT *obj, int nobj);
static int is_record_file(char *filename);
static int is_direct_file(char *filename);
static int compute_mean_for_paper2(OBJECT *obj, int io, int jm);
static int good_quadrant(MEASURE *me);

int main(int argc, char *argv[])
{
double scale_10, scale_20, theta0, sign;
char filein[60], fileout[60];
int i_filename, i_date, i_filter, i_eyepiece, i_rho, i_drho, i_theta, i_dtheta; 
int i_notes;
int comments_wanted, publi_mode, work_mode;
FILE *fp_in, *fp_out;

if(argc == 7 && *argv[4]) argc = 5;
if(argc == 7 && *argv[3]) argc = 4;
if(argc == 7 && *argv[2]) argc = 3;
if(argc == 7 && *argv[1]) argc = 2;
if(argc != 5)
  {
  printf(" Syntax: latex_calib scale_10,scale_20,sign,theta0 in_file out_file comments_wanted,publi_mode \n");
  printf(" out_rho = in_rho * scale_10  (or in_rho * scale_20 if 20 mm eyepiece)\n");
  printf(" out_theta = in_theta * sign + theta0 \n");
  printf(" Example: runs latex_calib 0.03181,0.07438,1.,90.04 astrom04g.tex astrom04g_calib.tex 0,1  (without comments, publication mode)\n");
  printf(" Example: runs latex_calib 0.03181,0.07438,1.,90.04 astrom04g.tex astrom04g_calib.tex 1,0  (with comments, no publication)\n");
  exit(-1);
  }
else
  {
  sscanf(argv[1],"%lf,%lf,%lf,%lf", &scale_10, &scale_20, &sign, &theta0);
  strcpy(filein,argv[2]);
  strcpy(fileout,argv[3]);
  sscanf(argv[4],"%d,%d", &comments_wanted, &publi_mode);
  }
  printf(" OK: out_rho = in_rho * %g (for 10 mm eyepiece)\n", scale_10);
  printf(" OK: out_rho = in_rho * %g (for 20 mm eyepiece)\n", scale_20);
  printf(" OK: out_theta = in_theta * %g + %g \n", sign, theta0);

printf(" OK: filein=%s fileout=%s comments_wanted=%d publi_mode=%d\n",
       filein,fileout,comments_wanted,publi_mode);

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
fprintf(fp_out,"%%%% Automatic conversion with latex_calib\n");
fprintf(fp_out,"%%%% JLP / Version of 08/06/2005 \n");
fprintf(fp_out,"%%%% out_rho = in_rho * %g (for 10 mm eyepiece)\n", scale_10);
fprintf(fp_out,"%%%% out_rho = in_rho * %g (for 20 mm eyepiece)\n", scale_20);
fprintf(fp_out,"%%%% out_theta = in_theta * %g + %g \n", sign, theta0);
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

if(publi_mode) {
fprintf(fp_out,"%%%% Version 16/06/2005 for Paper II (Merate) \n");
fprintf(fp_out,"%%%% Publication mode: only direct measurement is kept when rho > 0.3\" \n");
fprintf(fp_out,"%%%% Weighted mean between direct and recorded measurements when rho <= 0.3\" \n");
fprintf(fp_out,"%%%% Minimum error for rho: 0.1 pixel  or 0.5%% \n");
fprintf(fp_out,"%%%% Minimum error for theta: 0.3 degree \n");
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
  work_mode = 1;
/* Version "publi", to generate a laTeX array formatted for publication:
* sort out the objects, compute mean values, perform calibration, etc.
*/ 
  if(publi_mode)
  latex_calib_publi(fp_in,fp_out, scale_10, scale_20, theta0, sign, 
                    i_filename, i_date, i_filter, i_eyepiece, i_rho, i_drho, 
                    i_theta, i_dtheta, i_notes, comments_wanted, work_mode);
/* Version "copy", to copy all the contents, and simply convert the 
* measurements to arcseconds and degrees relative to North 
*/ 
  else
  latex_calib_copy(fp_in,fp_out, scale_10, scale_20, theta0, sign, i_date, 
              i_eyepiece, i_rho, i_drho, i_theta, i_dtheta, comments_wanted);

fclose(fp_in);
fclose(fp_out);
return(0);
}
/*************************************************************************
*
* INPUT:
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
* comments_wanted: 1 if comments are wanted, 0 otherwise
*************************************************************************/
static int latex_calib_copy(FILE *fp_in, FILE *fp_out, double scale_10, 
                            double scale_20, double theta0, double sign, 
                            int i_date, int i_eyepiece, int i_rho, int i_drho, 
                            int i_theta, int i_dtheta, int comments_wanted)
{
char b_in[NMAX], b_data[NMAX];
int inside_array, line_is_opened, status;
int line_with_object_name;
char *pc, *pc1;

inside_array = 0;
line_is_opened = 0;

while(!feof(fp_in))
{
/* Maximum length for a line will be 170 characters: */
  if(fgets(b_in,170,fp_in))
  {
  line_with_object_name = 0;
  b_in[169] = '\0';
    if(!strncmp(b_in,"\\begin{tabular}",15)){
       inside_array = 1;
#ifdef DEBUG
printf(" OK: >%s< inside_array=%d\n", b_in, inside_array);
#endif
        }
    else if(!strncmp(b_in,"\\end{tabular}",13)){
       inside_array = 0;
#ifdef DEBUG
printf(" OK: >%s< inside_array=%d\n", b_in, inside_array);
#endif
       }
    else if(inside_array && b_in[0] != '%' && strncmp(b_in,"\\hline",6)) {
       if(!line_is_opened) {
         strcpy(b_data, b_in);
         }
/* Fill the data array with the next line */
       else {
/* Look for the first zero (end of string marker) in data buffer */
         b_data[119] = '\0';
         pc1 = b_data;
         while(*pc1) pc1++; 
         pc1--; 
/* Then copy the second line from there*/
         strcpy(pc1, b_in);
         }

/* Check if this line is ended with "\\": */
         line_is_opened = 1;
         pc = b_data;
         while(*pc) {
           if(!strncmp(pc,"\\\\",2)){
             line_is_opened = 0;
             pc += 2; *pc = '\n'; pc++; *pc = '\0';
             break;
             }
           pc++;
           } 
#ifdef DEBUG
printf(" Data line: >%s<\n", b_data);
printf(" line_is_opened=%d\n", line_is_opened);
#endif
     if(!line_is_opened) { 
       status = calib_data_copy(b_data, b_in, scale_10, scale_20, theta0, sign, 
                                i_date, i_eyepiece, i_rho, i_drho, i_theta, 
                                i_dtheta);
       }
    } 

    if(comments_wanted ) {
      if(!line_is_opened) fputs(b_in,fp_out);
      } else {
      if(!line_is_opened && b_in[0] != '%') fputs(b_in,fp_out);
      }
  } /* EOF if fgets() */
} /* EOF while loop */
return(0);
}
/*************************************************************************
* Publication mode
*
* INPUT:
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
* comments_wanted: 1 if comments are wanted, 0 otherwise
*************************************************************************/
static int latex_calib_publi(FILE *fp_in, FILE *fp_out, double scale_10, 
                             double scale_20, double theta0, double sign, 
                             int i_filename, int i_date, int i_filter,
                             int i_eyepiece, int i_rho, int i_drho,
                             int i_theta, int i_dtheta, int i_notes,
                             int comments_wanted, int work_mode)
{
OBJECT *obj;
int *index_obj;
int nobj;

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

read_measures(fp_in, comments_wanted, obj, &nobj, i_filename, i_date, 
              i_filter, i_eyepiece, i_rho, i_drho, i_theta, i_dtheta, i_notes);

#ifdef DEBUG
printf("Returned by read_measures: nobj = %d\n", nobj);
#endif

calibrate_measures(obj, nobj, scale_10, scale_20, theta0, sign);

/* Sort the objects according to their Right Ascension: */
ra_sort_objects(obj, index_obj, nobj);

/* Compute mean values for the objects with rho < 0.3" (Merate-Paper II) */
mean_for_paper2(obj, nobj);

write_publi_table(fp_out, comments_wanted, obj, index_obj, nobj);

compute_statistics(fp_out, obj, nobj);

free(index_obj);
free(obj);
return(0);
}
/*****************************************************************************
* Write the LateX array in a new format 
*
*****************************************************************************/
static int write_publi_table(FILE *fp_out, int comments_wanted, OBJECT *obj, 
                             int *index_obj, int nobj)
{
MEASURE *me;
int good_q, nm, io, nlines, jj, nl_max = 50;
char qflag[20], asterisc[20], exclam[20], q_error[1];
register int i, j;

/* Portrait: nl_max=50 */
/* Landscape: nl_max= 30 */
nl_max = 30;
strcpy(asterisc,"\\rlap{$^*$}");
strcpy(exclam,"\\rlap{!!}");

fprintf(fp_out,"\\def\\idem{''} \n");

nlines = -1;

for(i = 0; i < nobj; i++) {
  io = index_obj[i];
  nm = (obj[io]).nmeas;
/* New table header in publi_mode */
  if(nlines > nl_max) {
    fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
    fprintf(fp_out,"\\hline \n");
    fprintf(fp_out,"\\end{tabular} \n");
    fprintf(fp_out,"\\end{table} \n");
  }
  if((nlines == -1) || (nlines > nl_max)) {
    fprintf(fp_out,"\\begin{table} \n");
    fprintf(fp_out,"\\tabcolsep=1mm \n");
    fprintf(fp_out,"\\begin{tabular}{cccccccccccc} \n");
    fprintf(fp_out,"\\hline \n");
    fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
    fprintf(fp_out,"WDS & Name & ADS & Epoch & Fil. & Eyep. & $\\rho$ & $\\sigma_\\rho$ & $\\theta$ & $\\sigma_\\theta$ & Orb. & Notes \\\\ \n");
    fprintf(fp_out,"& &       &       &        & (mm)& (arcsec) & (arcsec) & (deg.) & (deg.) & & \\\\ \n");
    fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
    fprintf(fp_out,"\\hline \n");
    fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
    nlines = 0;
  }

  jj = 0;
/************* Loop on measurements: **********************/
  for(j = 0; j < nm; j++) {
    me = &(obj[io]).measure[j];
/* BOF case not_flagged */
    if(!me->flagged_out) {
/*
  printf("i= %d io=%d nm=%d j=%d flag=%d\n", i, io, nm, j, me->flagged_out);
  printf(" WDS=%s rho=%.2f drho=%.2f theta=%.2f dtheta=%.2f eyepiece=%d\n", 
           obj[io].wds, me->rho, me->drho, me->theta, me->dtheta, me->eyepiece);
*/
  good_q = good_quadrant(me);
  if(good_q == 1) strcpy(qflag,asterisc);
  else if(good_q == -1) {
     strcpy(qflag,exclam);
     printf(" write_publi_table/Quadrant incompatible with theta value for %s!\n",obj[io].wds);
     }
  else strcpy(qflag," ");

  *q_error = (me->dquadrant == 1) ? '?' : ' ';

  if(jj == 0) {
/*** Format for first line of an object */
  if(me->rho != NO_DATA) 
#ifdef Q_IN_NOTES
    fprintf(fp_out,"%s & %s & %s & %.3f & %s & %d & %.3f & %.3f & %.1f%s & %.1f & %d & %s q=%d%c\\\\\n", 
         obj[io].wds, obj[io].name, obj[io].ads, me->epoch, me->filter, 
         me->eyepiece, me->rho, me->drho, me->theta, qflag, me->dtheta, 
         obj[io].orbit, me->notes, me->quadrant, *q_error);
#else
    fprintf(fp_out,"%s & %s & %s & %.3f & %s & %d & %.3f & %.3f & %.1f%s & %.1f & %d & %s \\\\\n", 
         obj[io].wds, obj[io].name, obj[io].ads, me->epoch, me->filter, 
         me->eyepiece, me->rho, me->drho, me->theta, qflag, me->dtheta, 
         obj[io].orbit, me->notes);
#endif
  else
    fprintf(fp_out,"%s & %s & %s & %.3f & %s & %d & \\nodata & \\nodata & \\nodata & \\nodata & %d & %s \\\\\n", 
         obj[io].wds, obj[io].name, obj[io].ads, me->epoch, me->filter, 
         me->eyepiece, obj[io].orbit, me->notes);
  } /* EOF j == 0 */
/*** Format for subsequent lines of an object */
  else {
  if(me->rho != NO_DATA) 
#ifdef Q_IN_NOTES
    fprintf(fp_out,"\\idem & \\idem & \\idem & %.3f & %s & %d & %.3f & %.3f & %.1f%s & %.1f & %d & %s q=%d%c\\\\\n", 
         me->epoch, me->filter, 
         me->eyepiece, me->rho, me->drho, me->theta, qflag, me->dtheta, 
         obj[io].orbit, me->notes, me->quadrant, *q_error);
#else
    fprintf(fp_out,"\\idem & \\idem & \\idem & %.3f & %s & %d & %.3f & %.3f & %.1f%s & %.1f & %d & %s\\\\\n", 
         me->epoch, me->filter, 
         me->eyepiece, me->rho, me->drho, me->theta, qflag, me->dtheta, 
         obj[io].orbit, me->notes);
#endif
  else
    fprintf(fp_out,"\\idem & \\idem & \\idem & %.3f & %s & %d & \\nodata & \\nodata & \\nodata & \\nodata & %d & %s \\\\\n", 
         me->epoch, me->filter, me->eyepiece, obj[io].orbit, me->notes);
  } /* EOF j != 0 */
 nlines++;
 jj++;
 } /* EOF case not_flagged */
 } /* EOF loop on j */
} /* EOF loop on i */

 fprintf(fp_out,"& & & & & & & & & & & \\\\ \n");
 fprintf(fp_out,"\\hline \n");
 fprintf(fp_out,"\\end{tabular} \n");
 fprintf(fp_out,"\\end{table} \n");

return(0);
}
/*****************************************************************************
* Read the measurements and object parameters from the input file 
*
*****************************************************************************/
static int read_measures(FILE *fp_in, int comments_wanted, OBJECT *obj, 
                         int *nobj, int i_filename, int i_date, 
                         int i_filter, int i_eyepiece, int i_rho, 
                         int i_drho, int i_theta, int i_dtheta, int i_notes)
{
char b_in[NMAX], b_data[NMAX];
char wds_name[40], official_name[40], ads_name[40]; 
int inside_array, line_is_opened, status, orbit, line_to_reject;
char *pc, *pc1;

inside_array = 0;
line_to_reject = 0;
line_is_opened = 0;
wds_name[0] = '\0';
official_name[0] = '\0';
ads_name[0] = '\0';
orbit = 0;
*nobj = 0;

while(!feof(fp_in))
{
/* Maximum length for a line will be 170 characters: */
  if(fgets(b_in,170,fp_in))
  {
  b_in[169] = '\0';
    if(!strncmp(b_in,"\\begin{tabular}",15)){
       inside_array = 1;
       strcpy(b_in,"\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|} \n");
        }
    else if(!strncmp(b_in,"\\end{tabular}",13)){
       inside_array = 0;
       }
    else if(inside_array && b_in[0] != '%' && strncmp(b_in,"\\hline",6)) {
       if(!line_is_opened) {
         strcpy(b_data, b_in);
         }
/* Fill the data array with the next line */
       else {
/* Look for the first zero (end of string marker) in data buffer */
         b_data[119] = '\0';
         pc1 = b_data;
         while(*pc1) pc1++; 
         pc1--; 
/* Then copy the second line from there*/
         strcpy(pc1, b_in);
         }

/* Check if this line is ended with "\\": */
         line_is_opened = 1;
         pc = b_data;
         while(*pc) {
           if(!strncmp(pc,"\\\\",2)){
             line_is_opened = 0;
             pc += 2; *pc = '\n'; pc++; *pc = '\0';
             break;
             }
           pc++;
           } 
     if(!line_is_opened) { 
#ifdef DEBUG
printf(" Data line: >%s<\n", b_data);
#endif
       status = check_measure(b_data, i_eyepiece, i_rho, i_drho, 
                              i_theta, i_dtheta);
       
       if(!status) {
           add_new_measure(b_data, obj, *nobj, i_filename, i_date, i_filter, 
                           i_eyepiece, i_rho, i_drho, i_theta, i_dtheta, 
                           i_notes);
           } else {
/* Try to add a new object: */
              status = add_new_object(b_data, obj, nobj);
           }
       } /* EOF !line_is_opened */
    } 
  } /* EOF if fgets() */
} /* EOF while loop */
return(0);
}
/**************************************************************************
*
**************************************************************************/
static int add_new_object(char *b_data, OBJECT *obj, int *nobj)
{ 
OBJECT *ob;
char wds_name[40], official_name[40], ads_name[40], *pc, buffer[40];
int status, orbit, ih, im;

ob = &obj[*nobj];
status = read_object_name(b_data, wds_name, official_name, ads_name, &orbit);

if(!status) {
  strcpy(ob->wds, wds_name);
  strcpy(ob->ads, ads_name);
  strcpy(ob->name, official_name);
  ob->nmeas = 0;
  ob->orbit = orbit;
  (*nobj)++;

/* Decode Right Ascension from WDS name: */
  ob->ra = 0.;
  strcpy(buffer, wds_name);
  pc = buffer;
  while(*pc && (*pc != '+') && (*pc != '-')) pc++;
  if((*pc == '+') || (*pc == '-')) {
    *pc = '\0';
    sscanf(buffer,"%2d%3d", &ih, &im);
/* Divide "im" by 600, since "im" is measured in tenths of minutes: */
    ob->ra = (double)ih + ((double)im)/600.;
    }
/* Decode Declination from WDS name: */
  ob->dec = 0.;
  pc = wds_name;
  while(*pc && (*pc != '+') && (*pc != '-')) pc++;
  if((*pc == '+') || (*pc == '-')) {
    strcpy(buffer,pc);
    sscanf(buffer,"%d", &ih);
    ob->dec = ih;
    }
  
#ifdef DEBUG
  printf(" New object added: nobj = %d\n", *nobj);
  printf(" name= %s %s %s orbit=%d", ob->wds, ob->name, ob->ads, ob->orbit);
  printf(" ra= %f dec=%d\n", ob->ra, ob->dec);
#endif
  }

return(status);
}
/*************************************************************************
*
 22388+4419 = HO 295 AB & ADS 16138 & 2004. & & & & & & & orb \\
*************************************************************************/
static int read_object_name(char *b_data, char *wds_name, char *official_name,
                            char *ads_name, int *orbit)
{
char *pc, buff[60];
int istat, status, icol;
status = 1;

/* First step: decodes WDS name: */
icol = 1;
istat = read_svalue(b_data, buff, icol); 
if(istat) {
#ifdef DEBUG
  printf("read_object_name/Error reading column #%d\n", icol);
#endif
  status = -1;
  return(status);
  }

pc = buff;
buff[59]='\0';
/* Look for the second item (official name) starting after a "=" symbol: */
while(*pc && *pc != '=') pc++;
  if(*pc == '=') {
  pc++;
  strcpy(official_name, pc);
  status = 0;
  }

/* Look for the first item ending with a "=" symbol: */
if(!status) {
  pc = buff;
  while(*pc && *pc != '=') pc++;
    if(*pc == '=') {
    *pc = '\0';
    strcpy(wds_name, buff);
    status = 0;
    *orbit = 0;
    }

/* Reads ADS name in 2nd column: */
icol = 2;
ads_name[0] = '\0';
istat = read_svalue(b_data, buff, icol); 
if(!istat){
 pc = buff;
 while(*pc && strncmp(pc,"ADS",3)) pc++;
 if(!strncmp(pc,"ADS",3)){
    pc +=3;
    strcpy(ads_name, pc);
    }
 }

/* Then tries to read "orb" in column notes: */
icol = 10;
istat = read_svalue(b_data, buff, icol); 
if(!istat) { 
   pc = buff; 
   while(*pc && *pc != 'o') pc++;
   if(*pc == 'o') {
     if(!strncmp(pc,"orb",3)) *orbit = 1;
     }
   }
} /* EOF !status case */
return(status);
}
/*************************************************************************
* Check whether input line is a measurement
*
* INPUT:
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
*************************************************************************/
static int check_measure(char *b_data, int i_eyepiece, int i_rho, int i_drho, 
                         int i_theta, int i_dtheta)
{
int status;
double ww; 

status = read_fvalue(b_data, &ww, i_rho); 
if(!status) status = read_fvalue(b_data, &ww, i_drho); 
if(!status) status = read_fvalue(b_data, &ww, i_theta); 
if(!status) status = read_fvalue(b_data, &ww, i_dtheta); 
if(!status) status = read_fvalue(b_data, &ww, i_eyepiece); 

return(status);
}
/*************************************************************************
*
* INPUT:
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
* wds_name, official_name, ads_name : object designation in various catalogues  
* orbit: flag, set to one if orbit is known, 0 otherwise
*************************************************************************/
static int calibrate_measures(OBJECT *obj, int nobj, double scale_10, 
                              double scale_20, double theta0, double sign)
{
MEASURE *me;
double scale;
int nm;
register int i, j;

for(i = 0; i < nobj; i++) {
  nm = (obj[i]).nmeas;
  for(j = 0; j < nm; j++) {
    me = &(obj[i]).measure[j];

/* Scale according to eyepiece: */
    if(me->eyepiece == 10) scale = scale_10;
     else scale = scale_20;
 
   if(me->rho != NO_DATA) { 
    me->rho *= scale;
    me->drho *= scale;
    me->theta = me->theta * sign + theta0;
    while(me->theta < 0.) me->theta += 360.;
    while(me->theta >= 360.) me->theta -= 360.;
    }
   } /* EOF loop on j */
} /* EOF loop on i */

return(0);
}
/*************************************************************************
*
* INPUT:
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
* wds_name, official_name, ads_name : object designation in various catalogues  
* orbit: flag, set to one if orbit is known, 0 otherwise
*************************************************************************/
static int calib_data_copy(char *b_data, char *b_out, double scale_10, 
                      double scale_20, double theta0, double sign, int i_date, 
                      int i_eyepiece, int i_rho, int i_drho, int i_theta, 
                      int i_dtheta)
{
int status, eyepiece;
double epoch, rho, drho, theta, dtheta, scale;
char date[30];

 status = decode_data(b_data, date, &epoch, &rho, &drho, &theta, &dtheta, 
                      &eyepiece, i_date, i_eyepiece, i_rho, i_drho, i_theta, 
                      i_dtheta);
 if(status) {
  strcpy(b_out, b_data);
  return(1);
  }

 if(eyepiece == 10) scale = scale_10;
 else scale = scale_20;
 
if(rho != NO_DATA) { 
 rho *= scale;
 drho *= scale;
 theta = theta * sign + theta0;
 while(theta < 0.) theta += 360.;
 while(theta >= 360.) theta -= 360.;
 dtheta = dtheta;

/* rho and drho with 3 decimals */
status = write_fvalue(b_data, b_out, rho, i_rho, 3); 
if(!status) {
   strcpy(b_data, b_out);
   status = write_fvalue(b_data, b_out, drho, i_drho, 3); 
   }
/* theta and dtheta with 1 decimal */
if(!status) {
   strcpy(b_data, b_out);
   status = write_fvalue(b_data, b_out, theta, i_theta, 1); 
   }
if(!status) {
   strcpy(b_data, b_out);
   status = write_fvalue(b_data, b_out, dtheta, i_dtheta, 1); 
   }
} /* EOF !NO_DATA */

/* epoch with 3 decimals */
if((epoch > 0) && (!status)) {
   strcpy(b_data, b_out);
   status = write_fvalue(b_data, b_out, epoch, i_date, 3); 
   }
if(status) {
 printf("calib_data/Fatal error, updating array!\n");
 exit(-1);
 }

return(0);
}
/**************************************************************************
*
* INPUT:
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
*
* OUTPUT:
* rho, drho, theta, dtheta
* eyepiece
* epoch
*
**************************************************************************/
static int add_new_measure(char *b_data, OBJECT *obj, int nobj, int i_filename, 
                           int i_date, int i_filter, int i_eyepiece, int i_rho,
                           int i_drho, int i_theta, int i_dtheta, int i_notes)
{
MEASURE *me;
double epoch, rho, drho, theta, dtheta, ww;
char notes[40], filter[10], date[30], filename[40];
int eyepiece, quadrant, dquadrant, status, nm;

/* Return if no object has been entered yet: */
if(nobj <= 0) return(-1);

eyepiece = 0;
quadrant = 0;
dquadrant = 0;
rho = 0.;
drho = 0.;
theta = 0.;
dtheta = 0.;
epoch = 0.;
filter[0] = '\0';
notes[0] = '\0';

/* Read date: */
epoch = 0.;
status = read_epoch_value(b_data, date, &epoch, i_date); 

status = read_fvalue(b_data, &rho, i_rho); 
if(!status) status = read_fvalue(b_data, &drho, i_drho); 
if(!status) status = read_fvalue(b_data, &theta, i_theta); 
if(!status) status = read_fvalue(b_data, &dtheta, i_dtheta); 
if(!status) {
  status = read_fvalue(b_data, &ww, i_eyepiece); 
  eyepiece = (int)(ww+0.5);
  }

/* Read non-compulsory parameters */
   read_svalue(b_data, filename, i_filename);
   read_svalue(b_data, filter, i_filter);
   read_svalue(b_data, notes, i_notes);
/* Read the quadreant value if present, and removes "Q=" from notes */
   read_quadrant(notes, &quadrant, &dquadrant);

/* Store data and increase nber of measurements for this object */
if(!status) {
   nm = (obj[nobj-1]).nmeas;
   me = &(obj[nobj-1]).measure[nm];

   me->rho = rho; 
/* Minimum value for rho error: 0.1 pixel or 0.5% */
   drho = MAXI(drho, 0.1); 
   me->drho = MAXI(drho, rho*0.005); 
   me->theta = theta; 
/* Minimum value for theta error: 0.3 degree */
   me->dtheta = MAXI(dtheta, 0.3); 
   me->eyepiece = eyepiece; 
   me->quadrant = quadrant; 
   me->dquadrant = dquadrant; 
   me->epoch = epoch; 
   me->flagged_out = 0;
   strcpy(me->filename, filename);
   strcpy(me->filter, filter);
   strcpy(me->notes, notes);
   strcpy(me->date, date);
   ((obj[nobj-1]).nmeas)++;

#ifdef DEBUG
     printf(" nobj=%d, new measure successfully added (nm=%d)\n", nobj,
              (obj[nobj-1]).nmeas);
     printf(" rho=%.2f drho=%.2f theta=%.2f dtheta=%.2f eyep.=%d Q=%d dQ=%d notes=>%s<\n", 
                   me->rho, me->drho, me->theta, me->dtheta, me->eyepiece,
                   me->quadrant, me->dquadrant, me->notes);
#endif
   }
#ifdef DEBUG
else
  printf("add_new_measure/Failure adding new measurement \n");
#endif

return(status);
}
/**************************************************************************
*
* INPUT:
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
*
* OUTPUT:
* rho, drho, theta, dtheta
* eyepiece
* epoch
*
**************************************************************************/
static int decode_data(char *b_data, char *date, double *epoch, double *rho, 
                       double *drho, double *theta, double *dtheta, int *eyepiece,
                       int i_date, int i_eyepiece, int i_rho, int i_drho, 
                       int i_theta, int i_dtheta)
{
double ww;
int status;

*eyepiece = 0;
*rho = *drho = *theta = *dtheta = 0;

/* Read date: */
*epoch = 0.;
status = read_epoch_value(b_data, date, epoch, i_date); 

status = read_fvalue(b_data, rho, i_rho); 
if(!status) status = read_fvalue(b_data, drho, i_drho); 
if(!status) status = read_fvalue(b_data, theta, i_theta); 
if(!status) status = read_fvalue(b_data, dtheta, i_dtheta); 
if(!status) {
  status = read_fvalue(b_data, &ww, i_eyepiece); 
  *eyepiece = (int)(ww+0.5);
  }

#ifdef DEBUG
if(!status) printf(" rho=%.2f drho=%.2f theta=%.2f dtheta=%.2f eyepiece=%d\n", 
                   *rho, *drho, *theta, *dtheta, *eyepiece);
#endif

return(status);
}
/**************************************************************************
* Read decimal value in column #icol from b_data string
*
**************************************************************************/
static int read_dvalue(char *b_data, int *value, int icol) 
{
int ival, status;
char buff[40];

*value = 0.;
status = read_svalue(b_data, buff, icol); 
if(!status) { 
   ival = sscanf(buff, "%d", value);
printf("read_dvalue/buff=>%s< value=%d ival=%d\n", buff, *value, ival);
/*
*/
   if(ival <= 0) status = 1;
  }

return(status);
}
/**************************************************************************
* Read double value in column #icol from b_data string
*
**************************************************************************/
static int read_fvalue(char *b_data, double *value, int icol) 
{
int ival, status;
char buff[40], nodata[40];

*value = 0.;
status = read_svalue(b_data, buff, icol); 
if(!status) { 
   sscanf(buff, "%s", nodata);
   if(!strncmp(nodata,"\\nodata",7)) { 
/*
     printf("read_fvalue: nodata found! \n");
*/
     *value = NO_DATA;
     }
   else {
   ival = sscanf(buff, "%lf", value);
   if(ival <= 0) { 
/*
      printf("read_fvalue/buff=>%s< value=%.2f ival=%d\n", buff, *value, ival);
*/
      status = 1;
      }
   }
  }
return(status);
}
/**************************************************************************
* Read epoch value from date information in column #icol from b_data string
*
* Input format: dd/mm/yy, e.g. 12/2/2004 or 01/12/1998 or 31/06/2002
* Output: epoch, as a fraction of year, e.g. 2004.234
*
**************************************************************************/
static int read_epoch_value(char *b_data, char *date, double *epoch , int icol) 
{
int ival, status, dd, mm, iyy;
double yy, time, jdate, jdate_jan01;

/* Read date: */
date[0] = '\0';
status = read_svalue(b_data, date, icol); 

/* Removes the blanks in "date" string: */
sscanf(date, "%s", date);

if(!status) { 
   ival = sscanf(date, "%d/%d/%d", &dd, &mm, &iyy);
/*
printf("read_epoch_value/date=>%s< dd=%d mm=%d iyy=%d ival=%d\n", 
        date, dd, mm, iyy, ival);
*/
   if(ival < 3) status = 1;
  }

if(!status) { 
/* Assume observations at 10:30 pm, local time, i.e., 8:30 (U.T) in summer */
/* Assume observations at 9:30 pm, local time, i.e., 8:30 (U.T) in winter */
time = 20.5;
yy = (double)iyy;
status = julian(yy, mm, dd, time, &jdate);

/* 1st of January at 0.00: */
julian(yy, 1, 1, 0., &jdate_jan01);

*epoch = yy + (jdate - jdate_jan01)/365.25;

/*
printf("read_epoch_value/jdate=%f jdate_jan01=%f epoch=%f\n", 
        jdate, jdate_jan01, *epoch); 
*/
}

return(status);
}
/**************************************************************************
* Read string value in column #icol from b_data string
*
**************************************************************************/
static int read_svalue(char *b_data, char *value, int icol) 
{
int ic, status, column_is_found;
char buff[NMAX], data[NMAX], *pc;

strcpy(data, b_data);

pc = data;
data[NMAX-1] = '\0';
column_is_found = 0;
ic = 1;
buff[0] = '\0';
while(*pc) {
  if(ic == icol) {
    column_is_found = 1; 
    strcpy(buff,pc);
    break;
    }
  if(*pc == '&') { 
    ic++;
    }
  pc++;
  }
*pc = '\0';
/* Return if column not found, or empty */
if(!buff[0]) return(-1); 

/* Otherwise go on analysis: */
status = 1;
buff[NMAX-1] = '\0';
pc = buff;
while(*pc) {
  if(*pc == '&' || !strncmp(pc,"\\\\",2)) {
    *pc = '\0'; 
    *value = '\0';
    strcpy(value,buff);
    if(*value) status = 0;
    break;
    }
  pc++;
  }

return(status);
}
/**************************************************************************
* Write a double value in column #icol to b_out string
*
**************************************************************************/
static int write_fvalue(char *b_data, char *b_out, double value, int icol, 
                         int nber_of_decimals) 
{
int ic, column_is_found, istart, iend;
char data[360], *pc;
register int i;

strcpy(data, b_data);

pc = data;
column_is_found = 0;
ic = 1;
i = 0;
while(*pc) {
  if(ic == icol && !column_is_found) {
    column_is_found = 1; 
    istart = i;
    }
  else if(ic == icol+1) {
    iend = i-1;
    break;
    }
  if(*pc == '&') { 
    ic++;
    }
  pc++;
  i++;
  }
/* Return if column not found, or empty */
if(istart == 0 || iend == istart) return(-1); 

strcpy(b_out, b_data);
switch(nber_of_decimals) {
  case 1:
    sprintf(&b_out[istart],"%8.1f ",value);
    break;
  case 2:
    sprintf(&b_out[istart],"%8.2f ",value);
    break;
  case 3:
  default:
    sprintf(&b_out[istart],"%8.3f ",value);
    break;
  case 4:
    sprintf(&b_out[istart],"%9.4f ",value);
    break;
  }
strcpy(&b_out[istart+9],&b_data[iend]);

/*
printf("update_value/from >%s< to >%s< (value=%.2f)\n", b_data, b_out, value);
*/

return(0);
}
/*********************************************************************
* Subroutine JULIAN to compute the Julian day of an observation:
* (from "cel_meca.c")
*
* The Julian day begins at Greenwich mean noon (at 12 U.T.)
*
* Here also the Gregorian calendar reform is taken into account.
* Thus the day following 1582 October 4 is 1582 October 15.
*
* The B.C. years are counted astronomically. Thus the year
* before the year +1 is the year 0.
*
* Input:
* AA, MM, IDD, TIME : year,month, day, time of the observation
* DJUL : Julian day
**********************************************************************/
int julian(double aa, int mm, int idd, double time, double *djul)
{
double day1, year1, date_obs, date_reform;
long month1, ia1, ib1;

  day1 = time/24. + (double)idd;
/* First the year after the 1st March ... */
  if(mm > 2) {
     year1 = aa;
     month1 = mm;
    }
   else {
     year1 = aa - 1;
     month1 = mm + 12;
    }

/* Then check if after the Gregorian reform: */
    date_obs = aa + ((int)(275 * mm / 9)
               - 2. * (int) ((mm + 9) / 12) + idd - 30 ) / 365.;
    date_reform = 1582. + 289./365.;
    if(date_obs >= date_reform) {
         ia1 = (int) (year1 / 100.);
         ib1 = 2 - ia1 + (int) (((double)ia1)/4.);
       }
    else
         ib1 = 0;

/* Now final formula: */
      *djul = (int)(365.25 * year1) + (int)(30.6001 * (month1 + 1))
              + day1 + 1720994.5 + ib1;

return(0);
}
/***************************************************************************
* ra_sort_objects
* Calling routine  int JLP_QSORT_INDX(double *array, int *index, int *nn)
* INPUT:
*  array[nn]: array to be sorted
*
* OUTPUT:
*  array[nn]: sorted array
*  index[nn]: array giving the index of the input array,
*             to sort other arrays in the same way if necessary
*             (f.i. array2[i] := array2[index[i]])
****************************************************************************/
static int ra_sort_objects(OBJECT *obj, int *index_obj, int nobj)
{
double *ra;
int j1, j2, nswap;
register int i;

ra = (double *)malloc((nobj) * sizeof(double));

for(i = 0; i < nobj; i++) ra[i] = obj[i].ra;

JLP_QSORT_INDX(ra, index_obj, &nobj);

/* Final check on declination criterion too: */
nswap = 1;
while(nswap > 0) {
nswap = 0;
for(i = 0; i < nobj-1; i++) {
 j1 = index_obj[i];
 j2 = index_obj[i+1];
/* Swap indices if declination is not sorted out: */
 if((obj[j1].ra == obj[j2].ra) && (obj[j1].dec > obj[j2].dec)) {
      index_obj[i] = j2;
      index_obj[i+1] = j1;
      nswap++;
#ifdef DEBUG
      printf("ra1=%f ra2=%f dec1=%d dec2=%d\n", obj[j1].ra, obj[j2].ra, 
              obj[j1].dec, obj[j2].dec);
#endif
   } 
 } /* EOF loop on i */
#ifdef DEBUG
printf("Sorting declination now: nswap=%d for this iteration\n",nswap);
#endif
} /* EOF while loop */

free(ra);
return(0);
}
/****************************************************************************
* Quicksort (Cf. "C: The complete reference", Herbert Schildt)
*
* INPUT:
*  array[nn]: array to be sorted
*
* OUTPUT:
*  array[nn]: sorted array
*  index[nn]: array giving the index of the input array, 
*             to sort other arrays in the same way if necessary
*             (f.i. array2[i] := array2[index[i]])
****************************************************************************/
static int JLP_QSORT_INDX(double *array, int *index, int *nn)
{
register int i;

/* Initialization of index array: */
for(i = 0; i < *nn; i++) index[i] = i;

if(*nn < 2) return(0);
 qs2(array, index, 0, (*nn)-1);

return(0);
}
/**************************************************************************
* The Quicksort, with index array 
***************************************************************************/
static void qs2(double *array, int *index, int left, int right)
{
register int i, j;
int iy;
double x, y;

i = left; j = right;
/* Take the element in the middle as reference to partition the array: */
x = array[(left+right)/2];

/* Put the elements < x to the left and those > x to the right: */ 
do {
  while(array[i] < x && i < right) i++;
  while(x < array[j] && j > left) j--;
  if(i <= j) {
/* Exchange array[i] and array[j]: */
    y = array[i];
    array[i] = array[j];
    array[j] = y;
/* Exchange index[i] and index[j]: */
    iy = index[i];
    index[i] = index[j];
    index[j] = iy;
    i++; j--;
  }
} while(i<=j);

if(left < j) qs2(array, index, left, j);
if(i < right) qs2(array, index, i, right);

return;
}
/*************************************************************************
* Compute mean values for the objects with rho < 0.3" (Merate-Paper II) 
**************************************************************************/
static int mean_for_paper2(OBJECT *obj, int nobj)
{
MEASURE *me, *me_prev;
int rec_dir, rec_dir_prev, nm, no_data_prev, no_data;
register int i, j;


for(i = 0; i < nobj; i++) {
  nm = (obj[i]).nmeas;
  me_prev = NULL;
  for(j = 0; j < nm; j++) {
    me = &(obj[i]).measure[j];
/* rec_dir = 1  if recorded file
*          = -1 if direct file
*          = 0 otherwise
*/
    if(is_record_file(me->filename) == 1)
      rec_dir = 1;
    else if(is_direct_file(me->filename) == 1)
      rec_dir = -1;
    else
      rec_dir = 0;
    no_data = ((me->rho == NO_DATA) || (me->theta == NO_DATA)) ? 1 : 0; 
/* Store current parameters (they become "previous" parameters in next 
* iteration if same object) */
    if(me_prev == NULL) {
      rec_dir_prev = rec_dir;
      no_data_prev = no_data;
      me_prev = me;
      } else {
/* Process measurements when two measurements refer to same epoch,
* with direct/recorded files: */ 
      if((me->epoch == me_prev->epoch) 
         && (me->filter[0] == me_prev->filter[0])
         && (me->eyepiece == me_prev->eyepiece)
         && (rec_dir * rec_dir_prev) == -1) {
/* If two measurements have been made: */
        if(!no_data  && !no_data_prev) {
/* Handle quadrant: */
          if(me->quadrant == -1 && me_prev->quadrant == 0) me_prev->quadrant = -1; 
          if(me->quadrant == 0 && me_prev->quadrant == -1) me->quadrant = -1; 
          if(me->quadrant > 0 && me_prev->quadrant <= 0) {
              me_prev->quadrant = me->quadrant; 
              me_prev->dquadrant = me->dquadrant; 
              }
          if(me_prev->quadrant > 0 && me->quadrant <= 0) {
              me->quadrant = me_prev->quadrant; 
              me->dquadrant = me_prev->dquadrant; 
              }
/* Handle rho and theta: */
          if(me->rho > 0.3) {
/* if rec and dir_prev */
            if(rec_dir == 1 && rec_dir_prev == -1) me->flagged_out = 1;
/* if dir and rec_prev */
            if(rec_dir == -1 && rec_dir_prev == 1) me_prev->flagged_out = 1;
            }
          else compute_mean_for_paper2(obj, i, j);
        } /* EOF !no_data && !no_data_prev */
/* Special handling when NO_DATA were present for both measurements: */
/* Removes the second measurement (NB: the 2nd is the same as the first): */
        else if(no_data  && no_data_prev) me->flagged_out = 1;
/* Case when only one mesurement over two files: */
        else {
         if(no_data) me->flagged_out = 1;
         else if(no_data_prev) me_prev->flagged_out = 1;
         }
/* Second file was linked to first, I neutralize it for further inspection */
      me_prev = NULL;
       } /* EOF epoch == epoch_prev */
/* Second file was not linked to first, it may be linked to next: */
  else {
       rec_dir_prev = rec_dir;
       no_data_prev = no_data;
       me_prev = me;
       }
     } /* EOF me_prev != NULL */

/* Go to next measurement */
  } /* EOF loop on j (measurements) */
/* Go to next object */
} /* EOF loop on i  (objects) */

return(0);
}
/************************************************************************
* Check if recorded file,
* i.e. filename with _Rr _Vr or _Sfr
*
************************************************************************/
static int is_record_file(char *filename)
{
char *pc;
int is_record_file;

is_record_file = 0;

/* Check if recorded file: */
pc = filename;
while(*pc && strncmp(pc,"_Vr",3) && strncmp(pc,"_Rr",3)
      && strncmp(pc,"_Sfr",4) && strncmp(pc,"_s.f.r",6)) pc++;
if(!strncmp(pc,"_Vr",3) || !strncmp(pc,"_Rr",3)
      || !strncmp(pc,"_Sfr",4) || !strncmp(pc,"_s.f.r",6)) is_record_file = 1; 

return(is_record_file);
}
/************************************************************************
* Check if direct file,
* i.e. filename with _Rd _Vd or _Sfd
*
************************************************************************/
static int is_direct_file(char *filename)
{
char *pc;
int is_direct_file;

is_direct_file = 0;

/* Check if direct file: */
pc = filename;
while(*pc && strncmp(pc,"_Vd",3) && strncmp(pc,"_Rd",3)
      && strncmp(pc,"_Sfd",4) && strncmp(pc,"_s.f.d",6)) pc++;
if(!strncmp(pc,"_Vd",3) || !strncmp(pc,"_Rd",3)
      || !strncmp(pc,"_Sfd",4) || !strncmp(pc,"_s.f.d",6)) is_direct_file = 1; 

return(is_direct_file);
}
/************************************************************************
*
* io: object index in obj
* jm: measurement index in obj[io].measure
************************************************************************/
static int compute_mean_for_paper2(OBJECT *obj, int io, int jm)
{
MEASURE *me1, *me2;
double w1, w2, sum;

me1 = &(obj[io]).measure[jm-1];
me2 = &(obj[io]).measure[jm];
w1 = (me2->drho / (me1->drho + me2->drho))
     + (me2->dtheta / (me1->dtheta + me2->dtheta));
w2 = (me1->drho / (me1->drho + me2->drho))
     + (me1->dtheta / (me1->dtheta + me2->dtheta));
sum = w1 + w2;
  if(sum == 0) {
     printf("compute_mean_forp_paper2/Fatal: rms error is null for io=%d jm=%d\n", io, jm);
     exit(-1);
     }
w1 /= sum; 
w2 /= sum; 
/*
printf(" WDS=%s ADS=%s NAME=%s io=%d jm=%d dt1=%f dt2=%f dr1=%f dr2=%f w1=%f w2 =%f \n", 
       (obj[io]).wds, (obj[io]).ads, (obj[io]).name, io, jm, me1->dtheta, 
       me2->dtheta, me1->drho, me2->drho, w1, w2);
*/
me1->rho = me1->rho * w1 + me2->rho * w2;
me1->theta = me1->theta * w1 + me2->theta * w2;
me1->drho = sqrt((SQUARE(me1->drho) * w1 + SQUARE(me2->drho) * w2));
me1->dtheta = sqrt((SQUARE(me1->dtheta) * w1 + SQUARE(me2->dtheta) * w2));

me2->flagged_out = 1;
return(0);
}
/***************************************************************************
* Read the quadrant value if present, and removes "Q=*" from notes 
*
* INPUT:
*  notes: string from notes column of input file 
*
* OUTPUT
* quadrant: 0 if Q was not present
*           -1 if Q=?
*           1, 2, 3 or 4 if Q=1, Q=2, Q=3 or Q=4
*  notes: same as input string, but without "Q=*" 
*
***************************************************************************/
static int read_quadrant(char *notes, int *quadrant, int *dquadrant)
{
char *pc, buffer[40];
int status, k;
status = 0;

*quadrant=0;
*dquadrant=0;

pc = notes;
k = 0;
while(*pc && strncmp(pc,"Q=",2)) {pc++; buffer[k++] = *pc;}

if(!strncmp(pc,"Q=",2)){
 k--;
 pc +=2;
 switch(*pc) {
   case '1':
    *quadrant=1;
    break;
   case '2':
    *quadrant=2;
    break;
   case '3':
    *quadrant=3;
    break;
   case '4':
    *quadrant=4;
    break;
   case '?':
   default:
    *quadrant=-1;
    break;
   }
/* Case when Q=1?, Q=2?, Q=3? or Q=4? */
 if(*pc) pc++;
 if(*pc == '?') {pc++; *dquadrant = 1;}

/* Normal syntax is f.i. "Q=4," */
 while(*pc) {pc++; buffer[k++] = *pc;}
 buffer[k++] = '\0';
 }

/* Update notes, after removing quadrant indication: */
strcpy(notes,buffer);

return(status);
}
/*********************************************************************************
* Check if theta value is compatible with the quadrant value
*********************************************************************************/
static int good_quadrant(MEASURE *me)
{
double theta;
int quad, good_quad;

good_quad = 0;
quad = me->quadrant;
theta = me->theta;

switch(quad)
  {
   case 1:
     if(theta > -1. && theta < 91.) good_quad = 1;
     else good_quad = -1;
     break;
   case 2:
     if(theta > 89. && theta < 181.) good_quad = 1;
     else good_quad = -1;
     break;
   case 3:
     if(theta > 179. && theta < 271.) good_quad = 1;
     else good_quad = -1;
     break;
   case 4:
     if(theta > 269. && theta <= 360.) good_quad = 1;
     else good_quad = -1;
     break;
   default:
     break;
  }

return(good_quad);
}
/*********************************************************************
*
*********************************************************************/
static int compute_statistics(FILE *fp_out, OBJECT *obj, int nobj)
{
MEASURE *me;
int nm, nmeas, no_detected, no_data;
int nquad, nquad_uncert, is_quad, is_bad_quad, nbad_quad;
register int i, j;

/* Compute number of measurements: */
nmeas = 0;
no_detected = 0;
nquad = 0;
nquad_uncert = 0;
nbad_quad = 0;

for(i = 0; i < nobj; i++) {
 nm = (obj[i]).nmeas;
   for(j = 0; j < nm; j++) {
   me = &(obj[i]).measure[j];
     if(!me->flagged_out) {
     nmeas ++;
     no_data = ((me->rho == NO_DATA) || (me->theta == NO_DATA)) ? 1 : 0; 
     no_detected += no_data;
     is_quad = (me->quadrant > 0) ? 1 : 0;
     nquad += is_quad;
     nquad_uncert += me->dquadrant;
     is_bad_quad = (good_quadrant(me) == -1) ? 1 : 0;
     nbad_quad += is_bad_quad;
     }
   }
 }

printf(" Number of objects: %d \n", nobj);
fprintf(fp_out, " Number of objects: %d \n \n", nobj);

printf(" Number of observations: %d with %d measurements and %d cases of no detection \n ", 
        nmeas, nmeas - no_detected, no_detected);
fprintf(fp_out, " Number of observations: %d with %d measurements and %d cases of no detection \n \n ", 
        nmeas, nmeas - no_detected, no_detected);

printf("Quadrant was determined for %d measurements (including %d uncertain determinations)\n", 
        nquad, nquad_uncert);
fprintf(fp_out, "Quadrant was determined for %d measurements (including %d uncertain determinations)\n \n", 
        nquad, nquad_uncert);

printf("Warning: %d quadrant values are inconsistent with theta values in this table! \n",
         nbad_quad);
fprintf(fp_out,"Warning: %d quadrant values are inconsistent with theta values in this table! \n \n",
         nbad_quad);
return(0);
}

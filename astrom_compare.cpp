/*************************************************************************
* astrom_compare
* Program to compare the measurements of two astrom files
*
* Format of input files:
& 120109\_ads684ab\_Rd\_8\_a & 12/01/2009 & R & 20 & 11.71 & 0.10 & -19.37 & 0.9 & Q=? EP=2009.0852 \\
*
Possible keywords in the notes:
EP=2004.1323 (epoch)
Q=2         (Quadrant with restricted triple-correlation)
LQ=3        (Quadrant with long integration)
*
* JLP
* Version 31/03/2020
*************************************************************************/
#include "astrom_utils1.h" 
#include "astrom_utils2.h" 

static int astrom_compare_meas(char *filein1, char *filein2, char *fileout);
static int astrom_compare_files(FILE *fp_in1, char *filein2, FILE *fp_out, 
                       int i_filename, int i_date, int i_filter,
                       int i_eyepiece, int i_rho, int i_drho,
                       int i_theta, int i_dtheta, int i_notes,
                       int comments_wanted, char *filein);
static int astrom_read_measures_cmp(FILE *fp_in, char *filein2,
                         int comments_wanted, 
                         OBJECT *obj, int *nobj, int i_filename, int i_date,
                         int i_filter, int i_eyepiece, int i_rho,
                         int i_drho, int i_theta, int i_dtheta, int i_notes,
                         int with_wds_data);
static int astrom_get_meas_from_file(char *filein2, char *obj_name1,
                       double epoch1, 
                       int i_filename, int i_date, int i_filter,
                       int i_eyepiece, int i_rho, int i_drho,
                       int i_theta, int i_dtheta, int i_notes,
                       double *rho2, double *err_rho2, double *theta2,
                       double *err_theta2);

#define DEBUG

/*********************************************************************/
int main(int argc, char *argv[])
{
char filein1[64], filein2[64], fileout[64];

/* If command line with "runs" */
if(argc == 7){
 if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
 else argc = 1;
 }


if(argc != 4)
  {
  printf(" Syntax: astrom_compare in_astrom_file1 in_astrom_file2 out_file\n");
  printf(" Example: runs astrom_compare astrom05a.tex astrom05b.tex cmp_5a_5b.txt \n");
  exit(-1);
  }
else
  {
  strcpy(filein1,argv[1]);
  strcpy(filein2,argv[2]);
  strcpy(fileout,argv[3]);
  }

printf(" OK: filein1=%s filein2=%s fileout=%s \n", filein1, filein2, fileout);

/* Scan the file and add epoch from FITS autocorrelation files:
* (in "astrom_utils2.c")
*/ 
  astrom_compare_meas(filein1, filein2, fileout); 

return(0);
}
/*************************************************************************
*
*
*************************************************************************/
static int astrom_compare_meas(char *filein1, char *filein2, char *fileout)
{
FILE *fp_in1, *fp_out;
int i_filename, i_date, i_filter, i_eyepiece, i_rho, i_drho, i_theta, i_dtheta;
int i_notes, comments_wanted, publi_mode, input_with_header;


if((fp_in1 = fopen(filein1,"r")) == NULL) {
  fprintf(stderr, " Fatal error opening input file1 %s \n", filein1);
  return(-1);
  }

if((fp_out = fopen(fileout,"w")) == NULL)
{
printf(" Fatal error opening output file %s \n",fileout);
return(-1);
}
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(fp_out,"%%%% astrom_compare file1=%s and file2=%s\n", filein1, filein2);
fprintf(fp_out,"%%%% JLP / Version of 31/03/2020 \n");
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

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
  astrom_compare_files(fp_in1, filein2, fp_out, 
                     i_filename, i_date, i_filter, i_eyepiece, i_rho, i_drho,
                     i_theta, i_dtheta, i_notes, comments_wanted, filein1);
fclose(fp_in1);
fclose(fp_out);
return(0);
}
/*************************************************************************
*
* INPUT:
* i_filename: column nber of the filename used for this measurement
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
* comments_wanted: 1 if comments are wanted, 0 otherwise
*
*************************************************************************/
static int astrom_compare_files(FILE *fp_in1, char *filein2, FILE *fp_out, 
                       int i_filename, int i_date, int i_filter,
                       int i_eyepiece, int i_rho, int i_drho,
                       int i_theta, int i_dtheta, int i_notes,
                       int comments_wanted, char *filein1)
{
OBJECT *obj1, *obj2;
int *index_obj1, tabular_only;
int i, nobj1 = 0, nobj2 = 0, with_wds_data;

if((obj1 = (OBJECT *)malloc(NOBJ_MAX * sizeof(OBJECT))) == NULL) {
  printf("astrom_calib_publi/Fatal error allocating memory space for OBJECT: nobj_max=%d\n",
          NOBJ_MAX);
  exit(-1);
  }
// Initialize the number of measurements to zero
for(i = 0; i < NOBJ_MAX; i++) (obj1[i]).nmeas = 0;

if((index_obj1 = (int *)malloc((NOBJ_MAX) * sizeof(int))) == NULL) {
  printf("astrom_calib_publi/Fatal error allocating memory space for index_obj: nobj_max=%d\n",
          NOBJ_MAX);
  exit(-1);
  }
with_wds_data = 0;
astrom_read_measures_cmp(fp_in1, filein2, comments_wanted, obj1, &nobj1, 
                         i_filename, i_date, i_filter, i_eyepiece, i_rho, 
                         i_drho, i_theta, i_dtheta, i_notes, with_wds_data);

#ifdef DEBUG
printf("Returned by astrom_read_measures: nobj1 = %d nmeas=%d rho=%.2f theta=%.2f\n", 
        nobj1, obj1[nobj1-1].nmeas, obj1[nobj1-1].measure[0].rho, 
        obj1[nobj1-1].measure[0].theta);
#endif

/* Sort the objects according to their name: */
astrom_name_sort_objects(obj1, index_obj1, nobj1);

/* For big tables, should set this parameter to 1: */
tabular_only = 0;
astrom_write_publi_table_gili(fp_out, comments_wanted, obj1, index_obj1, nobj1,
                         tabular_only);

// Number of objects, quadrants, etc...
astrom_compute_statistics(fp_out, obj1, nobj1, filein1);

free(index_obj1);
free(obj1);
return(0);
}
/*****************************************************************************
* Read the measurements and object parameters from the input file
*
* INPUT:
* i_filename: column nber of the filename used for this measurement
* i_eyepiece: column nber of eyepiece focal length information
* i_rho: column nber with rho values
* i_drho: column nber with drho values
* i_theta: column nber with theta values
* i_dtheta: column nber with dtheta values
* i_notes: column nber with the notes
* nobj: number of objects already entered into *obj
* with_wds_data: flag set to one if wds data are in the input file
*
* OUTPUT:
* obj: OBJECT structure
* nobj: total number of objects entered into *obj
*****************************************************************************/
static int astrom_read_measures_cmp(FILE *fp_in, char *filein2,
                         int comments_wanted, 
                         OBJECT *obj, int *nobj, int i_filename, int i_date,
                         int i_filter, int i_eyepiece, int i_rho,
                         int i_drho, int i_theta, int i_dtheta, int i_notes,
                         int with_wds_data)
{
char b_in[NMAX], b_data[NMAX];
char wds_name[40], discov_name[40], ads_name[40];
int inside_array, line_is_opened, status, orbit, line_to_reject, i_obj;
int input_with_header;
char *pc, *pc1;
char obj_name1[64];
double epoch1, rho2, err_rho2, theta2, err_theta2;

#ifdef DEBUG
printf("astrom_read_measures: input/nobj=%d\n", *nobj);
#endif

/* JLP 2014: automatic check with the first line */
/* Read first line: */
  fgets(b_in,170,fp_in);
/* input_with_header: 1 if input Latex file has a full header
*                    0 if no header at all
*/
 if(!strncmp(b_in,"% FILE_WITH_HEADER", 18)) input_with_header = 1;
 else input_with_header = 0;

/* Assume we are inside the array if (input_with_header == 0): */
if(input_with_header == 0)
  inside_array = 1;
 else
  inside_array = 0;

line_to_reject = 0;
line_is_opened = 0;
wds_name[0] = '\0';
discov_name[0] = '\0';
ads_name[0] = '\0';
orbit = 0;
while(!feof(fp_in))
{
/* Maximum length for a line will be 170 characters: */
  if(fgets(b_in,170,fp_in))
  {
  b_in[169] = '\0';
/* NEW/2009: I remove ^M (Carriage Return) if present: */
  pc = b_in;
  while(*pc) {
  if(*pc == '\r') *pc = ' ';
  pc++;
  }

    if(!strncmp(b_in,"\\begin{tabular}",15)){
       inside_array = 1;
       strcpy(b_in,"\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|} \n");
        }
    else if(!strncmp(b_in,"\\end{tabular}",13)){
       inside_array = 0;
       }
    else if(inside_array && (b_in[0] != '%' && b_in[1] != '%')
            && strncmp(b_in,"\\hline",6)) {
       if(!line_is_opened) {
         strcpy(b_data, b_in);
/* Fill the data array with the next line */
       } else {
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
printf("\n astrom_read_measures/New line: >%s<\n", b_data);
printf(" astrom_read_measures/Trying to add a new object, current nobj=%d (wds_data=%d)\n",
       *nobj, with_wds_data);
#endif
/* Try to add a new object and add the parameters of this object to obj: */
       if(with_wds_data == 1)
            status = astrom_add_new_object_with_wds_data(b_data, obj, nobj,
                                                         i_notes);
       else
            status = astrom_add_new_object_without_wds_data(b_data, obj, nobj);
       if(status == 0) {
/* Index of current object in OBJECT structure array: */
          i_obj = *nobj - 1;
          (obj[i_obj]).nmeas = 0;
          }


/* Try to add a new measure: */
// JLP2020 it was != 0 : why ??
       if((status == 0) && (*nobj > 1)) {
         status = astrom_check_measure(b_data, i_eyepiece, i_rho, i_drho,
                                       i_theta, i_dtheta);
#ifdef DEBUG
printf("astrom_read_measures/Adding new measurement for object #i_obj=%d nm=%d (nobj=%d, status=%d)\n",
         i_obj, (obj[i_obj]).nmeas, *nobj, status);
#endif
          if(status == 0) {
           astrom_add_new_measure(b_data, obj, i_obj, i_filename, i_date,
                           i_filter, i_eyepiece, i_rho, i_drho, i_theta,
                           i_dtheta, i_notes, comments_wanted);
           printf("astrom_read_measures/ nmeas=%d rho=%.2f theta=%.2f\n", 
                  obj[i_obj].measure[0].rho, (obj[i_obj]).nmeas, 
                  obj[i_obj].measure[0].theta);
           strcpy(obj_name1, (obj[i_obj]).name);
           epoch1 = (obj[i_obj]).measure[0].epoch;
           status = astrom_get_meas_from_file(filein2, obj_name1,
                       epoch1, i_filename, i_date, i_filter,
                       i_eyepiece, i_rho, i_drho, i_theta, i_dtheta, i_notes,
                       &rho2, &err_rho2, &theta2, &err_theta2);
            } // case status from check_measure = 0
         }

       } /* EOF !line_is_opened */
    } // EOF line tabular
  } /* EOF if fgets() */
} /* EOF while loop */
return(0);
}
/*********************************************************************
*
*********************************************************************/
static int astrom_get_meas_from_file(char *filein2, char *obj_name1,
                       double epoch1, 
                       int i_filename, int i_date, int i_filter,
                       int i_eyepiece, int i_rho, int i_drho,
                       int i_theta, int i_dtheta, int i_notes,
                       double *rho2, double *err_rho2, double *theta2,
                       double *err_theta2)
{
FILE *fp_in2;

if((fp_in2 = fopen(filein2,"r")) == NULL) {
  fprintf(stderr, " Fatal error opening input file2 %s \n", filein2);
  return(-1);
  }

fclose(fp_in2);
return(0);
}

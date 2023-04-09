/*************************************************************************
* Program convert_for_WDS 
* To convert a calibrated table to be suited to direct includion in the WDS 
*
00004+2749  &  TDS1238  & 2013.963 & 1 & 0.834 & 0.007 & 266.2\rlap{$^*$} & 0.8 &  &\\
00005+2031  &  COU444  & 2013.924 & 1 & 0.693 & 0.007 & 39.5  & 0.3 & 3.25 &\\
00010+2721  &  DAM361Aa,Ab & 2013.927 & 1 & 1.442 & 0.007 & 358.3  & 0.3 & 2.26 & NDp\\
*
wds000.new:00004+2749       2013.963   q 86.2    0.8      0.834    0.007      .     .       .     .                0.8   1 Gii2022  S    7
wds000.new:00005+2031       2013.924     39.5    0.3      0.693    0.007      .     .      3.25   .                0.8   1 Gii2022  S  X 7
wds000.new:00010+2721 Aa,Ab 2013.927    358.3    0.3      1.442    0.007      .     .      2.26   .                0.8   1 Gii2022  S  X 7
* JLP
* Version 02/05/2009
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>   /* exit() */
#include <ctype.h>   // isdigit() 
#include <math.h>
#include <string.h>
#include "latex_utils.h" // jlp: latex_read_fvalue...
#include "jlp_string.h"    // jlp_trim_string, jlp_compact_string
#include "tex_calib_utils.h" //extract_companion_from 

/* Maximum number of columns to be extracted: */
#define IMAX 10 
 
/*
#define DEBUG
*/

static int jlp_convert_table_for_WDS(char *filein,char *fileout);
static int jlp_convert_line_for_WDS(char *in_line, char *out_line, 
                                    int line_length);
static int jlp_latex_to_ascii(FILE *fp_in, FILE *fp_out, int ix, int *iy,
                              int ncols, int icol_name);
static int jlp_latex_table_to_ascii(FILE *fp_in, FILE *fp_out, int ix, int *iy,
                              int ncols, int icol_name);
int JLP_RDLATEX_TABLE(int ix, int iy, float *xx, float *yy, int *npts,
                      int idim, char *filename);

int main(int argc, char *argv[])
{
char filein[60], fileout[60];
int status;

  printf("latex_to_ascii/ JLP/ Version 20/12/2010\n");
  printf("Note that this program can handle multiple LaTeX tables\n\n");

if(argc == 7 && *argv[4]) argc = 5;
if(argc == 7 && *argv[3]) argc = 4;
if(argc == 7 && *argv[2]) argc = 3;
if(argc == 7 && *argv[1]) argc = 2;
if(argc != 3)
  {
  printf("Error: argc=%d\n\n", argc);
  printf("Syntax: convert_for_WDS calibrated_table out_ascii_file \n");
  exit(-1);
  }
else
  {
  strcpy(filein,argv[1]);
  strcpy(fileout,argv[2]);
  }

printf(" OK: filein=%s fileout=%s \n", 
         filein, fileout);

/* Scan the file and make the conversion: */
status = jlp_convert_table_for_WDS(filein,fileout);

return(status);
}
/***********************************************************************
*
*************************************************************************/
static int jlp_convert_table_for_WDS(char *filein,char *fileout)
{
char in_line[512], out_line[512];
int iline;
FILE *fp_in, *fp_out;

if((fp_in = fopen(filein,"r")) == NULL)
{
printf(" Fatal error opening input file %s \n",filein);
return(-1);
}

if((fp_out = fopen(fileout,"w")) == NULL)
{
printf(" Fatal error opening output file %s \n",fileout);
fclose(fp_in);
return(-1);
}
iline = 0;
while(!feof(fp_in))
{
  if(fgets(in_line, 512, fp_in))
  {
  iline++;
  if(isdigit(in_line[0]) != 0) {
   jlp_convert_line_for_WDS(in_line, out_line, 512);
   fprintf(fp_out,"%s\n", out_line);
   }
  }
}

fclose(fp_in);
fclose(fp_out);
return(0);
}
/***********************************************************************
* Example:
00004+2749  &  TDS1238  & 2013.963 & 1 & 0.834 & 0.007 & 266.2\rlap{$^*$} & 0.8 &  &\\
00005+2031  &  COU444  & 2013.924 & 1 & 0.693 & 0.007 & 39.5  & 0.3 & 3.25 &\\
00010+2721  &  DAM361Aa,Ab & 2013.927 & 1 & 1.442 & 0.007 & 358.3  & 0.3 & 2.26 & NDp\\
*
wds000.new:00004+2749       2013.963   q 86.2    0.8      0.834    0.007      .     .       .     .                0.8   1 Gii2022  S    7
wds000.new:00005+2031       2013.924     39.5    0.3      0.693    0.007      .     .      3.25   .                0.8   1 Gii2022  S  X 7
wds000.new:00010+2721 Aa,Ab 2013.927    358.3    0.3      1.442    0.007      .     .      2.26   .                0.8   1 Gii2022  S  X 7
*************************************************************************/
static int jlp_convert_line_for_WDS(char *in_line, char *out_line, 
                                    int line_length)
{
int status = 0, verbose_if_error = 0, icol, object_len0, i;
char wds_name0[64], object_name0[64], comp_name0[64], buffer[54], q_flag[1], *pc;
double epoch0, rho0, drho0, theta0, dtheta0, dmag0;
/*
* icol=1: wds_name
* icol=2: object_name including companion_name
* icol=3: epoch 
* icol=4: ibin 
* icol=5: rho 
* icol=6: drho 
* icol=7: theta
* icol=8: dtheta
* icol=9: dmag
*/
  icol = 1;
  status = latex_read_svalue(in_line, buffer, icol);
// Remove $ if present (for negative values: $-$):
  pc = buffer;
  i = 0;
  while(*pc) {
    if(*pc != '$') wds_name0[i++] = *pc; 
    pc++;
    }
  wds_name0[i] = '\0';
  jlp_compact_string(wds_name0, 64);
  icol = 2;
  status = latex_read_svalue(in_line, buffer, icol);
  extract_companion_from_name(buffer, object_name0, comp_name0, &object_len0);
  icol = 3;
  status = latex_read_dvalue(in_line, &epoch0, icol, verbose_if_error);
  icol = 5;
  status = latex_read_dvalue(in_line, &rho0, icol, verbose_if_error);
  status = latex_read_svalue(in_line, buffer, icol);
  if(strstr(buffer,"rlap") != NULL) 
    q_flag[0] = 'q';
  else
    q_flag[0] = ' ';
  icol = 6;
  status = latex_read_dvalue(in_line, &drho0, icol, verbose_if_error);
  icol = 7;
  status = latex_read_dvalue(in_line, &theta0, icol, verbose_if_error);
  icol = 8;
  status = latex_read_dvalue(in_line, &dtheta0, icol, verbose_if_error);
  icol = 9;
  status = latex_read_dvalue(in_line, &dmag0, icol, verbose_if_error);
  if(status != 0) dmag0 = 0.;
  if(dmag0 > 0.) {
    sprintf(out_line, "wds000.new:%s %5.5s %8.3f  %c%6.1f    %2.1f     %6.3f    %5.3f      .     .      %4.2f   .                0.8   1 Gii2022  S  X 7",
            wds_name0, comp_name0, epoch0, q_flag[0], theta0, dtheta0, rho0, 
            drho0, dmag0); 
    } else {
    sprintf(out_line, "wds000.new:%s %5.5s %8.3f  %c%6.1f    %2.1f     %6.3f    %5.3f      .     .       .     .                0.8   1 Gii2022  S    7",
            wds_name0, comp_name0, epoch0, q_flag[0], theta0, dtheta0, rho0, 
            drho0); 
    }
return(0);
}
/*************************************************************************
*
* INPUT:
* ix, iy1, iy2, iy3: column numbers for xx and yy1, yy2, yy3, yy4 
*************************************************************************/
static int jlp_latex_table_to_ascii(FILE *fp_in, FILE *fp_out, int ix, int *iy,
                                    int ncols, int icol_name) 
{
float xx, yy[IMAX];
char in_line[NMAX], b_data[NMAX], name[80];
int inside_array, line_is_opened, status, iline, i;
int verbose_if_error = 0;
char *pc, *pc1;

inside_array = 0;
line_is_opened = 0;

/* Possibility of storing the index in first column, when ix=0 */
if(ix <= 0) xx = 0.;

iline = 0;
while(!feof(fp_in))
{
/* Maximum length for a line will be NMAX/2 = 180 characters: */
  if(fgets(in_line, NMAX/2, fp_in))
  {
  iline++;
  in_line[NMAX/2] = '\0';
/* WARNING: since 2007, \begin{tabular*}
*                  instead of \begin{tabular}:
*/
    if(!strncmp(in_line,"\\begin{tabular",14)){
       inside_array = 1;
#ifdef DEBUG
printf(" OK: >%s< inside_array=%d\n", in_line, inside_array);
#endif
        }
/* WARNING: since 2007, \end{tabular*}
*                  instead of \end{tabular}:
*/
    else if(!strncmp(in_line,"\\end{tabular",12)){
       inside_array = 0;
#ifdef DEBUG
printf(" OK: >%s< inside_array=%d\n", in_line, inside_array);
#endif
       }
    else if(inside_array && in_line[0] != '%' && strncmp(in_line,"\\hline",6)) {
       if(!line_is_opened) {
         strcpy(b_data, in_line);
/* Fill the data array with the next line */
       } else {
/* Look for the first zero (end of string marker) in data buffer */
         b_data[NMAX/2] = '\0';
         pc1 = b_data;
         while(*pc1) pc1++; 
         pc1--; 
/* Then copy the second line from there*/
         strcpy(pc1, in_line);
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

/* Processed data when line is closed: */
     if(!line_is_opened) {

#ifdef DEBUG
       printf("%s\n", b_data);
#endif

/* Read iy[0] column first (since if ix=0, reading of xx is always good) */
       status = latex_read_fvalue(b_data, &yy[0], iy[0], verbose_if_error);
/* If successful reading of yy1, continue and try to read ix column: */
       if((status == 0) || (status == 3)) {
/* If ix is not zero, read the value of this column  and load it to xx */
         if(ix > 0) {
           status = latex_read_fvalue(b_data, &xx, ix, verbose_if_error);
/* Otherwise, increment xx which is simply the current index */
         } else {
           xx++;
           status = 0;
         }
/* Read next y columns now */
       for(i = 1; i < ncols && (status != 0 || status != 3); i++) {
         status = latex_read_fvalue(b_data, &yy[i], iy[i], verbose_if_error);
        }
       if(icol_name && (status != 0 && status != 3)) status = latex_read_svalue(b_data, name, icol_name);
       if(status == 0 || status == 3) {
         fprintf(fp_out,"%.3f %.3f", xx, yy[0]); 
         for(i = 1; i < ncols; i++) fprintf(fp_out," %.3f", yy[i]); 

         if(icol_name)
           fprintf(fp_out,"%s \n", name); 
         else
           fprintf(fp_out,"\n"); 

         }
       } /* EOF !status */
      } /* EOF case line_is_closed (i.e., !line_is_opened) */
    } /* EOF case inside array and line not starting with \hline or % */
  } /* EOF successful reading of buffer from the input file */
} /* EOF while !feof(fp_in) */
return(0);
}
/*************************************************************************
* Input file: plain ascii table (possibly without header)
*
* INPUT:
* ix, iy1, iy2, iy3: column numbers for xx and yy1, yy2, yy3, yy4 
*************************************************************************/
static int jlp_latex_to_ascii(FILE *fp_in, FILE *fp_out, int ix, int *iy,
                              int ncols, int icol_name) 
{
float xx, yy[IMAX];
char in_line[NMAX], b_data[NMAX], name[80];
int status, iline, i, verbose_if_error = 0;
char *pc, *pc1;

#ifdef DEBUG
printf(" ncols=%d\n", ncols);
for(i = 0; i < ncols; i++) printf("i=%d, iy=%d \n", i, iy[i]);
#endif

iline = -1;
xx = 0.;
while(!feof(fp_in))
{
/* Maximum length for a line will be NMAX/2 = 512 characters: */
  if(fgets(in_line, NMAX/2, fp_in))
  {
  iline++;
  in_line[NMAX/2] = '\0';
  strcpy(b_data, in_line);
#ifdef DEBUG
printf(" iline=%d Data line: >%s<\n", iline, b_data);
#endif

/* Read iy[0] column first (since if ix=0, reading of xx is always good) */
       status = latex_read_fvalue(b_data, &yy[0], iy[0], verbose_if_error);
/* If successful reading of yy1, continue and try to read ix column: */
       if(status == 0) {
/* If ix is not zero, read the value of this column  and load it to xx */
         if(ix > 0) {
           status = latex_read_fvalue(b_data, &xx, ix, verbose_if_error);
         } else {
/* Otherwise, xx which is simply incremented */
           xx++;
         }
       }
/* Read next y columns now */
       for(i = 1; (i < ncols) && (status == 0); i++) {
         status = latex_read_fvalue(b_data, &yy[i], iy[i], verbose_if_error);
        }
       if((icol_name != 0) && (status == 0)) status = latex_read_svalue(b_data, name, icol_name);
       if(status == 0) {
#ifdef DEBUG_
       for(i = 1; i < ncols; i++) printf(" Output: i=%d yy=%f\n", i, yy[i]);
#endif
         fprintf(fp_out,"%.8g %.8g", xx, yy[0]); 
         for(i = 1; i < ncols; i++) fprintf(fp_out," %.8g", yy[i]); 
         if(icol_name)
           fprintf(fp_out,"%s \n", name); 
         else
           fprintf(fp_out,"\n"); 
         } /* EOF !status */
   } /* EOF successful reading of buffer from the input file */
} /* EOF while !feof(fp_in) */
return(0);
}
/**************************************************************************
* JLP_RDLATEX_TABLE
* Interface with Fortran programs (not finished yet)
*
* INPUT
* ix, iy: column numbers for xx and yy arrays 
* idim: maximum number of data points
* filename: name of the input LaTeX file
*
* OUTPUT:
* xx, yy: data points
* npts: number of data points
*
***************************************************************************/
int JLP_RDLATEX_TABLE(int ix, int iy, float *xx, float *yy, int *npts,
                      int idim, char *filename)
{
FILE *fp;

if((fp = fopen(filename,"r")) == NULL) {
  return(-1);
  }

fclose(fp);
return(0);
}

/*************************************************************************
* Program latex_to_ascii
* To convert a LaTeX table to an ASCII list 
*
* JLP
* Version 02/05/2009
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>   /* exit() */
#include <math.h>
#include <string.h>
#include "latex_utils.h" // jlp: latex_read_fvalue...

/* Maximum number of columns to be extracted: */
#define IMAX 10 
 
/*
#define DEBUG
*/

static int jlp_latex_to_ascii(FILE *fp_in, FILE *fp_out, int ix, int *iy,
                              int ncols, int icol_name);
static int jlp_latex_table_to_ascii(FILE *fp_in, FILE *fp_out, int ix, int *iy,
                              int ncols, int icol_name);
int JLP_RDLATEX_TABLE(int ix, int iy, float *xx, float *yy, int *npts,
                      int idim, char *filename);

int main(int argc, char *argv[])
{
char filein[60], fileout[60];
int ix, iy[IMAX], icol_name, ncols, i;
FILE *fp_in, *fp_out;

  printf("latex_to_ascii/ JLP/ Version 20/12/2010\n");
  printf("Note that this program can handle multiple LaTeX tables\n\n");

if(argc == 7 && *argv[4]) argc = 5;
if(argc == 7 && *argv[3]) argc = 4;
if(argc == 7 && *argv[2]) argc = 3;
if(argc == 7 && *argv[1]) argc = 2;
if(argc != 4 && argc != 5)
  {
  printf("Error: argc=%d\n\n", argc);
  printf("Syntax: latex_to_ascii in_latex_table out_ascii_file ix,iy1,iy2,...,iy10 [icol_name]\n");
  printf("\n(Enter simply ix,iy1 for 2 columns and ix,iy1,iy2, for 3 columns)\n");
  printf("\n(Enter 0,ix to generate following list: (index, ix column)\n");
  exit(-1);
  }
else
  {
  strcpy(filein,argv[1]);
  strcpy(fileout,argv[2]);
  ix = 0;
  for(i = 0; i < IMAX; i++) iy[i] = 0;
  icol_name = 0;
  ncols = sscanf(argv[3],"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d", &ix, &iy[0], &iy[1], &iy[2], &iy[3], &iy[4], &iy[5], &iy[6], &iy[7], &iy[8], &iy[9]);
  ncols--;
  if(argc == 5) sscanf(argv[4],"%d", &icol_name);
  }

printf(" OK: filein=%s fileout=%s ncols=%d\n ix=%d iy1, iy2,...=%d %d %d %d %d %d %d %d %d %d icol_name=%d\n", 
         filein, fileout, ncols, ix, iy[0], iy[1], iy[2], iy[3], iy[4], iy[5], iy[6], 
         iy[7], iy[8], iy[9], icol_name);

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
fprintf(fp_out,"%% From %s ix=%d iy1,2,...=", filein, ix);
for(i = 0; i < ncols; i++) fprintf(fp_out," %d", iy[i]);
fprintf(fp_out," icol_name=%d\n", icol_name);

/* Scan the file and make the conversion: */
jlp_latex_table_to_ascii(fp_in,fp_out, ix, iy, ncols, icol_name);

fclose(fp_in);
fclose(fp_out);
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

/*************************************************************************
* Program latex_to_gili_residuals
* To convert a LaTeX table to the list of residuals 
*
* JLP
* Version 30/10/2019
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>   /* exit() */
#include <math.h>
#include <string.h>

/* Maximum length for one line will be NMAX=180 characters: */
/* Two lines are allowed */
#define NMAX 360

/* Maximum number of columns to be extracted: */
#define IMAX 10 
 
/*
#define DEBUG
*/

static int read_fvalue(char *b_data, float *value, int icol);
static int read_svalue(char *b_data, char *value, int icol);
static int jlp_latex_to_gili_residuals(FILE *fp_in, FILE *fp_out, 
                                       int icol_name);
static int decode_gili_residual(char *b_data, double *theta, double *rho, 
                                char *orbit_ref);

int main(int argc, char *argv[])
{
char filein[60], fileout[60];
int icol_name, i;
FILE *fp_in, *fp_out;

  printf("latex_to_gili_residuals/ JLP/ Version 30/10/2019\n");

if(argc == 7 && *argv[4]) argc = 5;
if(argc == 7 && *argv[3]) argc = 4;
if(argc == 7 && *argv[2]) argc = 3;
if(argc == 7 && *argv[1]) argc = 2;
if(argc != 3 && argc != 4)
  {
  printf("Error: argc=%d\n\n", argc);
  printf("Syntax: latex_to_gili_residuals in_latex_table out_ascii_file\n");
  exit(-1);
  }
else
  {
  strcpy(filein,argv[1]);
  strcpy(fileout,argv[2]);
  }

// Third column of ascii input file contains the name of the binary
icol_name = 3;
printf(" OK: filein=%s fileout=%s icol_name=%d\n", filein, fileout, icol_name);

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
fprintf(fp_out,"%% From %s ", filein);
fprintf(fp_out," icol_name=%d\n", icol_name);

/* Scan the file and make the conversion: */
jlp_latex_to_gili_residuals(fp_in,fp_out, icol_name);

fclose(fp_in);
fclose(fp_out);
return(0);
}
/*************************************************************************
*
* INPUT:
*************************************************************************/
static int jlp_latex_to_gili_residuals(FILE *fp_in, FILE *fp_out, 
                                       int icol_name) 
{
double theta, rho;
char in_line[NMAX], b_data[NMAX], name[80], orbit_ref[64];
int inside_array, line_is_opened, status, iline, i;
char *pc, *pc1;

inside_array = 0;
line_is_opened = 0;

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
       printf("%s\n", b_data);
       status = decode_gili_residual(b_data, &theta, &rho, orbit_ref);
       if(status == 0) {
// Decode the icol_name column and put it at the end of each line:
       if(icol_name && !status) status = read_svalue(b_data, name, icol_name);
       if(!status) {
         fprintf(fp_out," %.2f %.3f %s ", theta, rho, orbit_ref); 
         if(icol_name)
           fprintf(fp_out," %s \n", name); 
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
/**************************************************************************
* Decode current line (b_data) 
* and extract the residuals theta, rho and orbit_ref 
*
**************************************************************************/
static int decode_gili_residual(char *b_data, double *theta, double *rho, 
                                char *orbit_ref)
{
// 12th line contains the residuals
int status = -1, icol = 12;
char cvalue[64];

// Example:
// 19.0  0.021 Hrt2000c
  status = read_svalue(b_data, cvalue, icol);
  if(status == 0) {
    printf(" icol=%d cvalue= %s \n", icol, cvalue);
    if(sscanf(cvalue," %lf %lf", theta, rho) == 2) {
    printf(" theta=%f rho=%f \n", *theta, *rho);
    sscanf(cvalue," %*f %*f %s", orbit_ref);
    } else {
    status = -1;
    }
  }
return(status);
}
/**************************************************************************
* Read float value in column #icol from b_data string
*
**************************************************************************/
static int read_fvalue(char *b_data, float *value, int icol)
{
int nval, status;
char buff[40], *pc;

*value = 0.;
status = read_svalue(b_data, buff, icol);

/* Clear string if it is empty*/
if(!status) {
pc = buff;
while(*pc){
   if(*pc != ' ') break;
   pc++;
   }
if(*pc == '\0' && *(--pc) == ' ') buff[0] = '\0';
/* */

if(!*buff) status = -1;
}
 
if(!status) {
   nval = sscanf(buff, "%f", value);
   if(nval != 1) {
      fprintf(stderr, "read_fvalue/buff=>%s< value=%.2f nval=%d\n", 
              buff, *value, nval);
      status = 1;
      }
  }
if(!status) printf(">%s< value=%f\n", buff, *value);

return(status);
}
/**************************************************************************
* Read string value in column #icol from b_data string
*
**************************************************************************/
static int read_svalue(char *b_data, char *value, int icol)
{
int ic, status, column_is_found;
char buff[NMAX], data[NMAX], *pc, prev_pc[1];

strcpy(data, b_data);

pc = data;
data[NMAX-1] = '\0';
column_is_found = 0;
ic = 1;
buff[0] = '\0';
*prev_pc = '\0';
while(*pc) {
  if(ic == icol) {
    column_is_found = 1;
    strcpy(buff,pc);
    break;
    }
/* JLP 2008: Should avoid "\&" that may be present (bibliography for instance)
*/
  if(*pc == '&' && *prev_pc != '\\') {
    ic++;
    }
  *prev_pc = *pc;
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


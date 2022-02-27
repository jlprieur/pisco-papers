/*************************************************************************
* Program gaia_add_parallax
* 
* To process output gaia parallax and update Latex array with measurements 
*
\hline 
00022+2705 & 000210.18+270455.6 & BU733AB & G5Vb+K5V & 5.83 & 8.9 & & & 2.00 \\
*
* JLP
* Version 22/01/2022
*************************************************************************/
#include <stdio.h>
#include <stdlib.h> // exit(-1)
#include <string.h>
#include <ctype.h>  // isalpha(), isdigit()
#include "latex_utils.h" // latex_get_column_item(), latex_remove_column()
#include "jlp_string.h"  // (in jlplib/jlp_fits/ ) jlp_cleanup_..

#include "latex_utils.h"  // latex_get_column_item()
//#include "astrom_utils1.h" 
//#include "astrom_utils2.h" 

#define MAX_NLINES 1024
#define MAX_LENGTH 1024

#define DEBUG

static int csv_gaia_read_parallax(char *gaia_csv_data, int in_target_col,
                                  int in_plx_col, 
                                  int in_err_plx_col, double *plx_array,
                                  double *err_plx_array, char *target_label, 
                                  int *n_plx_array, int max_n_array);
static int gaia_add_parallax(FILE *fp_in, FILE *fp_out,
                             double *plx_array, double *err_plx_array, 
                             char *target_label, int n_plx_array,
                             int out_plx_col, int out_err_plx_col);
static int jlp_csv_get_item(char *in_line, int icol, char *item);
static int jlp_remove_dollar(char *in_string, char *out_string, int len_string); 
static int jlp_find_target_index(char *target_label, char *target_name, 
                                 int ntargets, int *indx);

int main(int argc, char *argv[])
{
int i, status, in_target_col, in_plx_col, out_plx_col; 
int in_err_plx_col, out_err_plx_col, n_plx_array; 
char gaia_csv_data[128], filein[128], fileout[128], target_label[MAX_NLINES*32];
FILE *fp_in, *fp_out;
double plx_array[MAX_NLINES], err_plx_array[MAX_NLINES];

/* If command line with "runs" */
if(argc == 7){
 if(*argv[5]) argc = 6;
 else if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
 else argc = 1;
 }

if(argc != 5)
  {
  printf(" Syntax: gaia_add_parallax in_gaia_data in_latex_file out_latex_file in_target_col,in_plx_col,in_err_plx_col,out_plx_col,out_err_plx_col\n");
  exit(-1);
  }
else
  {
  strcpy(gaia_csv_data,argv[1]);
  strcpy(filein,argv[2]);
  strcpy(fileout,argv[3]);
  sscanf(argv[4], "%d,%d,%d,%d,%d", &in_target_col, &in_plx_col, 
                 &in_err_plx_col, &out_plx_col, &out_err_plx_col);
  }

printf(" OK: gaia_csv_data=%s \n", gaia_csv_data);
printf(" OK: filein=%s fileout=%s \n", filein, fileout);
printf("in_target_col=%d, in_plx_col=%d out_plx_col=%d \n", 
        in_target_col, in_plx_col, out_plx_col);
printf("in_err_plx_col=%d out_err_plx_col=%d \n", in_err_plx_col, out_err_plx_col);

if((fp_in = fopen(filein,"r")) == NULL) {
  printf(" Fatal error opening input file %s \n",filein);
  exit(-1);
  }

if((fp_out = fopen(fileout,"w")) == NULL) {
  printf(" Fatal error opening output file %s \n",fileout);
  fclose(fp_in);
  exit(-1);
  }

csv_gaia_read_parallax(gaia_csv_data, in_target_col, in_plx_col, in_err_plx_col,
                   plx_array, err_plx_array, target_label, &n_plx_array, 
                   MAX_NLINES);

gaia_add_parallax(fp_in, fp_out, plx_array, err_plx_array, target_label,
                  n_plx_array, out_plx_col, out_err_plx_col);

fclose(fp_in);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input csv gaia data array
************************************************************************/
static int csv_gaia_read_parallax(char *gaia_csv_data, int in_target_col,
                                  int in_plx_col, 
                                  int in_err_plx_col, double *plx_array,
                                  double *err_plx_array, char *target_label, 
                                  int *n_plx_array, int max_n_array)
{
FILE *fp_gaia;
char in_line[MAX_LENGTH], in_item[MAX_LENGTH];
int status, iline, i;

if((fp_gaia = fopen(gaia_csv_data,"r")) == NULL) {
  printf(" Fatal error opening input file %s \n", fp_gaia);
  exit(-1);
  }

iline = -1;
while(!feof(fp_gaia)) {
  if(fgets(in_line, MAX_LENGTH, fp_gaia)) {
// Header is starting with a character which is not a digit:
     if(isdigit(in_line[0]) != 0) {
       iline++;
       plx_array[iline] = 0.;
       err_plx_array[iline] = 0.;
#ifdef DEBUG
       printf("\n OK0: iline=%d in_line=%s", iline, in_line);
#endif
// WDSxxx in in_target_col:
       status = jlp_csv_get_item(in_line, in_target_col, in_item);
       if(status == 0) {
         sprintf(&target_label[32 * iline], in_item);
         }
// parallax in in_plx_col:
       status = jlp_csv_get_item(in_line, in_plx_col, in_item);
       if(status == 0) {
         sscanf(in_item, "%lf", &plx_array[iline]);
         }
// parallax_error in in_err_plx_col:
       status = jlp_csv_get_item(in_line, in_err_plx_col, in_item);
       if(status == 0) {
         sscanf(in_item, "%lf", &err_plx_array[iline]);
         }
#ifdef DEBUG
       printf("OK: iline=%d target=%s plx_array=%f err_plx_array=%f \n", 
             iline, &target_label[32*iline], plx_array[iline], 
             err_plx_array[iline]);
#endif
       if(iline == (max_n_array-1) ) break;
       } // if isdigit...
    } // if(fgets...
  } // while(!feof...)
*n_plx_array = iline;

fclose(fp_gaia);

return(0);
} 
/************************************************************************
* Scan the input table and make the modifications
*
* INPUT:
* fp_in: pointer to the input file containing the input table
* fp_out: pointer to the output Latex file
*
*************************************************************************/
static int gaia_add_parallax(FILE *fp_in, FILE *fp_out,
                             double *plx_array, double *err_plx_array, 
                             char *target_label, int n_plx_array,
                             int out_plx_col, int out_err_plx_col)
{
char in_line[MAX_LENGTH], in_line3[MAX_LENGTH];
char buffer[MAX_LENGTH], *pc, plx_string[MAX_LENGTH], err_plx_string[MAX_LENGTH];
char target_str[128];
int iline, status, verbose_if_error = 0, string_len;
int i, in_line_length, indx;

iline = -1;
while(!feof(fp_in)) {
  if(fgets(in_line, MAX_LENGTH, fp_in)) {
  if(isdigit(in_line[0]) != 0) {
    iline++;
// Remove the end of line '\n' from input line:
//    jlp_cleanup_string(in_line, MAX_LENGTH);

// Get the name in the 1st column 
   status = latex_get_column_item(in_line, buffer, 1, verbose_if_error);
   jlp_remove_dollar(buffer, target_str, 64);
   status = jlp_find_target_index(target_label, target_str, n_plx_array, &indx);

   in_line_length = MAX_LENGTH;

   sprintf(plx_string, "\\nodata ");
   sprintf(err_plx_string, "\\nodata ");
   if((status == 0) && (indx >= 0)) {
     printf("indx=%d TARGET: %s \n", indx, &target_label[32 * indx]);
     if(plx_array[indx] > 0.0) {
       sprintf(plx_string, "%.3f ", plx_array[indx]);
       sprintf(err_plx_string, "%.3f ", err_plx_array[indx]);
       }
    }
    string_len = strlen(plx_string);
    status = latex_set_column_item(in_line, in_line_length, plx_string,
                                   string_len, out_plx_col, verbose_if_error); 
    string_len = strlen(err_plx_string);
    status = latex_set_column_item(in_line, in_line_length, err_plx_string,
                                   string_len, out_err_plx_col, 
                                   verbose_if_error); 
   } // isdigit

// Save to output file:
   fprintf(fp_out, "%s", in_line);

  } /* EOF if fgets */
 } /* EOF while ... */
printf("gaia_add_parallax: %d lines sucessfully read and processed\n",
        iline);
return(0);
}
/*************************************************************************
*
*************************************************************************/
static int jlp_csv_get_item(char *in_line, int icol, char *item)
{
char *pc, *pc1, buffer[MAX_LENGTH], item_buffer[MAX_LENGTH];
int status = -1, nitem, istart;

strcpy(item, "");
strcpy(buffer, in_line);
pc = buffer;
istart = 0;
nitem = 1;
while(*pc) {
// istart is setto 1 when nitem is good for the first time:
 if((nitem == icol) && (istart == 0)) {
  strcpy(item_buffer, pc);
  pc1 = item_buffer; 
  istart = 1;
  }
 if(*pc == ',') {
   if(istart != 0) {
     *pc1 = '\0'; 
     status = 0;
     break;
     }
   nitem++;
   }
 pc++;
 if(istart != 0) pc1++;
 }

strcpy(item, item_buffer);

#ifdef DEBUG
  printf("icol=%d item= >%s<\n", icol, item);
#endif
return(status);
}
/*************************************************************************
* Remove dollar in input string
*************************************************************************/
static int jlp_remove_dollar(char *in_string, char *out_string, int len_string)
{
int status = 0, i, j;

j = 0;
for(i = 0; i < len_string; i++) { 
  if(in_string[i] == '\0') 
    break; 
// Remove dollar and also blank space...
  else if((in_string[i] != '$') && (in_string[i] != ' ')) 
    out_string[j++] = in_string[i];
  }
out_string[j] = '\0';

#ifdef DEBUG
// printf("remove_dollar/ in: >%s< out= >%s<\n", in_string, out_string);
#endif
return(status);
}
/*************************************************************************
* 
*************************************************************************/
static int jlp_find_target_index(char *target_label, char *target_name, 
                                 int ntargets, int *indx)
{
char targt[128], buffer[128];
int i, status = -1;

*indx = -1;

jlp_compact_string(target_name, 128);

  for(i = 0; i < ntargets; i++) {
    sprintf(buffer, &target_label[32 * i]);
    jlp_compact_string(buffer, 128);
// "WDS23167+2959"
    strcpy(targt, &buffer[4]);
    targt[10] = '\0';
    if(!strcmp(targt, target_name)) {
#ifdef DEBUG
//   printf("jlp_find_target_index %s == %s \n", targt, target_name);
#endif
      *indx = i;
      status = 0;
      break;
      }
    }
return(status);
}

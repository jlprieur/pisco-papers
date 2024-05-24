/*************************************************************************
* Program gaia_plot_hr
* 
* Process input gaia data with parallax and t_eff and create latex graphic file
*
* JLP
* Version 08/02/2022
*************************************************************************/
#include <stdio.h>
#include <stdlib.h> // exit(-1)
#include <math.h>  // log10
#include <string.h>
#include <ctype.h>  // isalpha(), isdigit()
#include "latex_utils.h" // latex_get_column_item(), latex_remove_column()
#include "jlp_string.h"  // (in jlplib/jlp_fits/ ) jlp_cleanup_..

//#include "astrom_utils1.h" 
//#include "astrom_utils2.h" 

#define MAX_NLINES 1024
#define MAX_LENGTH 1024

#define DEBUG

static int csv_gaia_read_data(char *gaia_csv_data, int in_target_col,
                              int in_plx_col, int in_err_plx_col, 
                              int in_gmag_col, int in_absp_col, 
                              int in_teff_col, 
                              double *plx_array, double *err_plx_array, 
                              double *gmag_array, double *absp_array,
                              double *teff_array, char *target_label, 
                              int *n_plx_array, double rel_plx_error_maxi,
                              FILE *fp_out, int max_n_array);
static int gaia_add_parallax(FILE *fp_in, FILE *fp_out,
                             double *plx_array, double *err_plx_array, 
                             char *target_label, int n_plx_array,
                             int out_plx_col, int out_err_plx_col);
// Now in csv_utils.h csv_read_string():
static int jlp_csv_get_item(char *in_line, int icol, char *item);
static int jlp_remove_dollar(char *in_string, char *out_string, int len_string); 
static int jlp_find_target_index(char *target_label, char *target_name, 
                                 int ntargets, int *indx);
static int compute_abs_mag(double gmag, double absp, double plx, double *Gmag);

int main(int argc, char *argv[])
{
int i, status, in_target_col, in_plx_col, out_plx_col; 
int in_err_plx_col, in_gmag_col, in_absp_col, in_teff_col, n_plx_array; 
char gaia_csv_data[128], fileout_latex[128]; 
char target_label[MAX_NLINES*32];
FILE *fp_out;
double plx_array[MAX_NLINES], err_plx_array[MAX_NLINES];
double gmag_array[MAX_NLINES], absp_array[MAX_NLINES], teff_array[MAX_NLINES];
double rel_plx_error_maxi;

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
  printf(" Syntax: gaia_plot_hr in_gaia_data out_latex_file in_target_col,in_plx_col,in_err_plx_col,in_gmag_col,in_absp_col,in_teff_col rel_plx_error_maxi_max\n");
  exit(-1);
  }
else
  {
  strcpy(gaia_csv_data, argv[1]);
  strcpy(fileout_latex, argv[2]);
  sscanf(argv[3], "%d,%d,%d,%d,%d,%d", &in_target_col, &in_plx_col, 
                 &in_err_plx_col, &in_gmag_col, &in_absp_col, &in_teff_col);
  sscanf(argv[4], "%lf", &rel_plx_error_maxi); 
  }

printf(" OK: gaia_csv_data=%s \n", gaia_csv_data);
//if(argc != 1234) exit(-1);
printf(" OK: fileout_latex=%s \n", fileout_latex);
printf("in_target_col=%d, in_plx_col=%d in_err_plx_col=%d \n",
        in_target_col, in_plx_col, in_err_plx_col);
printf("in_gmag_col=%d in_absp_col=%d in_teff_col=%d \n", 
        in_gmag_col, in_absp_col, in_teff_col);
printf("rel_plx_error_maxi=%f \n", rel_plx_error_maxi); 

if((fp_out = fopen(fileout_latex,"w")) == NULL) {
  printf(" Fatal error opening output file %s \n",fileout_latex);
  exit(-1);
  }

status = csv_gaia_read_data(gaia_csv_data, in_target_col, in_plx_col, 
                            in_err_plx_col, in_gmag_col, in_absp_col, 
                            in_teff_col,
                            plx_array, err_plx_array, gmag_array, absp_array,
                            teff_array,
                            target_label, &n_plx_array, rel_plx_error_maxi,
                            fp_out, MAX_NLINES);


fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input csv gaia data array
************************************************************************/
static int csv_gaia_read_data(char *gaia_csv_data, int in_target_col,
                              int in_plx_col, int in_err_plx_col, 
                              int in_gmag_col, int in_absp_col, 
                              int in_teff_col, 
                              double *plx_array, double *err_plx_array, 
                              double *gmag_array, double *absp_array,
                              double *teff_array, char *target_label, 
                              int *n_plx_array, double rel_plx_error_maxi,
                              FILE *fp_out, int max_n_array)
{
FILE *fp_gaia;
char in_line[MAX_LENGTH], in_item[MAX_LENGTH], buffer[128], targt[128];
double Gmag, old_gmag, old_plx, relative_plx_error;
int status, iline, i;

// Prolog of the table
  fprintf(fp_out, "\\begin{table} \n");
  fprintf(fp_out, "\\begin{tabular}{llllll} \n");
  fprintf(fp_out, "\\hline \n");
  fprintf(fp_out, "WDS & $\\pi$ & err$_\\pi$ & $g_{mag}$ & $A$ & $T_{eff}$ & $G_{mag}$ \\\\ \n");
  fprintf(fp_out, " & (mas) & (mas) & (mag) & (mag) & (K) & (mag) \\\\ \n");
  fprintf(fp_out, "\\hline \n");

if((fp_gaia = fopen(gaia_csv_data,"r")) == NULL) {
  printf(" Fatal error opening input file %s \n", fp_gaia);
  exit(-1);
  }

iline = -1;
old_gmag = 0.;
while(!feof(fp_gaia)) {
  if(fgets(in_line, MAX_LENGTH, fp_gaia)) {
// Header is starting with a character which is not a digit:
     if(isdigit(in_line[0]) != 0) {
       iline++;
       plx_array[iline] = 0.;
       err_plx_array[iline] = 0.;
       gmag_array[iline] = 0.;
       absp_array[iline] = 0.;
       teff_array[iline] = 0.;
#ifdef DEBUG
       printf("\n OK0: iline=%d in_line=%s", iline, in_line);
#endif
// WDSxxx in in_target_col:
       targt[0] = '\0';
       status = jlp_csv_get_item(in_line, in_target_col, in_item);
       if(status == 0) {
         sprintf(&target_label[32 * iline], in_item);
         sprintf(buffer, in_item);
         jlp_compact_string(buffer, 128);
// "WDS23167+2959"
         strcpy(targt, &buffer[4]);
         targt[10] = '\0';
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
// gmag in in_gmag_col:
       status = jlp_csv_get_item(in_line, in_gmag_col, in_item);
       if(status == 0) {
         sscanf(in_item, "%lf", &gmag_array[iline]);
         }
// Absp in in_gmag_col:
       status = jlp_csv_get_item(in_line, in_absp_col, in_item);
       if(status == 0) {
         sscanf(in_item, "%lf", &absp_array[iline]);
         }
// teff in in_teff_col:
       status = jlp_csv_get_item(in_line, in_teff_col, in_item);
       if(status == 0) {
         sscanf(in_item, "%lf", &teff_array[iline]);
         }
#ifdef DEBUG
       printf("OK: iline=%d target=%s plx_array=%f err_plx_array=%f gmag=%f absp=%f teff=%f\n", 
             iline, &target_label[32*iline], plx_array[iline], 
             err_plx_array[iline], gmag_array[iline], absp_array[iline], 
             teff_array[iline]);
#endif
       relative_plx_error = 1000.;
       if((plx_array[iline] > 0.) && (err_plx_array[iline] > 0.)) {
         relative_plx_error = err_plx_array[iline] / plx_array[iline];
       }
// Save data if it is a new target and has a non null parallax:
       if((old_plx != plx_array[iline]) && (old_gmag != gmag_array[iline])
// Relative error smaller than 0.5 :
          && (relative_plx_error < rel_plx_error_maxi)
          && (gmag_array[iline] > 0.) && (teff_array[iline] > 0.) ) {
         status = compute_abs_mag(gmag_array[iline], absp_array[iline], 
                                  plx_array[iline], &Gmag);
         if(status == 0) {
           fprintf(fp_out, "%s & %.4f & %.4f & %.3f & %.3f & %.3f & %.3f \\\\ \n", 
                   targt, plx_array[iline], err_plx_array[iline], 
                   gmag_array[iline], absp_array[iline], teff_array[iline],
                   Gmag);
           }
         }
       old_plx = plx_array[iline];
       old_gmag = gmag_array[iline];
       if(iline == (max_n_array-1) ) break;
       } // if isdigit...
    } // if(fgets...
  } // while(!feof...)
*n_plx_array = iline;

// Epilog of the table
  fprintf(fp_out, "\\hline \n");
  fprintf(fp_out, "\\end{tabular} \n");
  fprintf(fp_out, "\\end{table} \n");

fclose(fp_gaia);

return(0);
} 
/**********************************************************************
* Compute absolute magnitude
**********************************************************************/
static int compute_abs_mag(double gmag, double absp, double plx, double *Gmag)
{
// Absolute magnitude M: 
// m - M = 5 log10(D) - 5
// D distance in parsecs : D = 1/pi
// Apparent magnitude: m = m_obs - A
// Interstellar absorption in the G band: A in magnitudes 
int status = -1;
 *Gmag = 0.;
 if(plx > 0.) {
// absp is positive:
printf("OKKK222: absp=%f\n", absp);
   *Gmag = gmag + absp - 5. * log10(1000. / plx) + 5.;
   status = 0;
   }
return(status);
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

/************************************************************************
* "modif_calib_table1.c"
*
* To modify the old formatted tables (MerI, MerII, MerIII, ... MerVI)
* by suppressing orbity column (11th column)
*
* JLP 
* Version 05/03/2020
*************************************************************************/
#include "jlp_catalog_utils.h"

#define DEBUG
#define DEBUG_1
/*
*/

static int modif_calib1(FILE *fp_calib, FILE *fp_out); 
static int remove_column(char *in_line, int i_col, int ilen);
static int remove_orbit_column(char *in_line, int ilen);

int main(int argc, char *argv[])
{
char calib_fname[80], out_fname[80];
FILE *fp_calib, *fp_out;
time_t t = time(NULL);

if(argc != 3) {
  printf("Syntax: modif_calib_table1 old_calib_table out_calib_table\n");
  return(-1);
}

strcpy(calib_fname, argv[1]);
strcpy(out_fname, argv[2]);

printf("OK: in_calib=%s output=%s \n", calib_fname, out_fname); 

/* Open input calibual table: */
if((fp_calib = fopen(calib_fname, "r")) == NULL) {
   fprintf(stderr, "modif_calib_table1/Fatal error opening calib table %s\n",
           calib_fname);
    return(-1);
  }

/* Open output table: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "modif_calib_table1/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the output Latex table: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "%% Modified table (orbit col. removed) from: %s \n%% Created by modif_calib_table1.c on %s", 
        calib_fname, ctime(&t));
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

/* Scan the input calibrated table 
*/
modif_calib1(fp_calib, fp_out); 

/* Close opened files:
*/
fclose(fp_calib);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input table and make the modifications 
*
* INPUT:
* fp_calib: pointer to the input file containing the input table 
* fp_out: pointer to the output Latex file
*
*************************************************************************/
static int modif_calib1(FILE *fp_calib, FILE *fp_out) 
{
char in_line[256], in_line2[256], in_line3[256]; 
char object_name[40], ads_name[40], discov_name[40], comp_name[40]; 
char buffer[256], companion[64], *pc;
int ads_nber;
int iline, norbits_found, status, verbose_if_error = 0;
int ref_slength = 60, nmax_orbits = 50;
register int i;

iline = 0;
while(!feof(fp_calib)) {
  if(fgets(in_line, 256, fp_calib)) {
    iline++;
// Good lines start with a digit (WDS names...)
// Lines starting with % are ignored
    if(in_line[0] != '%') { 
// Remove the end of line '\n' from input line:
      cleanup_string(in_line, 256);

// Remove orbit column:
      remove_orbit_column(in_line, 256);
// Remove 11th column:
//      remove_column(in_line, 11, 256);

// Save to output file:
    fprintf(fp_out, "%s\n", in_line);
    } else {
// Simply copy other lines:
    fprintf(fp_out, "%s", in_line);
    }/* EOF if % ... */
  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("modif_calib_table1: %d lines sucessfully read and processed\n", 
        iline);
return(0);
}
/***********************************************************************
*
* Remove the i_col th column
*
* INPUT:
*  ilen : length of in_line (=input line)
*  i_col : column number (from 1...) 
***********************************************************************/
static int remove_column(char *in_line, int i_col, int ilen)
{
char *pc, buffer1[256], buffer2[256];
int icol;

if(ilen > 256) {
  fprintf(stderr, "Fatal error: ilen=%d > 256! \n", ilen);
  exit(-1);
 }

if(in_line[0] == '\%') return(0);

strcpy(buffer1, in_line);
buffer1[255] = '\0';

/* Go to the _col th column: */
pc = buffer1;
icol = 0;
while(*pc) {
while(*pc && (*pc != '&') && 
      strncmp(pc, "\\cr", 3) && strncmp(pc, "\\\\", 2)) pc++;
if(*pc == '&') {
  icol++;
// If column number is OK, stop scanning buffer1 here and copy remaining 
// of the line to buffer2 (starting with "&")
  if(icol == i_col) {
// Skip "&" and copy remaining to buffer2:
   pc++;
   sprintf(buffer2, "%s", pc );
   buffer2[255] = '\0';
// Cut buffer1 after "&":
   *pc = '\0';
   break;
   }
// End of line:
  } else {
  *pc = '\0';
  }
  pc++;
}

/* DEBUG:
printf("FIRST/remove_column: icol=%d\n in: %s \n buffer1: >%s< buffer2: >%s<\n", 
        i_col, in_line, buffer1, buffer2);
*/

// Cut last column of buffer1:
if(i_col == 1) {
strcpy(buffer1, "");
} else {
icol = 0;
pc = buffer1;
while(*pc) {
 while(*pc && (*pc != '&') && 
      strncmp(pc, "\\cr", 3) && strncmp(pc, "\\\\", 2)) pc++;
 if(*pc == '&') {
   icol++;
// If column number stop buffer1 here and copy remaining of the line to buffer2
   if(icol == i_col - 1) {
// Cut after the "&" in buffer1:
    pc++;
    *pc = '\0';
    break;
    }
// End of line:
   } else {
   *pc = '\0';
   }
 }
 pc++;
}

/* DEBUG:
printf("SECOND/remove_column: buffer1: >%s< buffer2: >%s< \n out: %s %s\n", 
        buffer1, buffer2, buffer1, buffer2);
*/

sprintf(in_line, "%s %s", buffer1, buffer2);

return(0);
}
/***********************************************************************
*
* Remove the 11th column corresponding to the orbit column
***********************************************************************/
static int remove_orbit_column(char *in_line, int ilen)
{
int verbose_if_error = 0;
char *pc, buffer[256], notes[40];

if(ilen > 256) {
  fprintf(stderr, "Fatal error: ilen=%d > 256! \n", ilen);
  exit(-1);
 }

/* Get notes from column 12: */
if(latex_get_column_item(in_line, notes, 12, verbose_if_error)){
  fprintf(stderr, "Fatal error reading the notes: bad syntax of line=%s \n", in_line);
  exit(-1);
 }
/* Reduce length if full of ' ' ... */
trim_string(notes, 40);

strcpy(buffer, in_line);
/* Go to the end of Latex line (marked with "\\" or "\cr"): */
pc = buffer;
buffer[255] = '\0';
while(*pc && strncmp(pc, "\\cr", 3) && strncmp(pc, "\\\\", 2)) pc++;
if(!*pc) {
  fprintf(stderr, "Fatal error: bad syntax of line=%s \n", in_line);
  exit(-1);
 }
pc--;
while((pc != buffer) && (*pc != '&')) pc--;
pc--;
while((pc != buffer) && (*pc != '&')) pc--;
*pc = '\0';

sprintf(in_line, "%s & %s \\\\", buffer, notes);

return(0);
}

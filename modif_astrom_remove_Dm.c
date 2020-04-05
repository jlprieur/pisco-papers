/************************************************************************
* "modif_astrom_remove_Dm.c"
*
* To modify the astrom_15b.txt for example 
* (removes Dm=0.05+\-0.02 at the end of a line)
*  Dm=0.02+/-0.02\\ becomes \\
*
* JLP 
* Version 23/05/2018
*************************************************************************/
#include "jlp_catalog_utils.h"

#define DEBUG
#define DEBUG_1
/*
*/

static int modif_astrom2(FILE *fp_resid, FILE *fp_out);

int main(int argc, char *argv[])
{
char resid_fname[80], out_fname[80];
FILE *fp_resid, *fp_out;
time_t t = time(NULL);

if(argc != 3) {
  printf("Syntax: modif_astrom1 old_resid_table out_resid_table\n");
  return(-1);
}
strcpy(resid_fname, argv[1]);
strcpy(out_fname, argv[2]);

printf("OK: in_resid=%s output=%s \n", resid_fname, out_fname); 

/* Open input residual table: */
if((fp_resid = fopen(resid_fname, "r")) == NULL) {
   fprintf(stderr, "modif_resid_table1/Fatal error opening resid table %s\n",
           resid_fname);
    return(-1);
  }

/* Open output table: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "modif_resid_table1/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the output file: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "%% Modified file from: %s \n%% Created on %s", 
        resid_fname, ctime(&t));
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

/* Scan the input calibrated table and add the residuals 
*/
 modif_astrom2(fp_resid, fp_out); 

/* Close opened files:
*/
fclose(fp_resid);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input table and make the modifications 
*
* INPUT:
* fp_resid: pointer to the input file containing the input table 
* fp_out: pointer to the output Latex file
*
*************************************************************************/
static int modif_astrom2(FILE *fp_resid, FILE *fp_out) 
{
char in_line[256], in_line2[256], in_line3[256]; 
char buffer[256], my_string[64], *pc;
int ads_nber;
int iline, norbits_found, status, verbose_if_error = 0;
int ref_slength = 60, nmax_orbits = 50;
register int i;

strcpy(my_string, "Dm=");
iline = 0;
while(!feof(fp_resid)) {
  if(fgets(in_line, 256, fp_resid)) {
    iline++;
// Remove the end of line '\n' from input line:
    cleanup_string(in_line, 256);

  strcpy(in_line3, in_line);
  if(in_line3[0] == '&') {
      strcpy(buffer, in_line3);
// Search for substring "my_string" in string "in_line3":
      pc = strstr(in_line3, my_string); 
      if(pc != NULL) {
      *pc = '\\';
      pc++;
      *pc = '\\';
      pc++;
      *pc = '\0';
      }
   }

// Save to output file:
    fprintf(fp_out, "%s\n", in_line3);

  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("modif_resi_table1: %d lines sucessfully read and processed\n", 
        iline);
return(0);
}

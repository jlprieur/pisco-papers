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
#include "jlp_string.h"  // (in jlplib/jlp_fits/ ) jlp_cleanup_..

#define DEBUG
#define DEBUG_1
/*
*/

static int modif_astrom_removeDm(FILE *fp_in1, FILE *fp_out);

int main(int argc, char *argv[])
{
char in_file1[80], out_fname[80];
FILE *fp_in1, *fp_out;
time_t t0 = time(NULL);

if(argc != 3) {
  printf("Syntax: modif_astrom_remove_Dm in_astrom_file out_astrom_file\n");
  return(-1);
}
strcpy(in_file1, argv[1]);
strcpy(out_fname, argv[2]);

printf("OK: in_resid=%s output=%s \n", in_file1, out_fname); 

/* Open input astrom file: */
if((fp_in1 = fopen(in_file1, "r")) == NULL) {
   fprintf(stderr, "modif_resid_table1/Fatal error opening input file %s\n",
           in_file1);
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
        in_file1, ctime(&t0));
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

/* Scan the input calibrated table and add the residuals 
*/
 modif_astrom_removeDm(fp_in1, fp_out); 

/* Close opened files:
*/
fclose(fp_in1);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input table and make the modifications 
*
* INPUT:
* fp_in1: pointer to the input file containing the input table 
* fp_out: pointer to the output Latex file
*
*************************************************************************/
static int modif_astrom_removeDm(FILE *fp_in1, FILE *fp_out) 
{
char in_line[256], in_line2[256], in_line3[256]; 
char buffer[256], my_string[64], *pc;
int iline;

strcpy(my_string, "Dm=");
iline = 0;
while(!feof(fp_in1)) {
  if(fgets(in_line, 256, fp_in1)) {
    iline++;
// Remove the end of line '\n' from input line:
    jlp_cleanup_string(in_line, 256);

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
printf("modif_astrom_removeDm: %d lines sucessfully read and processed\n", 
        iline);
return(0);
}

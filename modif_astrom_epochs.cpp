/************************************************************************
* "modif_astrom_epochs.c"
*
* To modify the m1112_astrom.txt for example 
* correct the epochs
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

static int modif_astrom_epochs1(FILE *fp_in, FILE *fp_out);

int main(int argc, char *argv[])
{
char in_fname[80], out_fname[80];
FILE *fp_in, *fp_out;
time_t t = time(NULL);

if(argc != 3) {
  printf("Syntax: modif_astrom_epochs in_astrom_file out_astrom_file\n");
  return(-1);
}
strcpy(in_fname, argv[1]);
strcpy(out_fname, argv[2]);

printf("OK: input=%s output=%s \n", in_fname, out_fname); 

/* Open input residual table: */
if((fp_in = fopen(in_fname, "r")) == NULL) {
   fprintf(stderr, "modif_astrom_epochs/Fatal error opening resid table %s\n",
           in_fname);
    return(-1);
  }

/* Open output table: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "modif_astrom_epochs/Fatal error opening output file: %s\n",
           out_fname);
    return(-1);
   }

/* Header of the output file: */
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
fprintf(fp_out, "%% Modified file from: %s \n%% Created on %s", 
        in_fname, ctime(&t));
fprintf(fp_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");

/* Scan the input calibrated table and add the residuals 
*/
 modif_astrom_epochs1(fp_in, fp_out); 

/* Close opened files:
*/
fclose(fp_in);
fclose(fp_out);
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
static int modif_astrom_epochs1(FILE *fp_in, FILE *fp_out) 
{
char in_line[256], in_line3[256]; 
char buffer[256], epoch_string[64], title_string[64], *pc0, *pc, *pc1;
int iposition, iline, status, ival;
double new_epoch, old_epoch;
register int i;

old_epoch = 2011;
new_epoch = 2011;
strcpy(epoch_string, "EP=");
strcpy(title_string, "0 & & & & & & &");
iline = 0;
while(!feof(fp_in)) {
  if(fgets(in_line, 256, fp_in)) {
    iline++;
// Remove the end of line '\n' from input line:
    jlp_cleanup_string(in_line, 256);

  strcpy(in_line3, in_line);
  if(in_line3[0] == '&') {
      pc0 = strstr(in_line3, title_string); 
      if(pc0 != NULL) {
        printf("OK, title with null epoch found: >%s<\n", pc0);
        iposition = pc0 - in_line3;
        printf("iposition=%d\n", iposition);
        strcpy(buffer, in_line3);
        sprintf(&buffer[iposition], "%d & & & & & & & \\\\", (int)old_epoch);
// Copy back to in_line3
        strcpy(in_line3, buffer);
      } else {
// Search for substring "epoch_string" in string "in_line3":
      pc = strstr(in_line3, epoch_string); 
      if(pc != NULL) {
// Position of EP= :
        iposition = pc - in_line3;
// End of line:
        pc1 =  &in_line3[iposition + 3];
        while(*pc1 && (isdigit(*pc1) || *pc1 == '.')) pc1++;
// Read the epoch value:
        ival = sscanf(pc, "EP=%lf", &new_epoch);
        if(ival == 1) {
// Replace the epoch by the old epoch if epoch is dummy:
          if(new_epoch == 0.) {
            strcpy(buffer, in_line3);
            sprintf(&buffer[iposition], "EP=%9.4f", old_epoch);
            strcpy(&buffer[iposition + 12], pc1);
// Copy back to in_line3
            strcpy(in_line3, buffer);
// Load a new value of the old epoch if epoch is OK:
          } else {
            old_epoch = new_epoch;
          }
        } // if ival == 1
      } // if pc != NULL
      } // else if pc0 != NULL
   }

// Save to output file:
    fprintf(fp_out, "%s\n", in_line3);

  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("modif_resi_table1: %d lines sucessfully read and processed\n", 
        iline);
return(0);
}

/*************************************************************************
* astrom_add_epoch
* Program to add epoch (as recorded in FITS header of autocorrelation) 
* to astrom files (in the notes column of each measurement).
* For instance: EP=2009.2354 (as a fractional year)
*
* Format of input files:
*  Example of input:
& 120109\_ads684ab\_Rd\_8\_a & 12/01/2009 & R & 20 & 11.71 & 0.10 & -19.37 & 0.9 & Q=? \\
*
*  Example of output:
& 120109\_ads684ab\_Rd\_8\_a & 12/01/2009 & R & 20 & 11.71 & 0.10 & -19.37 & 0.9 & Q=? EP=2009.0852 \\
*
Possible keywords in the notes:
EP=2004.1323 (epoch)
Q=2         (Quadrant with restricted triple-correplation)
LQ=3        (Quadrant with long integration)
*
* JLP
* Version 19/04/2010
*************************************************************************/
#include "astrom_utils1.h" 
#include "astrom_utils2.h" 

int main(int argc, char *argv[])
{
char filein[60], fileout[60], fits_directory[100];
FILE *fp_in, *fp_out;

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
  printf(" Syntax: astrom_add_epoch in_file out_file fits_directory\n");
  printf(" Example: runs astrom_add_epoch astrom05a.tex astrom05aa.tex /home/data/pisco_merate/2004-2008/ \n");
  exit(-1);
  }
else
  {
  strcpy(filein,argv[1]);
  strcpy(fileout,argv[2]);
  strcpy(fits_directory,argv[3]);
  }

printf(" OK: filein=%s fileout=%s \n", filein, fileout);
printf(" OK: fits_directory=%s\n", fits_directory);

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
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
fprintf(fp_out,"%%%% Epoch added using FITS files located in %s\n", fits_directory);
fprintf(fp_out,"%%%% JLP / Version of 17/04/2008 \n");
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

/* Scan the file and add epoch from FITS autocorrelation files:
* (in "astrom_utils2.c")
*/ 
  astrom_add_epoch_from_fits_file(fp_in, fp_out, fits_directory); 

fclose(fp_in);
fclose(fp_out);
return(0);
}

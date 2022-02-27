/*************************************************************************
* astrom_add_WDS_pdb
* With Pisco data base
* Program to add WDS name from discoverer's name
* and the last reported observations in WDS, to astrom files.
* For instance: WY=2006 WT=12 WR=2.4 (Year, Theta (deg), Rho (arcsec))
*
* Format of input files:
\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}
& Name  & Epoch & Filter & Eyep. & $\rho$ & $\sigma_\rho$ & $\theta$ &
$\sigma_\theta$ & comments \\
 &       &       &        & (mm)
& (pixels) & (pixels) & (deg.) & (deg.) & \\
\hline % -----------------------------------------------------------------------
%%
*  Example of input:
& ADS 10279 & 2004. & & & & & & & \\
*  Example of output:
16564+6502 = STF 2118 AB & ADS 10279 & 2004. & & & & & & & orb 
WY=1999 WT=187 WR=1.4 \\
%%
*
Keywords in the notes:
EP=2004.133 (epoch)
WR=127      (Rho, WDS)
WT=127      (Theta, WDS)
WY=2002     (Year of WDS observation) 
Q=2         (Quadrant with restricted triple-correplation)
LQ=3        (Quadrant with long integration)
*
* From astrom_add_WDS / version 28/06/2009
*
* JLP
* Version 24/05/2018
*************************************************************************/
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include "astrom_utils_pdb.h" 

int main(int argc, char *argv[])
{
char filein[60], fileout[60], wds_cat[128];
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
  printf(" Syntax: astrom_add_WDS_pdb in_file out_file WDS_cat ADS_WDS_crossref\n");
  printf(" Example: runs astrom_add_WDS_pdb in_astrom15a.tex out_astrom15b.tex wdsweb_sum.txt\n");
  exit(-1);
  }
else
  {
  strcpy(filein, argv[1]);
  strcpy(fileout, argv[2]);
  strcpy(wds_cat, argv[3]);
  }

printf(" OK: filein=%s fileout=%s\n", filein, fileout);
printf(" OK: WDS_cat=%s\n", wds_cat);

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
fprintf(fp_out,"%%%% WDS names added using %s\n", wds_cat);
fprintf(fp_out,"%%%% JLP / Version of 19/09/2020 \n");
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

/* Scan the file and add various parameters from PISCO catalog:
* (in pscplib/astrom_utils_pdb.cpp)
*/ 
  astrom_add_WDS_from_discov(fp_in, fp_out, wds_cat); 

fclose(fp_in);
fclose(fp_out);
return(0);
}

/*************************************************************************
* astrom_add_WDS
* Program to add WDS name and discoverer's name
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
* JLP
* Version 28/06/2009
*************************************************************************/
#include "astrom_utils1.h" 
#include "astrom_utils2.h" 
#include "WDS_catalog_utils.h" 

int main(int argc, char *argv[])
{
char filein[60], fileout[60], pisco_cat[100], wds_cat[100];
FILE *fp_in, *fp_out;

/* If command line with "runs" */
if(argc == 7){
 if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
 else argc = 1;
 }


if(argc != 5)
  {
  printf(" Syntax: astrom_add_WDS in_file out_file PISCO_cat WDS_cat\n");
  printf(" Example: runs astrom_add_WDS astrom05a.tex astrom05aa.tex zeiss_doppie_new.cat wdsweb_sum.txt\n");
  exit(-1);
  }
else
  {
  strcpy(filein,argv[1]);
  strcpy(fileout,argv[2]);
  strcpy(pisco_cat,argv[3]);
  strcpy(wds_cat,argv[4]);
  }

printf(" OK: filein=%s fileout=%s\n", filein, fileout);
printf(" OK: pisco_cat=%s\n", pisco_cat);
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
fprintf(fp_out,"%%%% WDS names added using %s\n", pisco_cat);
fprintf(fp_out,"%%%% WDS data from %s\n", wds_cat);
fprintf(fp_out,"%%%% JLP / Version of 29/09/2008 \n");
fprintf(fp_out,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

/* Scan the file and add various parameters from PISCO catalog:
*/ 
  astrom_add_WDS(fp_in, fp_out, pisco_cat, wds_cat); 

fclose(fp_in);
fclose(fp_out);
return(0);
}

/************************************************************************
* "process_gili_table.c"
* To prepare the final LaTeX table from Gili's table 
*
* JLP 
* Version 22/09/2011
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>                  /* exit() */
#include <string.h>
#include <ctype.h>                   /* isprint... */
#include <math.h>
#include <time.h>                    /* date */
#include "jlp_trim.h"  /* (in jlplib/jlp_fits/ ) trim_string... */
#include "jlp_catalog_utils.h"       /* Routines used to read catalogs */
#include "WDS_catalog_utils.h"       /* Routines used to read WDS catalog */
#include "HIP_catalog_utils.h"       /* Routines used to read Hipparcos catalog */

#define DEGTORAD   (PI/180.00)

/*
#define DEBUG
#define DEBUG_1
*/

         
static int gili_table_process_line1(char *in_line, char *WDS_name, 
                                    int *has_an_orbit); 
static int write_header_gili_table(FILE *fp_out_cat);
static int write_end_gili_table(FILE *fp_out_cat);

int process_gili_table_main(char* in_gili_table, char *WDS_catalog,
                            char *HIC_catalog, char *HIP_catalog,
                            char *out_gili_table, double paral_rel_error);
int scan_and_process_gili_table(FILE *fp_in_cat, char *WDS_catalog, 
                                char *HIC_catalog, char *HIP_catalog,
                                FILE *fp_out_cat, double paral_rel_error);
static int process_gili_line(char *in_line, FILE *fp_out_cat,
                             char *WDS_name, char *WDS_catalog,
                             char *HIC_catalog, 
                             char *HIP_catalog, int *n_objects, int *n_HIP,
                             FILE *fp_out_data, double paral_rel_error);

int main(int argc, char *argv[])
{
char in_gili_table[80], out_gili_table[80]; 
char WDS_catalog[80], HIC_catalog[80], HIP_catalog[80];
double paral_rel_error;
int status, i;

if(argc != 7) {
  printf("Error: argc = %d\n", argc);
  for(i = 0; i < argc; i++) printf("argv[%d] = %s\n", i, argv[i]);
  printf("Syntax: process_gili_table out_table in_table WDS_catalog HIC_catalog HIP_main_catalog relative_error_on_parallax\n");
  printf("Example: process_gili_table tab.tex measures.tex wdsweb_summ.txt bds2ads2wds.txt wds2hds2hip.txt Hipparcos_HIC.dat Hipparcos_main.dat -1\n");
  printf("Enter a positive relative error to output absolute magnitudes of objects measured by Hipparcos\n");
  return(-1);
}
strcpy(out_gili_table, argv[1]);
strcpy(in_gili_table, argv[2]);
strcpy(WDS_catalog, argv[3]);
strcpy(HIC_catalog, argv[4]);
strcpy(HIP_catalog, argv[5]);
sscanf(argv[6], "%lf", &paral_rel_error);

printf("OK: in_cat=%s out_cat=%s\n", in_gili_table, out_gili_table);
printf("WDS=%s\n \n HIC=%s\n HIP=%s\n", WDS_catalog, HIC_catalog, HIP_catalog);
printf("paral_rel_error=%.3f\n", paral_rel_error);

/* Call process_gili_table_main that does the main job: */
status = process_gili_table_main(in_gili_table, WDS_catalog, 
                                 HIC_catalog, HIP_catalog, out_gili_table,
                                 paral_rel_error); 

return(status);
}
/************************************************************************
* process_gili_table_main
* main routine
*
* INPUT:
* in_gili_table: name of the input Latex catalog containing the list 
*             of PISCO targets 
* WDS_catalog: name of the WDS catalog
* HIC_catalog: name of the Hipparcos Input Catalog
* HIP_catalog: name of the Hipparcos-Tycho main catalog
* out_gili_table: name of modified catalog with the list of PISCO targets
*              and the WDS names 
*
*************************************************************************/
int process_gili_table_main(char* in_gili_table, char *WDS_catalog,
                            char *HIC_catalog, char *HIP_catalog, 
                            char *out_gili_table, double paral_rel_error)
{
FILE *fp_in_cat, *fp_out_cat; 
char buffer[64];
time_t tti = time(NULL);

/* Open input file containing the original Latex table */
if((fp_in_cat = fopen(in_gili_table, "r")) == NULL) {
   fprintf(stderr, "process_gili_table_main/Fatal error opening input table: %s\n",
           in_gili_table);
   return(-1);
  }

/* Open output file containing the modified Latex table */
if((fp_out_cat = fopen(out_gili_table, "w")) == NULL) {
   fprintf(stderr, "process_gili_table_main/Fatal error opening output table: %s\n",
           out_gili_table);
   fclose(fp_in_cat);
   return(-1);
  }

/* Header of the output Latex table: */
sprintf(buffer, ctime(&tti));
trim_string(buffer, 64);
fprintf(fp_out_cat, "%% From %s, modified by process_gili_table.c \n%% (%s)\n", 
        in_gili_table, buffer);

/* Scan the input file and compute the process_gili_table 
*/
scan_and_process_gili_table(fp_in_cat, WDS_catalog, HIC_catalog, HIP_catalog, 
                            fp_out_cat,paral_rel_error);

/* Close opened files:
*/
fclose(fp_in_cat);
fclose(fp_out_cat);
return(0);
}
/*********************************************************************
Ouput format:

\begin{table*}
\caption{}
\tabcolsep=1mm
\begin{tabular}{llrrlrlrllrrl}
\hline
WDS & Name & $(m_V)A$ & $\\Delta m_V$ & BD number & Epoch
& $\theta$ & $\rho$ & \\Delta m_F & Filter & Notes 
& \multicolumn{3}{c}{Residuals} \\
& & & & & & & & & & & & $\theta_{O-C}$ & $\rho_{O-C}$ & Reference\\
\hline
*********************************************************************/
static int write_header_gili_table(FILE *fp_out_cat)
{
fprintf(fp_out_cat, "\\afterpage{\\clearpage \n");
fprintf(fp_out_cat, "\\begin{table*} \n");
fprintf(fp_out_cat, "\\caption{Measurements of binaries with the Nice 76cm-refractor}\n");
fprintf(fp_out_cat, "\\tabcolsep=1mm \n");
fprintf(fp_out_cat, "\\begin{tabular}{llrrlrlrllrrl} \n");
fprintf(fp_out_cat, "\\hline \n\
WDS & Name & $m_{VA}$ & $\\Delta m_V$ & Epoch \n\
& $\\theta$ & $\\rho$ & $\\Delta m_F$ & Filter & Notes \n\
& $\\theta_{O-C}$ & $\\rho_{O-C}$ & Orbit\\\\\n\
& & mag & mag & & (\\degr) & (\\arcsec) & mag & \n\
& & (\\degr) & (\\arcsec) & \\\\\n");
fprintf(fp_out_cat, "\\hline \n");
return(0);
}
static int write_end_gili_table(FILE *fp_out_cat)
{
fprintf(fp_out_cat, "\\hline \n");
fprintf(fp_out_cat, "\\end{tabular} \n");
fprintf(fp_out_cat, "\\end{table*} } %%EOF afterpage \n \n");
fprintf(fp_out_cat, "\\addtocounter{table}{-1} \n");

return(0);
}
/***************************************************************************
* Scan the input table and format it for the output paper 
*
* INPUT:
* fp_in_cat: pointer to the input Latex table
* WDS_catalog: name of the input WDS catalog
* HIC_catalog: name of the Hipparcos Input catalog
* HIP_catalog: name of the Hipparcos-Tycho main catalog
* fp_out_cat: pointer to the output Latex table
***************************************************************************/
int scan_and_process_gili_table(FILE *fp_in_cat, char *WDS_catalog, 
                                char *HIC_catalog, char *HIP_catalog,
                                FILE *fp_out_cat, double paral_rel_error)
{
char WDS_name[40], outname[100];
char in_line[1024], buffer[200], *pc, *pc1;
int status, iline, n_objects, n_HIP, n_orbitals;
int has_an_orbit; 
register int i;
time_t ttime = time(NULL);
FILE *fp_out_data;

/* Open output file containing photometric data or HR 
* (depending on the sign of the error on the parallax,
*  entered in the command line) */
if(paral_rel_error < 0) {
 strcpy(outname, "gili_data_photom.dat");   /* NOT YET IMPLEMENTED ! */
 } else { 
 strcpy(outname, "gili_data_HR.dat");
}

if((fp_out_data = fopen(outname, "w")) == NULL){
   fprintf(stderr, "Fatal error opening output file (%s)\n", outname);
   exit(-1);
   }

/* Header of the output photometric file: */
if(paral_rel_error < 0) {
  fprintf(fp_out_data, "%% Photometric data created on %s",
         ctime(&ttime));
  fprintf(fp_out_data, "%% magV magV_A B-V Delta_mag discov_name\n");
/* Header of the output HR file: */
  } else {
  fprintf(fp_out_data, "%% HR diagram with paral_rel_error=%.3f\n%% Created on %s",
        paral_rel_error, ctime(&ttime));
  fprintf(fp_out_data, "%% magV B-V MagV WDS_name\n");
  }

/* Header of the output table file: */
write_header_gili_table(fp_out_cat);

/* Scan the input Latex table: */
iline = 0;
n_objects = 0;
n_HIP = 0;
n_orbitals = 0;
/* DEBUG: */
#ifdef DEBUG
printf("DEBUG/WARNING: limitation to the first 50 lines! \n");
while(!feof(fp_in_cat) && iline < 50) {
#else
while(!feof(fp_in_cat)) {
#endif
  pc = in_line; 
  if(fgets(buffer, 160, fp_in_cat)) {
    buffer[160] = '\0';
/* Eliminates all non printable characters: */
    cleanup_string(buffer, 160);
    iline++;
    if(isdigit(buffer[0])) {
      for(i = 0; i < 161 && buffer[i] != '\0' && buffer[0] != '\n'
                 && buffer[i] != '\r'; i++) *(pc++) = buffer[i];
      pc1 = buffer;
/* Look for end of line */
      while(*pc1){
        if(*pc1 == '\\') {
          pc1++;
          if(*pc1 == '\\') {
          break;
          }
        }
        pc1++;
      }
      pc1 = '\0';

         status = gili_table_process_line1(in_line, WDS_name, &has_an_orbit); 
         if(status) {
           fprintf(stderr, "scan_and_process/Fatal error/Bad object syntax in line #%d\n", iline);
           exit(-1);
           }
         if(has_an_orbit) n_orbitals++;
/* Before:
    fprintf(fp_out_cat, "%s\n", in_line);

Now this is done in "process_gili_line" with Hipparcos information at the end
*/
    process_gili_line(in_line, fp_out_cat, WDS_name, WDS_catalog, 
                      HIC_catalog, HIP_catalog, &n_objects, &n_HIP, 
                      fp_out_data, paral_rel_error);
/* Write header to cut the big table into pages of 60 lines: */
       if((iline != 0) && ((iline % 60) == 0)) {
         write_end_gili_table(fp_out_cat);
         write_header_gili_table(fp_out_cat);
       }
/* If commented line, simply copy this line to the ouput file: */
    } else {
      cleanup_string(buffer, 180);
      fprintf(fp_out_cat, "%s\n", buffer);
    }
 }
}

printf("End of process: %d lines read (%d objects modified, %d orbitals)\n", 
        iline, n_objects, n_orbitals);

printf("%d Hipparcos objects identified among the known binaries\n", 
        n_HIP );

write_end_gili_table(fp_out_cat);

fclose(fp_out_data);

return(0);
}
/***************************************************************************
*
21445+3933&  A  1447       &10.3 &10.3 & +38 4595 & 2008.689  & 241.1 &  0.337: 0.2     & IRC  &  LY                                \\\\  ncol=10
21439+2751&  HO  166       & 8.36& 8.25& +27 4145 & 2008.563  & 327.0 &  0.183 &i&      & IRC  & OR         &   -16.0& 0.001&Cou1958d \\  ncol=15
21435+4448&  HO  167       &10.36&11.52& +44 3916 & 2008.563  & 200.7 &  1.130 &i&      & IRC  &            &                      \\     ncol=13
***************************************************************************/
static int gili_table_process_line1(char *in_line, char *WDS_name, 
                                    int *has_an_orbit)
{
int icol, ncol, k, is_speckle, ival1, ival2;
char *pc, *pc1, items[32][64];
double theta, rho, mva, mvb;

/* Decode the line and save items into item array */
*has_an_orbit = 0;
pc = in_line;
icol = 0;
k = 0;
items[0][0] = '\0';
while(*pc) {
 if(*pc == '&' || *pc == '\\') {
   items[icol][k++] = '\0';
   icol++;
   k = 0;
   if(*pc == '\\') break;
   }
 else items[icol][k++] = *pc;
 pc++;
 }
ncol = icol;

/* Save WDS name (located in the first column) */
WDS_name[0]= '\0';
strcpy(WDS_name, items[0]);

/* Transform 4th column (mv_B) into Delta m_V */
ival1 = sscanf(&items[2][0],"%lf", &mva);
if(ival1 == 1) { 
   sprintf(&items[2][0], "%.1f", mva);
   ival2 = sscanf(&items[3][0],"%lf", &mvb);
   if(ival2 == 1) sprintf(&items[3][0], "%.1f", ABS(mva - mvb));
   else strcpy(&items[3][0], "\\nodata");
  } else {
   strcpy(&items[2][0], "\\nodata");
   strcpy(&items[3][0], "\\nodata");
  }

/* Removes heading and trailing blanks: */
for(icol = 0; icol < ncol; icol++) trim_string(&items[icol][0], 64);

if(ncol == 15 && !strncmp(&items[11][0], "OR", 2)) *has_an_orbit = 1;
 else *has_an_orbit = 0;

is_speckle = 0;
if((ncol == 15  || ncol == 13) && !strcmp(&items[8][0], "i")) is_speckle  = 1;

#ifdef DEBUG_1
printf(">%s< ncol=%d WDS=%s orbit=%d speckle=%d\n", in_line, ncol, 
       WDS_name, *has_an_orbit, is_speckle);
#endif

/* Error with orbitals: Delta mag is in the same column as $\rho$ */
if(ncol == 10) {
  ncol++;
  strcpy(&items[10][0], &items[9][0]);
  strcpy(&items[9][0], &items[8][0]);
  strcpy(&items[8][0], &items[7][0]);
  pc = &items[7][0];
  while(*pc && *pc != ' ') pc++;
  *pc = '\0';
  pc = &items[8][0];
  while(*pc && *pc != ' ') pc++;
  if(*pc == ' ') strcpy(&items[8][0], ++pc);
  else strcpy(&items[8][0], "");
/* Eliminates the last observation: */
  pc = &items[10][0];
  while(*pc && *pc != ' ') pc++;
  *pc = '\0';
}
/* Eliminates OR or OR- for orbitals: */
if(ncol == 15) {
  pc = &items[11][0];
  while(*pc) { 
  if(!strncmp(pc, "OR-", 3)) {*pc++ = ' '; *pc++ = ' '; *pc++ = ' '; break;}
  else if(!strncmp(pc, "OR", 2)) {*pc++ = ' '; *pc++ = ' '; break;}
  pc++;
  }
  trim_string(&items[11][0], 64);
/* Add ^Q if the orbits quadrant is not good */
  if(sscanf(&items[12][0],"%lf", &theta) == 1) {
   if(theta > 180.) theta -= 360.;
   if(theta < -180.) theta += 360.;
   if(ABS(theta - 180.) < 90.) {
      theta -= 180.;
      sprintf(&items[12][0],"$%.1f$\\rlap{$^Q$}", theta);
   } else if(ABS(theta + 180.) < 90.) {
      theta += 180.;
      sprintf(&items[12][0],"$%.1f$\\rlap{$^Q$}", theta);
/* Necessary for the minus sign: */
   } else {
     if(theta > 0)
      sprintf(&items[12][0],"$+%.1f$", theta);
     else
      sprintf(&items[12][0],"$%.1f$", theta);
   }
  }
/* Necessary for the minus sign: */
  if(sscanf(&items[13][0],"%lf", &rho) == 1) {
      if(rho > 0)
        sprintf(&items[13][0],"$+%.3f$", rho);
      else
        sprintf(&items[13][0],"$%.3f$", rho);
  }
}

/* Compact discoverer's name: */
  compact_string(&items[1][0], 64);

#ifdef DEBUG_1
for(icol = 0; icol < ncol; icol++) printf("#%d >%s<\n", icol, &items[icol][0]);
#endif

/* If not speckle interferometry, add "(lym)" in column #10 */
/* ns = non speckle measurement */
/* lym = Lucky imaging measurement */
  if(!is_speckle) sprintf(&items[10][0], "%s (lym) ", &items[10][0]);

/* If Lucky imaging witout any orbit add 3 empty columns 
* (theta_O-C, rho_O-C and reference of residuals) */
if(ncol == 11) {
  pc = in_line;
  for(icol = 0; icol < ncol; icol++) { 
/* Removes column 4 (BD number) */
    if(icol != 4) {
      pc1 = &items[icol][0];
      while(*pc1) *pc++ = *pc1++;
      if(icol != ncol - 1) {
        *pc++ = ' '; *pc++ = '&'; *pc++ = ' ';
        }
      }
    }
/* Add 3 empty columns (theta_O-C, rho_O-C and reference of residuals) */
   strcpy(pc, " & & & \\\\");
/*printf("DDD>%s<\n", in_line);
*/
} else if (ncol == 13) {
  pc = in_line;
  for(icol = 0; icol < ncol - 1; icol++) { 
/* Removes column 4 (BD number), column 8 (speckle) and last column (last observation from WDS) */
    if(icol != 4 && icol != 8) {
      pc1 = &items[icol][0];
      while(*pc1) *pc++ = *pc1++;
      if(icol != ncol - 2) {
        *pc++ = ' '; *pc++ = '&'; *pc++ = ' ';
        }
      }
    }
/* Add 3 empty columns (theta_O-C, rho_O-C and reference of residuals) */
   strcpy(pc, " & & & \\\\");
} else if (ncol == 15) {
  pc = in_line;
  for(icol = 0; icol < ncol; icol++) { 
/* Removes column 4 (BD number) and column 8 (speckle) */
    if(icol != 4 && icol != 8) {
      pc1 = &items[icol][0];
      while(*pc1) *pc++ = *pc1++;
      if(icol != ncol -1) {
        *pc++ = ' '; *pc++ = '&'; *pc++ = ' ';
        }
      }
    }
   strcpy(pc, " \\\\");
}

#ifdef DEBUG_1
printf(">%s<\n", in_line);
#endif
return(0);
}
/*************************************************************************
* Gili's table 
* Process the current line relative to an object
* and extract information from catalogues
* WDS0000+000 HIP00 PI_HIP=3.5(+/-2)mas \cr
*
* INPUT:
* in_line: current line of the input table to be processed
* fp_out_cat: pointer to the ouput table 
* WDS_name: WDS name of the object as read in gili's input table
*
* INPUT/OUTPUT
* n_objects: counter of the objects successfuly processed
* n_HIP: counter of the Hipparcos objects
*
*************************************************************************/
static int process_gili_line(char *in_line, FILE *fp_out_cat,
                             char *WDS_name, char *WDS_catalog, 
                             char *HIC_catalog, 
                             char *HIP_catalog, int *n_objects, int *n_HIP,
                             FILE *fp_out_data, double paral_rel_error)
{
char HIP_name[40], CCDM_name[40]; 
char *pc, buffer[128];
double alpha, delta, equinox, D_tolerance;
double paral, err_paral, V_mag, B_V_index, V_mag_abs;
int iline, status, has_an_orbit, found, found_in_Hip_cat;

/* Look for accurate coordinates in WDS catalog
* alpha, delta, equinox: coordinates of the object
*            (alpha in hours and delta in degrees)
*/
  read_coordinates_from_WDS_catalog(WDS_name, WDS_catalog, &alpha, &delta,
                                    &equinox, &found);
 if(!found) {
   fprintf(stderr,"\n process_gili_line/Warning!\n");
   fprintf(stderr, "read_coordinates/%s was not found in WDS catalog\n", 
           WDS_name);
   alpha = 0.; delta = 0.;
 }

/* Tolerance of the coordinates in degrees 
* used for searching for Hipparcos names.
* 0.1 arcminute is a good value: */
D_tolerance = 0.1/60.;

/* DEBUG: search in Hipparcos catalog */
 found_in_Hip_cat = 0;
/* No coordinates for some objects in WDS catalog, so I skip
* this part in that case: */
 if(alpha != 0. || delta != 0.) {
   search_object_in_HIC_catalog(HIC_catalog, alpha, delta, equinox,
                                HIP_name, CCDM_name, &V_mag, &B_V_index, 
                                D_tolerance, &found_in_Hip_cat);
   }
 if(found_in_Hip_cat) {
#ifdef DEBUG
   printf("%s was found in Hipparcos catalog (=HIP%s=CCDM%s) alpha=%f delta=%f\n", 
           WDS_name, HIP_name, CCDM_name, alpha, delta);
   printf("V=%.3f B-V=%.3f\n", V_mag, B_V_index);
#endif
   (*n_HIP)++;
   read_data_in_HIP_catalog(HIP_catalog, HIP_name, &paral, &err_paral, 
                            &found); 
/* (Parallax in mas) */
   if(found) {
#ifdef DEBUG
      printf("Paral=%.2f+/-%.2f\n", paral, err_paral);
#endif
/* compute M_V, B-V, and build a HR diagram if paral_rel_error > 0
*/
      if(paral > 0) {
      if(err_paral/paral < paral_rel_error) {
/* MV = mV - 2.5 log10 (d/d_10)^2
*  MV = mV - 5 log10 (d) + 5 log10(d_10)
*  MV = mV - 5 log10(d) + 5
* paral = 1/d        
* MV = mV + 5 log10(paral) + 5
*/
/* Paral in mas  so I multiply with 1/1000: */
        V_mag_abs = V_mag + 5. * log10(paral * 0.001) + 5.;
        compact_string(WDS_name, 40);
        fprintf(fp_out_data, "%.3f %.3f %.3f %s\n",
                V_mag, B_V_index, V_mag_abs, WDS_name);
      }
      } /* EOF case paral > 0 */
     } else {
      fprintf(stderr, "Update_2ndline/Fatal error: HIP=%s not in Hipparcos main catalog!\n", HIP_name);
      exit(-1);
     }
  }
#ifdef DEBUG_1
   else {
   printf("%s was not found in Hipparcos catalog (alpha=%f delta=%f)\n", 
           WDS_name, alpha, delta);
  } /* EOF if(!found_in_Hip_cat) */  
#endif


/* Accurate coordinates in last column of the new version of the WDS catalog:
000013.77+164056.6
*/

/* Copy the modified input line to ouput catalog: */
  strncpy(buffer, in_line, 128);
  pc = buffer;
  while(*pc && strncmp(pc, "\\cr", 3)) pc++;
  *pc = '\0';
/* Write input line and add Hipparcos information at the end of this line */
  fprintf(fp_out_cat, "%s", buffer);
  if(*HIP_name) {
    fprintf(fp_out_cat, "%%HIP%s", HIP_name);
    compact_string(CCDM_name, 20);
    if(*CCDM_name) fprintf(fp_out_cat, " CCDM%s", CCDM_name);
    if(V_mag != 100.) fprintf(fp_out_cat, " V=%.2f", V_mag);
    if(B_V_index != 100.) fprintf(fp_out_cat, " B-V=%.2f", B_V_index);
/* Parallax in mas: */
    if(paral > 0.) {
      if(err_paral > 0.) 
        fprintf(fp_out_cat, " PI=%.2f(%.2f)", paral, err_paral);
      else
        fprintf(fp_out_cat, " PI=%d", (int)paral);
      }
   }
/* Write "eof line": */
  fprintf(fp_out_cat, "\n");
  (*n_objects)++;

return(0);
}

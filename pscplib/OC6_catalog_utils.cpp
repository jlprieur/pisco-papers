/************************************************************************
* "OC6_catalog_utils.c"
* Set of routines used to read OC6 catalog
*
* JLP 
* Version 05/06/2015
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>    // exit()
#include <string.h>    // strcpy()
#include "jlp_numeric.h" // JLP_QSORT, MINI, MAXI, PI, etc
#include "jlp_catalog_utils.h"
#include "jlp_string.h"

/* The prototypes of routines included here
* are defined in "OC6_catalog_utils.h":
*/
#include "OC6_catalog_utils.h"

#ifndef MAXI
#define MAXI(a,b) ((a) < (b)) ? (b) : (a)
#endif

/*
#define DEBUG 
#define DEBUG_1 
*/
static int decode_authors_in_reference(char *reference0, char *decoded_authors,
                                       int length0);
static int decode_year_in_short_reference(char *author0, char *decoded_year,
                                    int length0);
static int remove_year_in_long_reference(char *reference1, int length1);

/***************************************************************************
* get_orbit_from_OC6_list
* Read input line from a list of selected orbits in OC6 format 
* and retrieve the orbital parameters and the measurements
* in OC6 ("Sixth Orbital Catalog") format (i.e., iformat = 2 or -2)
*
* Input lines contain
*  - name of object and orbital parameters in OC6 format (if iformat = 2 or -2) 
*  - the measurements (if iformat = -2)
*    otherwise, (if format = 2) the program will retrieve the measurements 
*      from the calibrated table.
*
* Example:
*
* With iformat = 2:
*
**** Old format (2016): 
000000.00+000000.0 00093+7943 STF   2          102    431    760   6.68   6.89    540.      y    .         0.995  a   .      110.1       .     171.2        .      1887.5     y    .       0.715     .       333.7       .  
   2000      3 n Hei1997  abcdefghikj.png
092059.40+381117.9 09210+3811 STF1338AB       7307  80441  45858   6.72   7.08    444.27    y    .         1.624  a   .       33.4       .     177.4        .      1983.69    y    .       0.247     .        83.6       .          1999 2002-03   3.4   n Sca2002b wds09210+3811b.png
**** New format (2018): 
111810.90+313144.9 11182+3132 STF1523A        8119  98231  55203   4.3     .        1.832   y    .         0.057  a   .       94.9       .     263.5   *    .      1986.495   y    .       0.53      .       143.0       .     2000 1994           9.0 n n Msn1995  wds11182+3132s.png
111810.90+313144.9 11182+3132 STF1523A        8119  98231  55203   4.3     .        1.834   y    .         0.054  a   .       91.        .     318.    *    .      1992.290   y    .       0.61      .       324.        .                         9.0 n y Hei1996b wds11182+3132t.png
111810.90+313145.0 11182+3132 STF1523AB       8119  98231  55203   4.33   4.80     60.72    y    .         3.278  a   .       56.10      .      97.78       .      1816.73    y    .       0.3777    .       134.37      .          1827           6.0 n n HJ_1833b
111810.90+313145.0 11182+3132 STF1523AB       8119  98231  55203   4.33   4.80     60.4596  y    .         2.290  a   .       52.26      .      95.01       .      1816.95    y    .       0.404     .       129.68      .          1836 1836-46   6.0 n n Mad1837

*
*
* INPUT:
* in_line: line from the list extracted from the OC6 catalog,
*           corresponding to the object
* iline: number of the corresponding line in input file
* is_master_file: flag set to 1 if master file ("orb6.master", 
*                          to 0 if OC6 file ("orb6orbits.txt")
* 
* OUTPUT:
* object designation, orbital parameters and arrays with the measurements
*
***************************************************************************/
int get_orbit_from_OC6_list(char *in_line, int iline, int is_master_file, 
                       char *WDS_name, char *ADS_name, char *discov_name, 
                       char *comp_name, char *object_name, char *author, 
                       double *Omega_node, double *omega_peri, double *i_incl, 
                       double *e_eccent, double *T_periastron, double *Period, 
                       double *a_smaxis, double *mean_motion, 
                       double *orbit_equinox, int *orbit_grade)
{
char buffer[80];
int nval;
double ww;

get_orbit_from_OC6_list_gili(in_line, iline, is_master_file, WDS_name, 
                             discov_name, comp_name, object_name, author, 
                             Omega_node, omega_peri, i_incl, e_eccent, 
                             T_periastron, Period, a_smaxis, mean_motion, 
                             orbit_equinox, orbit_grade);

strncpy(ADS_name, &in_line[45], 5);
ADS_name[5] = '\0';
/* Remove all blanks: */
jlp_compact_string(ADS_name, 40);

/* Restriction of the object name to ADS name or discoverer name */
if(ADS_name[0] != '\0' && ADS_name[0] != '.') 
      sprintf(object_name, "ADS %s", ADS_name);
else if (discov_name[0] != '\0') strcpy(object_name, discov_name);

#ifdef DEBUG_1
printf("DEBUG/Object=%s WDS=%s ADS=%s discov=%s comp=%s author=%s\n",
        object_name, WDS_name, ADS_name, discov_name, comp_name, author);
#endif

return(0);
}
/***************************************************************************
* get_orbit_from_OC6_list_gili
* Read input line from a list of selected orbits in OC6 format 
* and retrieve the orbital parameters and the measurements
* in OC6 ("Sixth Orbital Catalog") format (i.e., iformat = 2 or -2)
*
* Input lines contain
*  - name of object and orbital parameters in OC6 format (if iformat = 2 or -2) 
*  - the measurements (if iformat = -2)
*    otherwise, (if format = 2) the program will retrieve the measurements 
*      from the calibrated table.
*
* Example:
*
* With iformat = 2:
*
**** Old format (2016): 
000000.00+000000.0 00093+7943 STF   2          102    431    760   6.68   6.89    540.      y    .         0.995  a   .      110.1       .     171.2        .      1887.5     y    .       0.715     .       333.7       .  
   2000      3 n Hei1997  abcdefghikj.png
092059.40+381117.9 09210+3811 STF1338AB       7307  80441  45858   6.72   7.08    444.27    y    .         1.624  a   .       33.4       .     177.4        .      1983.69    y    .       0.247     .        83.6       .          1999 2002-03   3.4   n Sca2002b wds09210+3811b.png
**** New format (2018): 
111810.90+313144.9 11182+3132 STF1523A        8119  98231  55203   4.3     .        1.832   y    .         0.057  a   .       94.9       .     263.5   *    .      1986.495   y    .       0.53      .       143.0       .     2000 1994           9.0 n n Msn1995  wds11182+3132s.png
111810.90+313144.9 11182+3132 STF1523A        8119  98231  55203   4.3     .        1.834   y    .         0.054  a   .       91.        .     318.    *    .      1992.290   y    .       0.61      .       324.        .                         9.0 n y Hei1996b wds11182+3132t.png
111810.90+313145.0 11182+3132 STF1523AB       8119  98231  55203   4.33   4.80     60.72    y    .         3.278  a   .       56.10      .      97.78       .      1816.73    y    .       0.3777    .       134.37      .          1827           6.0 n n HJ_1833b
111810.90+313145.0 11182+3132 STF1523AB       8119  98231  55203   4.33   4.80     60.4596  y    .         2.290  a   .       52.26      .      95.01       .      1816.95    y    .       0.404     .       129.68      .          1836 1836-46   6.0 n n Mad1837

* http://www.astro.gsu.edu/wds/orb6/orb6format.txt
*
  column  format           description
    1    T1,2I2,F5.2,     epoch-2000 right ascension (hours, minutes, seconds).
         A1,2I2,f4.1      epoch-2000 declination (degrees, minutes, seconds).
    2    T20,A10          WDS designation (based on arcminute-accuracy epoch-
                          2000 coordinates).
    3    T31,A14          Discover designation and components, or other catalog 
                          designation.
    4    T46,I5           ADS (Aitken Double Star catalog) number.
    5    T52,I6           HD catalog number.
    6    T59,I6           Hipparcos catalog number.
    7    T67,F5.2,A1      Magnitude of the primary (usually V), and flag:
                            > = fainter than quoted magnitude
                            < = brighter than quoted magnitude 
                            v = variable magnitude
                            k = magnitude is in K-band or other infrared band
                            ? = magnitude is uncertain
    8    T74,F5.2,A1      Magnitude of the secondary (usually V), and flag:
                            > = fainter than quoted magnitude
                            < = brighter than quoted magnitude 
                            v = variable magnitude
                            k = magnitude is in K-band or other infrared band
                            ? = magnitude is uncertain
    9    T82,F11.6,A1     Period (P) and code for units:
                            m = minutes (not yet used!)
                            h = hours (not yet used!)
                            d = days
                            y = years
                            c = centuries (rarely used)
....

    25    T234,I1          Orbit grade, ranging from 1 ("definitive") to 5 
                          ("indeterminate"). Additionally, a grade of 8 is used
                          for interferometric orbits based on visibilities 
                          rather than rho and theta measures (hence not gradable
                          by the present scheme) and a grade of 9 indicates an
                          astrometric binary (also lacking rho and theta data).

*
*
* INPUT:
* in_line: line from the list extracted from the OC6 catalog,
*           corresponding to the object
* iline: number of the corresponding line in input file
* is_master_file: flag set to 1 if master file ("orb6.master", 
*                          to 0 if OC6 file ("orb6orbits.txt")
* 
* OUTPUT:
* object designation, orbital parameters and arrays with the measurements
*
***************************************************************************/
int get_orbit_from_OC6_list_gili(char *in_line, int iline, int is_master_file, 
                       char *WDS_name, char *discov_name, 
                       char *comp_name, char *object_name, char *author, 
                       double *Omega_node, double *omega_peri, double *i_incl, 
                       double *e_eccent, double *T_periastron, double *Period, 
                       double *a_smaxis, double *mean_motion, 
                       double *orbit_equinox, int *orbit_grade)
{
char buffer[80];
int nval;
double ww;

#ifdef DEBUG
printf("get_orbit_from_OC6_list_gili/iline=%d\n %s", iline, in_line);
#endif
/* Retrieve the object designation and the author of the orbit
*/
strncpy(WDS_name, &in_line[19], 10);
WDS_name[10] = '\0';

#ifdef DEBUG_1
printf("Discov_name: %.7s comp=%.7s \n", &in_line[30], &in_line[37]);
printf("Period=%.11s unit=%c \n", &in_line[80], in_line[92]);
printf("semi-major axis: a=%.20s unit=>%c< \n", &in_line[105], in_line[114]);
printf("Inclination i=%.8s (deg) \n", &in_line[125]);
printf("Omega Node=%.8s (deg) \n", &in_line[143]);
printf("Longitude of periastron, omega=%.8s (deg) \n", &in_line[205]);
printf("Orbit grade=%.2s \n", &in_line[233]);
#endif

/* Discover's name in fixed format with 7 characters (like with the WDS) */
strncpy(discov_name, &in_line[30], 7);
discov_name[7] = '\0';
strcpy(object_name, discov_name);

strncpy(comp_name, &in_line[37], 7);
comp_name[7] = '\0';
/* Remove all blanks: */
jlp_compact_string(comp_name, 10);

/* Period: */
strncpy(buffer, &in_line[80], 11);
buffer[11] = '\0';
if(sscanf(buffer, "%lf", Period) != 1) {
  fprintf(stderr, "\nWDS=%s Error reading period: buffer=%s\n", WDS_name, buffer);
  return(-1);
  }

/* Units: minutes, hours, days, centuries, or years */
if(in_line[92] == 'm') 
  *Period /= (365.25 * 24. * 60.);
else if(in_line[92] == 'h') 
  *Period /= (365.25 * 24.);
else if(in_line[92] == 'd') 
  *Period /= 365.25;
else if(in_line[92] == 'c') 
  *Period *= 100.;
else if(in_line[92] != 'y') 
  {
  fprintf(stderr, "\nError reading period of %s: unknown unit (P=%.20s unit=%c)\n", 
          WDS_name, &in_line[81], in_line[92]);
  return(-1);
  }

/* Semi-major axis: */
strncpy(buffer, &in_line[105], 9);
buffer[9] = '\0';
*a_smaxis = 0.;
if((nval = sscanf(buffer, "%lf", a_smaxis)) != 1) {
  fprintf(stderr, "\nError reading semi-major axis: nval=%d buffer=%s\n", 
          nval, buffer);
  fprintf(stderr, "iline=%d a_smaxis=%f inline=%s\n", iline, *a_smaxis, in_line);
  return(-1);
  }
/* Units: milliarcseconds or arcseconds */
if(in_line[114] == 'm') 
  *a_smaxis /= 1000.;
else if(in_line[114] != 'a') 
  {
  fprintf(stderr, "\nError reading semi-major axis: unknown unit (a=%.20s unit=%c)\n", 
          &in_line[105], in_line[114]);
  return(-1);
  }

/* Inclination (degrees) */
strncpy(buffer, &in_line[125], 8);
buffer[8] = '\0';
if(sscanf(buffer, "%lf", i_incl) != 1) {
  fprintf(stderr, "\nError reading inclination: buffer=%s\n", buffer);
  return(-1);
  }

/* Node, Omega (degrees) */
strncpy(buffer, &in_line[143], 8);
buffer[8] = '\0';
if(sscanf(buffer, "%lf", Omega_node) != 1) {
  fprintf(stderr, "\nError reading Omega (node): buffer=%s\n", buffer);
  return(-1);
  }

/* T of periastron passage */
strncpy(buffer, &in_line[162], 12);
buffer[12] = '\0';
if(sscanf(buffer, "%lf", T_periastron) != 1) {
  fprintf(stderr, "WDS_name=%s /error reading T periastron: buffer=%s\n", 
          WDS_name, buffer);
  return(-1);
  }
/* Units: modified Julian date or fractionnal Besselian year */
/* The time of periastron passage (T0) and code for units:
                            d = Julian date (-2,400,000 days)
                            m = modified Julian date (MJD = JD-2,400,000.5 days)
                            y = fractional Besselian year
*/

/* Orbit grade: */
strncpy(buffer, &in_line[233], 1);
buffer[1] = '\0';
*orbit_grade = -1;
if(sscanf(buffer, "%d", orbit_grade) != 1) {
/********** DEBUG:
  fprintf(stderr, "WDS_name=%s /error reading Orbit grade: buffer=%s\n", 
          WDS_name, buffer);
******/
  }

/*
printf("JLPPPDDEBUG: WDS_name=%s raw periastron = %f yr (unit=%c)\n", 
        WDS_name, *T_periastron, in_line[174]);
*/

/* Modified JD: JD - 2400000
*  reduced JD:  JD - 2400000.5
*/
if(in_line[174] == 'm')
  {
/* JLP 2015... */
  *T_periastron = 1900.0 + (*T_periastron - 15020.31352 - 0.5) / 365.242198781; 
  }
else if(in_line[174] == 'd')
  {
/* d = modified Julian date (JulianDate - 2,400,000 days) */
/* Bessel epoch = 1900.0 + (JulianDate - 2415020.31352) / 365.242198781 */
  *T_periastron = 1900.0 + (*T_periastron - 15020.31352) / 365.242198781; 
  }
else if(in_line[174] != 'y') 
  {
  fprintf(stderr, "\nError reading T_periastron: unknown unit (P=%.20s unit=%c)\n", 
          &in_line[162], in_line[174]);
  return(-1);
  }

/* Eccentricity */
strncpy(buffer, &in_line[187], 8);
buffer[8] = '\0';
if(sscanf(buffer, "%lf", e_eccent) != 1) {
  fprintf(stderr, "\nError in line: %s\n", in_line);
  fprintf(stderr, "Error reading eccentricity: buffer=%s\n", buffer);
  return(-1);
  }

/* omega (longitude of periastron) */
strncpy(buffer, &in_line[205], 8);
buffer[8] = '\0';
if(sscanf(buffer, "%lf", omega_peri) != 1) {
  fprintf(stderr, "\nError reading omega_periastron: buffer=%s\n", buffer);
  return(-1);
  }

/* Equinox */
strncpy(buffer, &in_line[223], 4);
buffer[4] = '\0';
/* Default value is 2000.0 */
if(sscanf(buffer, "%lf", orbit_equinox) != 1) {
  *orbit_equinox = 2000.0;
}

/* Reference of orbit (author) */
if(is_master_file) 
   strncpy(author, &in_line[251], 8);
  else
   strncpy(author, &in_line[237], 8);
author[8] = '\0';
/* Remove heading and trailing blanks: */
jlp_trim_string(author, 60);

/*
printf("JLPPPDDEBUG: WDS_name=%s periastron = %f yr (unit=%c, author=%s)\n", 
        WDS_name, *T_periastron, in_line[174], author);
*/

#ifdef DEBUG_1
printf("Object=%s WDS=%s discov=%s comp=%s author=%s\n",
        object_name, WDS_name, discov_name, comp_name, author);

 printf("Omega_node=%.3f omega_peri=%.3f incl=%.3f e=%.4f T=%.3f P=%.3f a=%.5f Equinox=%.3f\n", 
        *Omega_node, *omega_peri, *i_incl, *e_eccent, *T_periastron, 
        *Period, *a_smaxis, *orbit_equinox);
#endif

/* Conversion to radians: */
 *Omega_node *= DEGTORAD;
 *omega_peri *= DEGTORAD;
 *i_incl *= DEGTORAD;
 *mean_motion = (360.0 / *Period) * DEGTORAD;

return(0);
}
/***************************************************************************
* line_extraction_from_OC6_catalog ("orb6.master" in November 2009)
* Return the line in OC6 that corresponds to the input object name 
*
* INPUT:
* OC6_fname: name of master file or Sixth Orbit catalog 
* is_master_file: flag set to 1 if master file ("orb6.master", 
*                          to 0 if OC6 file ("orb6orbits.txt")
* ads_name[60]: ADS name of object ('\0' if not in ADS)
* discov_name[60]: discovery name of object 
* comp_name: name of the companion (e.g., AB, Aa, etc) 
* norbits_per_object: maximum number of orbits to be loaded to output list
*                     fror each object
* 
* OUTPUT:
* fp_out: pointer to the output ASCII file 
* found: flag set to one if at least one orbit was found for object
*
***************************************************************************/
int line_extraction_from_OC6_catalog(char *OC6_fname, int is_master_file,
                                     char *ads_name, 
                                     char *discov_name, char *comp_name, 
                                     FILE *fp_out, int *found, 
                                     int *candidate_found, 
                                     int norbits_per_object)
{
/* Line length is 278 + "\n" for OC6 catalog (orb6.master, november 2009)... */
/* Line length is 264 + "\n" for OC6 catalog (orb6orbits.txt, november 2009)... */
int iline, status, max_norbits = 1024, norbits, imin, nlines_in_header;
int OC6_comp_is_AB, comp_is_AB, discov_name_only;
char OC6_comp_really_compacted[40], comp_really_compacted[40];
char OC6_ads_name[60], OC6_comp_name[40], OC6_discov_name[40];
char line_buffer[300], compacted_ads_name[60], compacted_comp_name[40];
char compacted_discov_name[40];
/* Assume that the number of orbits for one object is always less than 1024: */
char object_orbits[300*1024];
FILE *fp_in;
int i;

*found = 0;
*candidate_found = 0;
norbits = 0;

strcpy(compacted_ads_name, ads_name);
jlp_compact_string(compacted_ads_name, 60);

if(compacted_ads_name[0] == '\0') discov_name_only = 1;
else discov_name_only = 0; 

strcpy(compacted_comp_name, comp_name);
jlp_compact_string(compacted_comp_name, 40);
jlp_really_compact_companion(compacted_comp_name, comp_really_compacted, 40);

strcpy(compacted_discov_name, discov_name);
jlp_compact_string(compacted_discov_name, 40);

#ifdef DEBUG
 printf("CURRENT OBJECT: ads_name=%s comp_name=%s discov_name=%s (discov_name_only=%d\n", 
ads_name, comp_name, discov_name, discov_name_only);
#endif

/* Open OC6 catalog: */
if((fp_in = fopen(OC6_fname,"r")) == NULL) {
  fprintf(stderr, "line_extraction_from_OC6_catalog/Fatal error opening %s\n",
          OC6_fname);
  exit(-1);
  }

iline = 0;

if(is_master_file)
  nlines_in_header = 4;
else
  nlines_in_header = 8;

while(!feof(fp_in)) {
/* Should read more than 
* 278 characters in master OC6 file
* 264 characters in non-master OC6 file
* in order to be sure to copy the complete line: */
  if(fgets(line_buffer, 280, fp_in)) {
    iline++;
    if(*found && !strncmp(line_buffer,"    ",4)) break; 
/* Skip the header and empty lines (that
* generally indicate the end of a given object):*/
    if(iline > nlines_in_header && strncmp(line_buffer,"    ",4)) {
/* Get the object name of each line */
     status = get_name_from_OC6_line(line_buffer, OC6_ads_name,
                                     OC6_discov_name, OC6_comp_name);
     if(status) {
     fprintf(stderr, "line_extraction_from_OC6_catalog/Error processing line #%d\n", iline); 
     return(-1);
     }
#ifdef DEBUG_1
if(iline < nlines_in_header+3) 
     printf("OK1found=%d iline=%d: >%s<\n ads=%s< disc=%s< comp=%s<\n", 
             *found, iline, line_buffer, OC6_ads_name, OC6_discov_name, 
             OC6_comp_name);
#endif
/* If "ads_name" was found, copy the full line containing 
* the orbital elements: */
     jlp_compact_string(OC6_ads_name, 60);
/* CASE 1 : no ADS name*/
    if(discov_name_only){
    if(!strcmp(OC6_discov_name, compacted_discov_name)) {
        for(i = 0; i < 300; i++) 
           object_orbits[i + 300*norbits] = line_buffer[i];
        *found = 1;
        norbits++;
        if(norbits >= max_norbits) {
            fprintf(stderr, "line_extraction_from_OC6_catalog/Fatal error: %d orbits found for a single object!\n", norbits);
            fprintf(stderr, "(OC6_discov_name=%s)\n", OC6_discov_name);
            exit(-1);
            }
      }
/* CASE 2 */
    } else if(!strcmp(compacted_ads_name, OC6_ads_name)
         && !strcmp(OC6_discov_name, compacted_discov_name)) {
/* Test on the companion names if present in the object name: 
*/
     jlp_compact_string(OC6_comp_name, 40);
     jlp_really_compact_companion(OC6_comp_name, OC6_comp_really_compacted, 40);
#ifdef DEBUG
printf("DEBUG/OC6_object=%s OC6_companion=%s| comp=%s| (really compacted: OC6 comp=%s comp=%s)\n", 
        OC6_ads_name, OC6_comp_name, comp_name, OC6_comp_really_compacted, 
        comp_really_compacted);
#endif
        comp_is_AB = 0;
        if((compacted_comp_name[0] == '\0') 
           || !strcmp(compacted_comp_name,"AB")
           || !strncmp(compacted_comp_name,"Aa-B",4)) comp_is_AB = 1;
        OC6_comp_is_AB = 0;
        if((OC6_comp_name[0] == '\0') || !strcmp(OC6_comp_name,"AB")
           || !strncmp(OC6_comp_name,"Aa-B",4)) OC6_comp_is_AB = 1;
           
        if((*compacted_comp_name != '\0' 
             && !strcmp(OC6_comp_name, compacted_comp_name))
/* If not mentionned in Latex calibrated table, should
* be either not mentioned in OC6 or equal to AB: */
           || (comp_is_AB && OC6_comp_is_AB)
           || (*comp_really_compacted != '\0' && 
               !strcmp(comp_really_compacted, OC6_comp_really_compacted))){
        for(i = 0; i < 300; i++) 
           object_orbits[i + 300*norbits] = line_buffer[i];
        *found = 1;
#ifdef DEBUG
        if(*found) {
          printf("From_OC6_cat/Object found now in OC6: >%s< >%s< >%s<\n", 
                  OC6_ads_name, OC6_discov_name, OC6_comp_name); 
          }
#endif
        norbits++;
        if(norbits >= max_norbits) {
            fprintf(stderr, "line_extraction_from_OC6_catalog/Fatal error: %d orbits found for a single object!\n", norbits);
            fprintf(stderr, "(OC6_ads_name=%s OC6_discov_name=%s)\n", 
                    OC6_ads_name, OC6_discov_name);
            exit(-1);
            }
        } else if(!(*found)) {
        printf("CURRENT OBJECT: ads_name=%s comp_name=%s discov_name=%s \n",
	 ads_name, comp_name, discov_name);
        printf("From_OC6_cat/Not yet found, possible candidate in OC6: >%s< >%s< >%s< (companion names look different though...)\n", 
               OC6_ads_name, OC6_discov_name, OC6_comp_name); 
         *candidate_found = 1;
        }
#ifdef DEBUG
  printf("DEBUG/iline=%d OC6_ads_name=>%s< OC6_comp_name=>%s<\n", 
          iline, OC6_ads_name, OC6_comp_name);
#endif
     } /* EOF if(!strcmp ads_name ...) */
    } /* EOF if line > nlines_header */
  } /* EOF if fgets */ 
 }

#ifdef DEBUG
printf("line_extraction_from_OC6_catalog: %d lines read and %d orbits found for current object\n", 
        iline, norbits);
#endif

/* number of orbits per object: 
* 0=all 1=last 2=last two orbits, etc.
*/
 
if(*found) {
if(norbits_per_object > 0) 
  imin = MAXI(norbits - norbits_per_object, 0);
else
  imin = 0;

for(i = norbits - 1; i >= imin; i--)
   fprintf(fp_out, "%s", &object_orbits[i*300]);
}

fclose(fp_in);
return(0);
}
/***************************************************************************
* line_extraction_from_OC6_catalog ("orb6.master" in November 2009)
* Return the line in OC6 that corresponds to the input object name 
*
* INPUT:
* OC6_fname: name of master file or Sixth Orbit catalog 
* is_master_file: flag set to 1 if master file ("orb6.master", 
*                          to 0 if OC6 file ("orb6orbits.txt")
* discov_name[60]: discovery name of object 
* comp_name: name of the companion (e.g., AB, Aa, etc) 
* norbits_per_object: maximum number of orbits to be loaded to output list
*                     fror each object
* 
* OUTPUT:
* fp_out: pointer to the output ASCII file 
* found: flag set to one if at least one orbit was found for object
*
***************************************************************************/
int line_extraction_from_OC6_catalog_gili(char *OC6_fname, int is_master_file,
                                     char *discov_name, char *comp_name, 
                                     FILE *fp_out, int *found, 
                                     int *candidate_found, 
                                     int norbits_per_object)
{
/* Line length is 278 + "\n" for OC6 catalog (orb6.master, november 2009)... */
/* Line length is 264 + "\n" for OC6 catalog (orb6orbits.txt, november 2009)... */
int iline, status, max_norbits = 1024, norbits, imin, nlines_in_header;
int OC6_comp_is_AB, comp_is_AB, discov_name_only;
char OC6_comp_really_compacted[40], comp_really_compacted[40];
char OC6_comp_name[40], OC6_discov_name[40];
char line_buffer[300], compacted_comp_name[40];
char compacted_discov_name[40];
/* Assume that the number of orbits for one object is always less than 1024: */
char object_orbits[300*1024];
FILE *fp_in;
int i;

*found = 0;
*candidate_found = 0;
norbits = 0;

discov_name_only = 1;

strcpy(compacted_comp_name, comp_name);
jlp_compact_string(compacted_comp_name, 40);
jlp_really_compact_companion(compacted_comp_name, comp_really_compacted, 40);

strcpy(compacted_discov_name, discov_name);
jlp_compact_string(compacted_discov_name, 40);

#ifdef DEBUG
 printf("(line_ext_OC6_cat_gili) CURRENT OBJECT: comp_name=%s discov_name=%s \n", 
 comp_name, discov_name);
#endif

/* Open OC6 catalog: */
if((fp_in = fopen(OC6_fname,"r")) == NULL) {
  fprintf(stderr, "line_extraction_from_OC6_catalog_gili/Fatal error opening %s\n",
          OC6_fname);
  exit(-1);
  }

iline = 0;

if(is_master_file)
  nlines_in_header = 4;
else
  nlines_in_header = 8;

while(!feof(fp_in)) {
/* Should read more than 
* 278 characters in master OC6 file
* 264 characters in non-master OC6 file
* in order to be sure to copy the complete line: */
  if(fgets(line_buffer, 280, fp_in)) {
    iline++;
    if(*found && !strncmp(line_buffer,"    ",4)) break; 
/* Skip the header and empty lines (that
* generally indicate the end of a given object):*/
    if(iline > nlines_in_header && strncmp(line_buffer,"    ",4)) {
/* Get the object name of each line */
     status = get_name_from_OC6_line_gili(line_buffer, OC6_discov_name, 
                                          OC6_comp_name);
     if(status) {
     fprintf(stderr, "line_extraction_from_OC6_catalog/Error processing line #%d\n", iline); 
     return(-1);
     }
#ifdef DEBUG_1
if(iline < nlines_in_header+3) 
     printf("OK1found=%d iline=%d: >%s<\n disc=%s< comp=%s<\n", 
             *found, iline, line_buffer, OC6_discov_name, OC6_comp_name);
#endif
    if(!strcmp(OC6_discov_name, compacted_discov_name)) {

/* Test on the companion names if present in the object name:
*/
     jlp_compact_string(OC6_comp_name, 40);
     jlp_really_compact_companion(OC6_comp_name, OC6_comp_really_compacted, 40);
#ifdef DEBUG_1
printf("DEBUG/OC6_object=%s OC6_companion=%s| comp=%s| (really compacted: OC6 comp=%s comp=%s)\n",
        OC6_discov_name, OC6_comp_name, comp_name, OC6_comp_really_compacted,
        comp_really_compacted);
#endif
        comp_is_AB = 0;
        if((compacted_comp_name[0] == '\0')
           || !strcmp(compacted_comp_name,"AB")
           || !strncmp(compacted_comp_name,"Aa-B",4)) comp_is_AB = 1;
        OC6_comp_is_AB = 0;
        if((OC6_comp_name[0] == '\0') || !strcmp(OC6_comp_name,"AB")
           || !strncmp(OC6_comp_name,"Aa-B",4)) OC6_comp_is_AB = 1;

        if((*compacted_comp_name != '\0'
             && !strcmp(OC6_comp_name, compacted_comp_name))
/* If not mentionned in Latex calibrated table, should
* be either not mentioned in OC6 or equal to AB: */
           || (comp_is_AB && OC6_comp_is_AB)
           || (*comp_really_compacted != '\0' &&
               !strcmp(comp_really_compacted, OC6_comp_really_compacted))){
        for(i = 0; i < 300; i++)
           object_orbits[i + 300*norbits] = line_buffer[i];
        *found = 1;
        norbits++;
        if(norbits >= max_norbits) {
            fprintf(stderr, "line_extraction_from_OC6_catalog/Fatal error: %d orbits found for a single object!\n", norbits);
            fprintf(stderr, "(OC6_discov_name=%s)\n", OC6_discov_name);
            exit(-1);
            }
         }
      }
    } /* EOF if line > nlines_header */
  } /* EOF if fgets */ 
 }

#ifdef DEBUG
printf("line_extraction_from_OC6_catalog: %d lines read and %d orbits found for current object\n", 
        iline, norbits);
#endif

/* number of orbits per object: 
* 0=all 1=last 2=last two orbits, etc.
*/
 
if(*found) {
if(norbits_per_object > 0) 
  imin = MAXI(norbits - norbits_per_object, 0);
else
  imin = 0;

for(i = norbits - 1; i >= imin; i--)
   fprintf(fp_out, "%s", &object_orbits[i*300]);
}

fclose(fp_in);
return(0);
}
/***************************************************************************
* get_name_from_OC6_line
* Extract the object and companion names from the line in OC6 
* that corresponds to the input object name 
*
*
* Example:
*
* 000000.00+000000.0 00093+7943 STF   2          102    431    760   6.68   6.89    540.      y    .         0.995  a   .      110.1       .     171.2        .      1887.5     y    .       0.715     .       333.7       .  
   2000      3 n Hei1997  abcdefghikj.png
*
*
* INPUT:
* in_line: full line of OC6 corresponding to the object
* 
* OUTPUT:
* OC6_ads_name: ADS name of object (NULL if not in ADS catalog) 
* OC6_comp_name: name of the companion (e.g., AB, Aa, etc) 
*
***************************************************************************/
int get_name_from_OC6_line(char *in_line, char *OC6_ads_name,
                           char *OC6_discov_name, char *OC6_comp_name) 
{
char WDS_name[40], ads_name[40], discov_name[40], comp_name[40];

/* Retrieve the object designation and the author of the orbit
*/
strncpy(WDS_name, &in_line[19], 10);
WDS_name[10] = '\0';

/* Discover's name in fixed format with 7 characters */
strncpy(discov_name, &in_line[30], 7);
discov_name[7] = '\0';
/* Remove all blanks: */
jlp_compact_string(discov_name, 10);
strcpy(OC6_discov_name, discov_name);

/* Companion name */
strncpy(comp_name, &in_line[37], 7);
comp_name[7] = '\0';
/* Remove all blanks: */
jlp_compact_string(comp_name, 10);
strcpy(OC6_comp_name, comp_name);

/* ADS name */
OC6_ads_name[0] = '\0';
strncpy(ads_name, &in_line[45], 5);
ads_name[5] = '\0';
/* Remove all blanks: */
jlp_compact_string(ads_name, 40);
if(ads_name[0] == '.') ads_name[0] = '\0'; 
if(ads_name[0] != '\0') sprintf(OC6_ads_name, "ADS %s", ads_name);

#ifdef DEBUG_1
  printf("WDS=%s OC6_ADS=%s OC6_discov=%s OC6_comp=%s\n",
          WDS_name, OC6_ads_name, OC6_discov_name, OC6_comp_name);
#endif

return(0);
}
/***************************************************************************
* get_name_from_OC6_line_gili
* Extract the object and companion names from the line in OC6 
* that corresponds to the input object name 
*
*
* Example:
*
* 000000.00+000000.0 00093+7943 STF   2          102    431    760   6.68   6.89    540.      y    .         0.995  a   .      110.1       .     171.2        .      1887.5     y    .       0.715     .       333.7       .  
   2000      3 n Hei1997  abcdefghikj.png
*
*
* INPUT:
* in_line: full line of OC6 corresponding to the object
* 
* OUTPUT:
* OC6_discov_name: name of the discover of the object
* OC6_comp_name: name of the companion (e.g., AB, Aa, etc) 
*
***************************************************************************/
int get_name_from_OC6_line_gili(char *in_line, char *OC6_discov_name, 
                                char *OC6_comp_name) 
{
char WDS_name[40], discov_name[40], comp_name[40];

/* Retrieve the object designation and the author of the orbit
*/
strncpy(WDS_name, &in_line[19], 10);
WDS_name[10] = '\0';

/* Discover's name in fixed format with 7 characters */
strncpy(discov_name, &in_line[30], 7);
discov_name[7] = '\0';
/* Remove all blanks: */
jlp_compact_string(discov_name, 10);
strcpy(OC6_discov_name, discov_name);

/* Companion name */
strncpy(comp_name, &in_line[37], 7);
comp_name[7] = '\0';
/* Remove all blanks: */
jlp_compact_string(comp_name, 10);
strcpy(OC6_comp_name, comp_name);

#ifdef DEBUG_1
  printf("WDS=%s OC6_OC6_discov=%s OC6_comp=%s\n",
          WDS_name, OC6_discov_name, OC6_comp_name);
#endif

return(0);
}
/************************************************************************
* Fill the reference table, with the references of "object_name"
*
* INPUT:
* object_name: e.g. ADS 2345 or COU 128
* author: e.g.  Sca1983c
* OC6_references_fname: name of filename with OC6 references
*
* OUPUT:
* refer0: e.g.: Scardia et al (1983c)
* refer1: e.g.: Scardia, M., Pansecchi, L., Sala, M., 1983, Apj 234, 56
*
************************************************************************/
int get_OC6_full_reference(char *object_name, char *author,
                           char *OC6_references_fname,
                           char *refer0, char *refer1)
{
int found;
char in_line0[130], author0[20], compacted_author[30];
#ifdef TTT
char decoded_authors[130], decoded_year[130];
#endif
char reference0[130], reference1[130];
FILE *fp_OC6_ref;

refer0[0] = '\0';
refer1[0] = '\0';

/* Removes all blanks contained in "author": */
strncpy(compacted_author, author, 30);
jlp_compact_string(compacted_author, 30);
strcpy(author, compacted_author);

/* Open OC6 references file: */
if((fp_OC6_ref = fopen(OC6_references_fname,"r")) == NULL) {
  fprintf(stderr, "get_OC6_full_reference/Fatal error error opening reference file >%s<\n",
          OC6_references_fname);
  exit(-1);
  }

/* Look for the reference associated to this author */
/* Example:
01234567890123456789012345678901234567890123456789012345678901234567890123456789
ABH  AbH2000b  F    Abt, H.A. & Corbally, C.J.                       ApJ 541, 841, 2000
*/
found = 0;
while(!feof(fp_OC6_ref)) {
  if(fgets(in_line0, 130, fp_OC6_ref)) {
   strncpy(author0, &in_line0[5], 9);
   jlp_compact_string(author0, 9);
    if(!strcmp(author0, compacted_author)){
         found = 1;
         strncpy(reference0, &in_line0[20], 48);
         reference0[47] = '\0';
/* Removes non-printable characters, heading, trailing and successive blanks */
         jlp_trim_string(reference0, 48);
         strcpy(reference1, &in_line0[68]);
/* Removes non-printable characters, heading, trailing and successive blanks */
         jlp_trim_string(reference1, 130);
         break;
         } 
  } /* EOF if(fgets) */
} /* EOF while feof */
 
if(found){
#ifdef TTT
  decode_authors_in_reference(reference0, decoded_authors, 130);
  decode_year_in_short_reference(author, decoded_year, 9);
  fprintf(fp_out_ref1, "%s = %s (%s),\n", author, decoded_authors, decoded_year); 
  fprintf(fp_out_ref2, "\\bibitem{%s} \n", author); 
  fprintf(fp_out_ref2, "\\fullref\\n%s\n%s = %s (%s),\n} \n", 
                       object_name, author, decoded_authors, decoded_year);
  remove_year_in_long_reference(reference1, 130);
  fprintf(fp_out_ref2, "%s : %s, %s\n\n", reference0, decoded_year, reference1);
#endif
strcpy(refer0, reference0);
strcpy(refer1, reference1);
  } else {
  fprintf(stderr, "reference not found for %s (%s)\n", 
          object_name, author);
  }

fclose(fp_OC6_ref);
return(0);
}
/************************************************************************
* Extract the authors of a publication from the long OC6 reference
* Example:
* Abt. H.A. & Biggs, E.S.
* Abt, H.A. & Levy, S.G.
* AbT, H.A. et al.
* Abt, H.A.
* Abt, H.A. Jeffers, H., Gibson, & Sandage, A.
* Astrographic Catalog (Greenwich)
************************************************************************/
static int decode_authors_in_reference(char *reference0, char *decoded_authors,
                                       int length0)
{
char author1[80], buffer[80], author2[80], etal[20];
char *pc;

strcpy(etal, "et al.\n");
decoded_authors[0] = '\0';

strncpy(author1, reference0, 40);
pc = author1;
author1[39] = '\0';
while(*pc && *pc != '.' && *pc != ',') pc++;
*pc = '\0';

author2[0] = '\0';
strncpy(buffer, reference0, 40);
pc = buffer;
buffer[39] = '\0';
buffer[79] = '\0';
while(*pc && *pc != '.' && *pc != ',') pc++;
while(*pc && *pc == ' ') pc++; 
if(*pc){
  if(isalpha(*pc) && *(pc+1) =='.') {
    pc+=2;
    if(isalpha(*pc) && *(pc+1) =='.') {
    pc+=2;
    strncpy(author2, pc, 40);
    }
  }
}

sprintf(decoded_authors, "%s %s", author1, author2);
return(0);
}
/************************************************************************
* Extract the year of publication from the short OC6 reference
* Example:
* Cou1973b   -> 1973b
************************************************************************/
static int decode_year_in_short_reference(char *author0, char *decoded_year,
                                          int length0)
{
char *pc;
decoded_year[0] = '\0';
author0[length0] = '\0';
pc = author0;
while(*pc && !isdigit(*pc)) pc++; 
if(isdigit(*pc)) strcpy(decoded_year, pc);
return(0);
}
/************************************************************************
* Removes the year of publication at the end of the long OC6 reference
* Example:
* Olevic, D. & Jovanovic, P. : Serbian AJ 163, 5, 2001
* -> Olevic, D. & Jovanovic, P. : Serbian AJ 163, 5
************************************************************************/
static int remove_year_in_long_reference(char *long_refer, int length1)
{
int i;
for(i = 0; i < length1; i++)
  if(long_refer[i] == '\0') break;

while((i > 0) && (long_refer[i] != ',')) i--; 
  if(long_refer[i] == ',') long_refer[i] = '\0'; 

return(0);
}
/****************************************************************************
*
* INPUT:
*  fp_out_ref1: pointer to output file with short references 
*  fp_out_ref2: pointer to output file with long references 
*  refer0, refer1: references extracted from OC6 reference file 
*  author: e.g. Scard1983c
*  object_name: e.g. ADS 2345 or COU 1342
*
****************************************************************************/
int save_references_for_residuals(FILE *fp_out_ref1, FILE *fp_out_ref2,
                                  char *object_name, char *author,
                                  char *refer0, char *refer1)
{
char decoded_authors[130], decoded_year[130];

if(fp_out_ref1 == NULL || fp_out_ref2 == NULL ) {
  fprintf(stderr, "save_references_for_residuals/Error: reference files not opened!\n");
  return(-1);
  }
  decode_authors_in_reference(refer0, decoded_authors, 130);
  decode_year_in_short_reference(author, decoded_year, 9);
  fprintf(fp_out_ref1, "%s = %s (%s),\n", author, decoded_authors, decoded_year);
  fprintf(fp_out_ref2, "\\bibitem{%s} \n", author);
  fprintf(fp_out_ref2, "\\fullref\\n%s\n%s = %s (%s),\n} \n",
                       object_name, author, decoded_authors, decoded_year);
  remove_year_in_long_reference(refer1, 130);
  fprintf(fp_out_ref2, "%s : %s, %s\n\n", refer0, decoded_year, refer1);
return(0);
}
/*************************************************************************
* Sort the references by alphabetic order 
*
* INPUT:
*  fp_out_ref1: pointer to output file with short references 
*  fp_out_ref2: pointer to output file with long references 
*  author[60*1024]: e.g. Scard1983c
*  object_name[60*1024]: e.g. ADS 2345 or COU 1342
*  refer1[130*1024], refer1[130*1024] : references extracted 
*                                       from the OC6 reference file 
*
**************************************************************************/
int sort_references(FILE *fp_out_ref1, FILE *fp_out_ref2, char *object_name, 
                    char *author, char *refer0, char *refer1, int n_names)
{
int k, length, j, status;
int index_ref[1024];
char buffer[60*1024];

if(fp_out_ref1 == NULL || fp_out_ref2 == NULL ) {
  fprintf(stderr, "sort_references/Error: reference files not opened!\n");
  return(-1);
  }

length = 60;

for(j = 0; j < 60 * 1024; j++) buffer[j] = author[j];

/* Sort "author" array by alphabetic order: */
JLP_QSORT_INDX_CHAR(buffer, &length, index_ref, &n_names);

for(j = 0; j < n_names; j++) {
   k = index_ref[j];
   status = save_references_for_residuals(fp_out_ref1, fp_out_ref2, &object_name[k * 60],
                                 &author[k * 60], &refer0[k * 130], 
                                 &refer1[k * 130]);
   if(status) return(-1);
   }
return(0);
}

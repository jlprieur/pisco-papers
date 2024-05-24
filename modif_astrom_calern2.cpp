/************************************************************************
* "modif_astrom_calern2.cpp"
* To modify the astrom file of measurements made by Gdpisco 
*
* 1. Removes Dm=0.05+\-0.02
* 2. If only EP=xxxx, add EJUL=xxxx:  EP=2017.0572 EJUL=2017.0555 
* 3. Add quadrant from LogFile.csv
* 4. Correct eyepiece if needed using LogFile.csv
* 5. Add unresolved objects
* 6. Correct all the bessellian epochs from the epochs of LogFile.csv
*
* JLP 
* Version 21/05/2024
*************************************************************************/
#include "jlp_catalog_utils.h"
#include "latex_utils.h" // jlp: latex_read_fvalue...
#include "csv_utils.h"  
#include "jlp_string.h"  // (in jlplib/jlp_fits/ ) jlp_cleanup_..
#include "jlp_fitsio.h"  // BesselToJulian

#define DEBUG
#define DEBUG_1

/***************************************************************************
* keywd: EP
* EP=2017.542
****************************************************************************/
static int delete_keywd_and_value_from_string(char *in_str, char *keywd,
                                             char *out_str)
{
char buffer[64], *pc2, *pc4, *pc_keywd;
int status = -1;

// To allow same argument (out_str=in_str)
 strcpy(buffer, in_str);
// Search for substring keywd in string in_str:
 pc_keywd = strstr(buffer, keywd);

// if not found, copy the full string and return from here
 if(pc_keywd == NULL){
   strcpy(out_str, in_str);
// if found, copy only the useful part: 
 } else {
   pc2 = buffer;
   pc4 = out_str;
   while(*pc2) {
     if(pc2 == pc_keywd) {
// Skip the in_str string until ' ' or '\\' is found
       while(*pc2  && (*pc2 != ' ') && (*pc2 != '\\')) pc2++;
       status = 0;
     } else {
// Copy in_str to out_str:
      *pc4 = *pc2;
      pc4++;
      pc2++;
      }
   }
// End character should be set to zero:
  *pc4 = '\0';
 }

return(status);
}
/***************************************************************************
* keywd: EP
* EP=2017.542
****************************************************************************/
static int read_keywd_dvalue_in_string(char *in_str, char *keywd, double *dvalue)
{
char *pc2, *pc3, *pc_keywd, kwd_arg[64];
int status = -1;

*dvalue = 0.;

// Search for substring keywd in string in_str:
 pc_keywd = strstr(in_str, keywd);
 if(pc_keywd != NULL){
   pc3 = kwd_arg;
   pc2 = in_str;
   while(*pc2) {
     if(pc2 == pc_keywd) {
// Skip the in_str string until '=' is found
          while(*pc2  && (*pc2 != '=')) pc2++;
          if(*pc2 == '=') pc2++;
// Skip the in_str string until ' ' or '\\' is found
          while(*pc2  && (*pc2 != ' ') && (*pc2 != '\\'))
            {
// Copy in_str to kwd_arg:
             *pc3 = *pc2;
             pc3++;
             pc2++;
            }
          *pc3 = '\0';
      sscanf(kwd_arg, "%lf", dvalue); 
      status = 0;
      return(status);
     }
   pc2++;
  }
 }

return(status);
}
/****************************************************************
* "2017,256" to  "2017.256" 
****************************************************************/
static int convert_coma_to_dot(char *in_str, char *out_str)
{
char *pc1, *pc2;
pc1 = in_str;
pc2 = out_str;
while(*pc1) {
  if(*pc1 != ',') *pc2 = *pc1;
  else *pc2 = '.';
  pc1++;
  pc2++;
  }
return(0);
}
static int modif_astrom_removeDm(char *in_file1, char *out_fname);
static int modif_astrom_Bessel_to_Julian(char *in_file1, char *out_fname);
static int modif_astrom_check_Bessel(char *in_file1, char *csv_LogFile,
                                     char *out_fname);
static int modif_astrom_add_quadrant(char *in_file1, char *csv_LogFile, 
                                     char *out_fname);
static int modif_astrom_check_eyepiece(char *in_file1, char *csv_LogFile, 
                                       char *out_fname);
static int modif_astrom_add_unresolved(char *in_file1, char *csv_LogFile, 
                                       char *out_fname);
static int decode_fname_calern(char *autoc_fname1, char *dble_name1, 
                               char *filter1, int *day1,
                               int *month1, int *year1);
static int read_from_LogFile(char *csv_LogFile, char *dble_name1, 
                             char *filter1, int day1, int month1, int year1,
                             int *eyepiece2, char *quad2, 
                             double *ep_bessel2, double *ep_julian2, 
                             int *ifound);
static int read_LogFile(char *csv_LogFile, char *dble_name2, 
                        char *filter2, char *quadrant2, double *ep_bessel2,
                        double *ep_julian2, int *day2, 
                        int *month2, int *year2, int *eyep2, int *nitems2);

int main(int argc, char *argv[])
{
char in_file1[128], csv_LogFile[128], out_fname[128];
int iopt;

/* If command line with "runs" */
if(argc == 7){
 if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
 else argc = 1;
 }

if(argc != 5) {
  printf("Syntax: modif_astrom_calern2 iopt in_astrom_file csv_LogFile out_astrom_file\n");
  return(-1);
}
sscanf(argv[1], "%d", &iopt);
strcpy(in_file1, argv[2]);
strcpy(csv_LogFile, argv[3]);
strcpy(out_fname, argv[4]);

printf("OK: iopt=%d in_file1=%s csv_LogFile=%s out_fname=%s \n", 
        iopt, in_file1, csv_LogFile, out_fname); 

// 1. Removes Dm=0.05+\-0.02
switch(iopt)
 {
  case 1:
    modif_astrom_removeDm(in_file1, out_fname); 
    break;
// 2. If only EP=xxxx, add EJUL=xxxx:  EP=2017.0572 EJUL=2017.0555 
  case 2:
    modif_astrom_Bessel_to_Julian(in_file1, out_fname);
    break;
  case 3:
    modif_astrom_add_quadrant(in_file1, csv_LogFile, out_fname);
    break;
  case 4:
    modif_astrom_check_eyepiece(in_file1, csv_LogFile, out_fname);
    break;
  case 5:
    modif_astrom_add_unresolved(in_file1, csv_LogFile, out_fname);
    break;
// 6. Correct all the bessellian epochs from the epochs of LogFile.csv
  case 6:
    modif_astrom_check_Bessel(in_file1, csv_LogFile, out_fname);
    break;
 }

return(0);
}
/************************************************************************
* Scan the input astrom file and make the modifications 
*
* INPUT:
*   in_file1: name of the input file 
*   out_fname: name of the output file
*
*************************************************************************/
static int modif_astrom_removeDm(char *in_file1, char *out_fname) 
{
char in_line[256], in_line2[256]; 
char buffer[256], Dm_string[64], *pc1, *pc2, *pc_Dm;
int iline;
FILE *fp_in1, *fp_out;
time_t t0 = time(NULL);

/* Open input astrom file: */
if((fp_in1 = fopen(in_file1, "r")) == NULL) {
   fprintf(stderr, "modif_astrom_calern2/Fatal error opening input file %s\n",
           in_file1);
    return(-1);
  }

/* Open output file: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "modif_astrom_calern2/Fatal error opening output file: %s\n",
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

strcpy(Dm_string, "Dm=");
iline = 0;
while(!feof(fp_in1)) {
  if(fgets(in_line, 256, fp_in1)) {
    iline++;
// Remove the end of line '\n' from input line:
//    jlp_cleanup_string(in_line, 256);

    strcpy(in_line2, in_line);
    if(in_line2[0] == '&') {
// Search for substring "Dm_string" in string "in_line2":
      pc_Dm = strstr(in_line2, Dm_string); 
      if(pc_Dm != NULL) {
        pc1 = buffer;
        pc2 = in_line2;
        while(*pc2) {
         if(pc2 == pc_Dm) {
// Skip the Dm string until ' ' or '\\' is found
          while(*pc2  && (*pc2 != ' ') && (*pc2 != '\\')) pc2++;
          } else {
// Copy in_line2 to buffer:
          *pc1 = *pc2;
          pc2++;
          pc1++;
          }
        } // EOF while(*pc2)
        *pc1 = '\0';
       strcpy(in_line2, buffer);
      } // EOF pc_Dm != NULL
    } //EOF in_line2[0] == '&')

// Save to output file:
    fprintf(fp_out, "%s", in_line2);

  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("modif_astrom_removeDm: %d lines successfully read and processed\n", 
        iline);

/* Close opened files:
*/
fclose(fp_in1);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input astrom file and make the modifications 
*
* INPUT:
*   in_file1: name of the input file 
*   out_fname: name of the output file
*
*************************************************************************/
static int modif_astrom_Bessel_to_Julian(char *in_file1, char *out_fname)
{
char in_line[256], in_line2[256]; 
char buffer[256], EP_string[64], EJUL_string[64], EP_arg[64];
char *pc1, *pc2, *pc3, *pc_EP, *pc_EJUL;
double BesselEpoch, JulianEpoch;
int iline;
FILE *fp_in1, *fp_out;
time_t t0 = time(NULL);

/* Open input astrom file: */
if((fp_in1 = fopen(in_file1, "r")) == NULL) {
   fprintf(stderr, "modif_astrom_calern2/Fatal error opening input file %s\n",
           in_file1);
    return(-1);
  }

/* Open output file: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "modif_astrom_calern2/Fatal error opening output file: %s\n",
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

// If only EP=xxxx, add EJUL=xxxx:  EP=2017.0572 EJUL=2017.0555 
strcpy(EP_string, "EP=");
strcpy(EJUL_string, "EJUL=");
iline = 0;
while(!feof(fp_in1)) {
  if(fgets(in_line, 256, fp_in1)) {
    iline++;
// Remove the end of line '\n' from input line:
//    jlp_cleanup_string(in_line, 256);

    strcpy(in_line2, in_line);
    if(in_line2[0] == '&') {
// Search for substring "EP_string" in string "in_line2":
      pc_EP = strstr(in_line2, EP_string); 
      pc_EJUL = strstr(in_line2, EJUL_string); 
      BesselEpoch = 0.;
      JulianEpoch = 0.;
      if((pc_EP != NULL) && (pc_EJUL == NULL)) {
// Look for EP_arg (i.e, all the string after EP, like "EP=2017.564 Q=3?")
        pc3 = EP_arg;
        pc1 = buffer;
        pc2 = in_line2;
        while(*pc2) {
         if(pc2 == pc_EP) {
// Skip the EP string until ' ' or '\\' is found
          while(*pc2  && (*pc2 != ' ') && (*pc2 != '\\')) 
            { 
// Copy in_line2 to EP_arg:
             *pc3 = *pc2;
// Copy in_line2 to buffer:
             *pc1 = *pc2;
             pc1++;
             pc2++;
             pc3++;
            }
// Stop EP_arg here:
          *pc3 = '\0';
           sscanf(EP_arg, "EP=%lf", &BesselEpoch);
// JLP_besselian_to_julian_epoch(double b_date, double *j_date);
           JLP_besselian_to_julian_epoch(BesselEpoch, &JulianEpoch);
           sprintf(pc1, " EJUL=%.4f", JulianEpoch);
           pc1 +=15;
          } else {
// Copy in_line2 to buffer:
          *pc1 = *pc2;
          pc2++;
          pc1++;
          }
        } // EOF while(*pc2)
        *pc1 = '\0';
       strcpy(in_line2, buffer);
      } // EOF pc_EP != NULL
    } //EOF in_line2[0] == '&')

// Save to output file:
    fprintf(fp_out, "%s", in_line2);

  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("modif_astrom_Julian: %d lines successfully read and processed\n", 
        iline);

/* Close opened files:
*/
fclose(fp_in1);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input astrom file and make the modifications 
*
* INPUT:
*   in_file1: name of the input file 
*   out_fname: name of the output file
*
*************************************************************************/
static int modif_astrom_add_quadrant(char *in_file1, char *csv_LogFile, 
                                     char *out_fname)
{
char in_line[256], in_line2[256], comments1[64]; 
char buffer[256], Quad_string[64], autoc_fname1[64];
char *pc1, *pc2, *pc_Q, quad2[64];
int status, iline, icol, ieyepiece, verbose_if_error=0; 
int eyepiece2, ifound, nfail, day1, month1, year1; 
double ep_bessel2, ep_julian2;
char dble_name1[64], filter1[64];
FILE *fp_in1, *fp_out;
time_t t0 = time(NULL);

/* Open input astrom file: */
if((fp_in1 = fopen(in_file1, "r")) == NULL) {
   fprintf(stderr, "modif_astrom_calern2/Fatal error opening input file %s\n",
           in_file1);
    return(-1);
  }

/* Open output file: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "modif_astrom_calern2/Fatal error opening output file: %s\n",
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

strcpy(Quad_string, "Q=");
iline = 0;
nfail = 0;
while(!feof(fp_in1)) {
  if(fgets(in_line, 256, fp_in1)) {
    iline++;
// Remove the end of line '\n' from input line:
    jlp_cleanup_string(in_line, 256);

    strcpy(in_line2, in_line);
    if(in_line2[0] == '&') {
// Look for eyepiece:
     icol = 5;
     status = latex_read_ivalue(in_line2, &ieyepiece, icol);
     if(status != 0) ieyepiece = -1;
#ifdef DEBUG_1
     printf("modif_astrom_add_quadrant/ >%s< eyepiece=%d\n", in_line2,ieyepiece);
#endif
// Look for filename:
     icol = 2;
     status = latex_read_svalue(in_line2, autoc_fname1, icol);
#ifdef DEBUG_1
     printf("modif_astrom_add_quadrant/filename=%s\n", autoc_fname1);
#endif
// Look for comments:
     icol = 10;
     status = latex_read_svalue(in_line2, comments1, icol);
#ifdef DEBUG_1
     printf("modif_astrom_add_quadrant/comments=%s\n", comments1);
#endif
  
     if(ieyepiece > -1) {
      status = decode_fname_calern(autoc_fname1, dble_name1, filter1, 
                                   &day1, &month1, &year1);
      read_from_LogFile(csv_LogFile, dble_name1, filter1, day1, month1, year1,
                        &eyepiece2, quad2, &ep_bessel2, &ep_julian2, &ifound);
      if(ifound == 0) {
        nfail++;
        fprintf(stderr, "Error: >%s< not found in LogFile (nfail=%d)\n", dble_name1, nfail);
        } else {
#ifdef DEBUG_1
     printf("dble_name1=%s filter1=%s eyepiece2=%d quad2=%s ifound=%d\n", 
             dble_name1, filter1, eyepiece2, quad2, ifound);
#endif
       }
// Search for substring "Quad_string" in string "in_line2":
      pc_Q = strstr(comments1, Quad_string); 
      if(pc_Q == NULL) {
        pc1 = buffer;
        pc2 = comments1;
        while(*pc2) 
         { 
// Copy comments1 to buffer:
          *pc1 = *pc2;
          pc2++;
          pc1++;
          } // EOF while(*pc2)
          if(ifound == 1) {
            sprintf(pc1, " Q=%s ", quad2);
          } else {
            sprintf(pc1, " QuadToBeFound ");
          }
        icol = 10;
        latex_set_column_item(in_line2, 256, buffer, 256, icol, verbose_if_error); 
      } // EOF pc_Q == NULL
     } // ieyepiece > -1
    } //EOF in_line2[0] == '&')

// Save to output file:
    fprintf(fp_out, "%s\n", in_line2);

  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("modif_astrom_quad: %d lines successfully read and processed, nfail=%d\n", 
        iline, nfail);

/* Close opened files:
*/
fclose(fp_in1);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input astrom file and make the modifications 
*
* INPUT:
*   in_file1: name of the input file 
*   out_fname: name of the output file
*
*************************************************************************/
static int modif_astrom_check_Bessel(char *in_file1, char *csv_LogFile, 
                                     char *out_fname)
{
char in_line[256], in_line2[256], comments1[64], autoc_fname1[64]; 
char buffer[256], EP_string[64], EJUL_string[64], EP_arg[64], quad2[64];
char *pc1, *pc2, *pc3, *pc_EP, *pc_EJUL;
double BesselEpoch, JulianEpoch;
double ep_bessel2, ep_julian2;
int status, iline, icol, ieyepiece, verbose_if_error=0; 
int eyepiece2, ifound, nfail, day1, month1, year1; 
char dble_name1[64], filter1[64], autoc_fname[64];
FILE *fp_in1, *fp_out;
time_t t0 = time(NULL);

/* Open input astrom file: */
if((fp_in1 = fopen(in_file1, "r")) == NULL) {
   fprintf(stderr, "modif_astrom_calern2/Fatal error opening input file %s\n",
           in_file1);
    return(-1);
  }

/* Open output file: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "modif_astrom_calern2/Fatal error opening output file: %s\n",
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

// If only EP=xxxx, add EJUL=xxxx:  EP=2017.0572 EJUL=2017.0555 
strcpy(EP_string, "EP=");
strcpy(EJUL_string, "EJUL=");
iline = 0;
nfail = 0;
while(!feof(fp_in1)) {
  if(fgets(in_line, 256, fp_in1)) {
    iline++;
// Remove the end of line '\n' from input line:
    jlp_cleanup_string(in_line, 256);

    strcpy(in_line2, in_line);
    if(in_line2[0] == '&') {
// Look for eyepiece:
     icol = 5;
     status = latex_read_ivalue(in_line2, &ieyepiece, icol);
     if(status != 0) ieyepiece = -1;
#ifdef DEBUG_1
     printf("modif_astrom_check_bessel/ >%s< eyepiece=%d\n", in_line2,ieyepiece);
#endif
// Look for filename:
     icol = 2;
     status = latex_read_svalue(in_line2, autoc_fname1, icol);
#ifdef DEBUG_1
     printf("modif_astrom_add_quadrant/filename=%s\n", autoc_fname1);
#endif
// Look for comments:
     icol = 10;
     status = latex_read_svalue(in_line2, comments1, icol);
#ifdef DEBUG_1
     printf("modif_astrom_add_quadrant/comments=%s\n", comments1);
#endif
  
     if(ieyepiece > -1) {
      status = decode_fname_calern(autoc_fname1, dble_name1, filter1, 
                                   &day1, &month1, &year1);
      read_from_LogFile(csv_LogFile, dble_name1, filter1, day1, month1, year1,
                        &eyepiece2, quad2, &ep_bessel2, &ep_julian2, &ifound);
      if(ifound == 0) {
        nfail++;
        fprintf(stderr, "Error: >%s< not found in LogFile (nfail=%d)\n", dble_name1, nfail);
        } else {
#ifdef DEBUG_1
     printf("dble_name1=%s filter1=%s eyepiece2=%d quad2=%s ifound=%d\n", 
             dble_name1, filter1, eyepiece2, quad2, ifound);
#endif
       }
// Read EP and EP_JUL if present
     read_keywd_dvalue_in_string(comments1, EP_string, &BesselEpoch);
     read_keywd_dvalue_in_string(comments1, EJUL_string, &JulianEpoch);
     delete_keywd_and_value_from_string(comments1, EP_string, comments1);
     delete_keywd_and_value_from_string(comments1, EJUL_string, comments1);
// Update comments1:
        pc1 = buffer;
        pc2 = comments1;
        while(*pc2) 
         { 
// Copy comments1 to buffer:
          *pc1 = *pc2;
          pc2++;
          pc1++;
          } // EOF while(*pc2)
        if(ifound == 1) {
// JLP_besselian_to_julian_epoch(double b_date, double *j_date);
          JLP_besselian_to_julian_epoch(ep_bessel2, &ep_julian2);
//
          if(ABS(BesselEpoch - ep_bessel2) < 0.01) {
           sprintf(pc1, " EP=%.4f EJUL=%.4f", 
                   ep_bessel2, ep_julian2);
           } else {
           sprintf(pc1, " EP0=%.4f EJ0=%.4f EP2=%.4f EJUL2=%.4f", 
                  BesselEpoch, JulianEpoch, ep_bessel2, ep_julian2);
           }
        } else {
          sprintf(pc1, " EpochToBeFound ");
        }
      icol = 10;
      latex_set_column_item(in_line2, 256, buffer, 256, icol, verbose_if_error); 
     } // ieyepiece > -1
    } //EOF in_line2[0] == '&')

// Save to output file:
    fprintf(fp_out, "%s\n", in_line2);

  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("modif_astrom_quad: %d lines successfully read and processed, nfail=%d\n", 
        iline, nfail);

/* Close opened files:
*/
fclose(fp_in1);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input astrom file and make the modifications 
*
* INPUT:
*   in_file1: name of the input file 
*   out_fname: name of the output file
*
*************************************************************************/
static int modif_astrom_check_eyepiece(char *in_file1, char *csv_LogFile, 
                                       char *out_fname)
{
char in_line[256], in_line2[256], buffer[256]; 
char autoc_fname1[64], quad2[64], dble_name1[64], filter1[64];
double ep_bessel2, ep_julian2;
int status, iline, icol, ieyepiece, verbose_if_error=0; 
int eyepiece2, ifound, nfail, ncorr, day1, month1, year1; 
FILE *fp_in1, *fp_out;
time_t t0 = time(NULL);

/* Open input astrom file: */
if((fp_in1 = fopen(in_file1, "r")) == NULL) {
   fprintf(stderr, "modif_astrom_calern2/Fatal error opening input file %s\n",
           in_file1);
    return(-1);
  }

/* Open output file: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "modif_astrom_calern2/Fatal error opening output file: %s\n",
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

iline = 0;
nfail = 0;
ncorr = 0;
while(!feof(fp_in1)) {
  if(fgets(in_line, 256, fp_in1)) {
    iline++;
// Remove the end of line '\n' from input line:
    jlp_cleanup_string(in_line, 256);

    strcpy(in_line2, in_line);
    if(in_line2[0] == '&') {
// Look for eyepiece:
     icol = 5;
     status = latex_read_ivalue(in_line2, &ieyepiece, icol);
     if(status != 0) ieyepiece = -1;
#ifdef DEBUG_1
     printf("modif_astrom_check_eyepiece/ >%s< eyepiece=%d\n", in_line2,ieyepiece);
#endif
// Look for filename:
     icol = 2;
     status = latex_read_svalue(in_line2, autoc_fname1, icol);
#ifdef DEBUG_1
     printf("modif_astrom_check_eyepiece/filename=%s\n", autoc_fname1);
#endif
  
     if(ieyepiece == 0) {
      status = decode_fname_calern(autoc_fname1, dble_name1, filter1, 
                                   &day1, &month1, &year1);
      read_from_LogFile(csv_LogFile, dble_name1, filter1, day1, month1, year1,
                        &eyepiece2, quad2, &ep_bessel2, &ep_julian2, &ifound);
      if(ifound == 0) {
        nfail++;
        fprintf(stderr, "Error: >%s< not found in LogFile (nfail=%d)\n", dble_name1, nfail);
        sprintf(buffer, "4444 ");
      } else {
        printf("WWW NOW correct ieyepiece=%d eyepiece2=%d for %s\n", 
               ieyepiece, eyepiece2, dble_name1);
        sprintf(buffer, "%d ", eyepiece2);
        ncorr++;
      }
      icol = 5;
      latex_set_column_item(in_line2, 256, buffer, 256, icol, 
                              verbose_if_error); 
     } // ieyepiece == 0
    } //EOF in_line2[0] == '&')

// Save to output file:
    fprintf(fp_out, "%s\n", in_line2);

  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("modif_astrom_check_eyepiece: %d lines successfully read and processed, nfail=%d ncorr=%d\n", 
        iline, nfail, ncorr);

/* Close opened files:
*/
fclose(fp_in1);
fclose(fp_out);
return(0);
}
/************************************************************************
* Scan the input astrom file and make the modifications 
*
* INPUT:
*   in_file1: name of the input file 
*   out_fname: name of the output file
*
*************************************************************************/
static int modif_astrom_add_unresolved(char *in_file1, char *csv_LogFile, 
                                       char *out_fname)
{
char in_line[256], in_line2[256], buffer[256]; 
char autoc_fname1[64], quad2[64], dble_name1[64], filter1[64];
char *pc_Unres, Unres_string[64], *pc, *pc1;
int status, iline, icol, ieyepiece, verbose_if_error=0; 
int eyepiece2, ifound, nunres, day1, month1, year1; 
double BesselEpoch, JulianEpoch;
double ep_bessel2, ep_julian2;
FILE *fp_in1, *fp_out;
time_t t0 = time(NULL);

/* Open input astrom file: */
if((fp_in1 = fopen(in_file1, "r")) == NULL) {
   fprintf(stderr, "modif_astrom_calern2/Fatal error opening input file %s\n",
           in_file1);
    return(-1);
  }

/* Open output file: */
  if((fp_out = fopen(out_fname, "w")) == NULL) {
    fprintf(stderr, "modif_astrom_calern2/Fatal error opening output file: %s\n",
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

iline = 0;
nunres = 0;
while(!feof(fp_in1)) {
  if(fgets(in_line, 256, fp_in1)) {
    iline++;
// Remove the end of line '\n' from input line:
    jlp_cleanup_string(in_line, 256);

    strcpy(in_line2, in_line);
    if(in_line2[0] == '%') {
// Search for substring "Quad_string" in string "in_line2":
// "Pas resolu"  "pas resolu"
      strcpy(Unres_string, "as ");
      pc_Unres = strstr(in_line2, Unres_string);
      if(pc_Unres != NULL) {
         printf("WWW: %s\n", in_line2);
         strncpy(autoc_fname1, &in_line2[3], 64);
         pc = autoc_fname1;
         while(*pc && *pc != ' ') pc++;
         *pc = '\0';
// Change _ to \_ for latex:
         strcpy(buffer, autoc_fname1);
         pc1 = buffer;
         pc = autoc_fname1;
         while(*pc) {
           if(*pc == '_') {
             *pc1 = '\\';
             pc1++;
             }
           *pc1 = *pc;
           pc1++;
           pc++;
           }
         *pc1 = '\0';
         strcpy(autoc_fname1, buffer);
         status = decode_fname_calern(autoc_fname1, dble_name1, filter1, 
                                      &day1, &month1, &year1);
         printf("VVV: %s filter=%s %d-%d-%d\n", autoc_fname1, filter1, 
                day1, month1, year1);
         read_from_LogFile(csv_LogFile, dble_name1, filter1, day1, month1, 
                           year1, &eyepiece2, quad2, &ep_bessel2, &ep_julian2, 
                           &ifound);
         printf("VVV: ifound=%d %s eyepiece=%d Q=%s\n", 
                 ifound, autoc_fname1, eyepiece2, quad2);
         if((ifound == 0) || (eyepiece2 <= 0)) eyepiece2 = 20;
         printf("& %s & %d & & & & & & & \\\\\n", dble_name1, year1);
         fprintf(fp_out, "& %s & %d & & & & & & & \\\\\n", dble_name1, year1);
// JLP_besselian_epoch(double aa, int mm, int idd, double time, double *b_date)
         JLP_besselian_epoch((double)year1, month1, day1, 0., &BesselEpoch);
// JLP_besselian_to_julian_epoch(double b_date, double *j_date);
         JLP_besselian_to_julian_epoch(BesselEpoch, &JulianEpoch);
         nunres++;
         printf("& %s & %02d/%02d/%d & %s & %d & \\nodata & \\nodata & \\nodata & \\nodata & NR EP=%.4f EJUL=%.4f \\\\\n", 
                autoc_fname1, day1, month1, year1, filter1, eyepiece2, 
                BesselEpoch, JulianEpoch);
         fprintf(fp_out, "& %s & %02d/%02d/%d & %s & %d & \\nodata & \\nodata & \\nodata & \\nodata & NR EP=%.4f EJUL=%.4f \\\\\n", 
                autoc_fname1, day1, month1, year1, filter1, eyepiece2, 
                BesselEpoch, JulianEpoch);
         }
    } //EOF in_line2[0] == '%')

// Save to output file:
    fprintf(fp_out, "%s\n", in_line2);

  } /* EOF if fgets */ 
 } /* EOF while ... */
printf("modif_astrom_add_unresolved: %d lines successfully read and processed, nunres=%d\n", 
        iline, nunres);

/* Close opened files:
*/
fclose(fp_in1);
fclose(fp_out);
return(0);
}
/***************************************************************************
* Use filename to get the date of observation, the star name and the filter
*
* Example: 101117\_rst5196\_R\_a
***************************************************************************/
static int decode_fname_calern(char *autoc_fname1, char *dble_name1, 
                               char *filter1, int *day1,
                               int *month1, int *year1)
{
char *pc1, *pc, date1[64];
int len1, istart, status = 0;

 strcpy(date1, autoc_fname1);
 pc1 = autoc_fname1;
 pc = date1;
 while(*pc1 && (*pc1 != '\\') && (*pc1 != '_')) {
  if(*pc1 != ' ') {
    *pc = *pc1;
    pc++;
    }
  pc1++;
  }
 *pc = '\0';
  len1 = strlen(date1);
  if(len1 == 6) {
    sscanf(date1, "%02d%02d%02d", day1, month1, year1);
    } else {
    fprintf(stderr," decode_fname_calern/fatal error: date1=%s< len1 = %d\n", 
            date1, len1); 
    exit(-1);
    }
  *year1 += 2000;
#ifdef DEBUG_1
  printf("decode_fname_calern/date: %d %d %d\n", *day1, *month1, *year1);
#endif

// Example: 101117\_rst5196\_R\_a
 strcpy(dble_name1, autoc_fname1);
 pc = dble_name1;
 istart = 0;
 while(*pc1) {
  if(((*pc1 == '\\') || (*pc1 == '_')) && (istart == 1)) break;
  if(istart > 0) { *pc = *pc1; pc++;}
  if((*pc1 == '_') && (istart == 0)) istart = 1;
  pc1++;
  }
 *pc = '\0';
  pc = dble_name1;
  while(*pc) {*pc = toupper(*pc); pc++;}
#ifdef DEBUG_1
  printf("decode_fname_calern/dble_name: %s\n", dble_name1);
#endif
 
  if(strstr(autoc_fname1, "\\_V") != NULL) strcpy(filter1, "V"); 
  else if(strstr(autoc_fname1, "\\_R") != NULL) strcpy(filter1, "R"); 
  else if(strstr(autoc_fname1, "\\_RL") != NULL) strcpy(filter1, "RL"); 
  else if(strstr(autoc_fname1, "\\_I") != NULL) strcpy(filter1, "I"); 
  else if(strstr(autoc_fname1, "\\_W") != NULL) strcpy(filter1, "W"); 
  else strcpy(filter1, " "); 
#ifdef DEBUG_1
  printf("decode_fname_calern/filter: %s\n", filter1);
#endif

return(status);
}
#define NMAXI 1024
/*********************************************************************
*
**********************************************************************/
static int read_from_LogFile(char *csv_LogFile, char *dble_name1, 
                             char *filter1, int day1, int month1, int year1,
                             int *eyepiece2, char *quad2, 
                             double *ep_bessel2, double *ep_julian2, 
                             int *ifound)
{
char dble_name2[64 * NMAXI], filter2[16 * NMAXI], quadrant2[16 * NMAXI]; 
char buffer[64], *pc;
int day2[NMAXI], month2[NMAXI], year2[NMAXI], eyep2[NMAXI], nitems2;
int i, status = 0;
double EP_Bessel2[NMAXI], EP_Julian2[NMAXI];
static int initialized0 = 0;

strcpy(quad2, " ");
*eyepiece2 = -1;
*ep_bessel2 = -1.;
*ep_julian2 = -1.;
*ifound = 0;

if(initialized0 == 0) {
 read_LogFile(csv_LogFile, dble_name2, filter2, quadrant2,
              EP_Bessel2, EP_Julian2, day2, month2, year2, eyep2, &nitems2);
#ifdef DEBUG1
   printf("read_from_Logfile/nitems2=%d \n", nitems2);
   for(i = 0; i < nitems2; i++) {
      printf("read_from_Logfile/i=%d dble_name2=%s filter2=%s ep_bessel=%.4f ep_julian=%.4f\n", 
         i, dble_name2[16*i], filter2[16 * i], EP_Bessel2[i], EP_Julian2[i]);
     }
#endif
 initialized0 = 1;
 }

jlp_compact_string(dble_name1, 64);
jlp_compact_string(filter1, 64);

//printf("read_from_Logfile/dble_name1=%s< filter1=%s<\n", dble_name1, filter1);

for(i = 0; i < nitems2; i++)
  {
    strcpy(buffer, &dble_name2[64 * i]);
    pc = buffer;
    while(*pc) {*pc = toupper(*pc); pc++;}
/*
if((dble_name1[0] == 'A') && (i == 15))
  printf("ZZZ i=%d Star=%s< filter=%s< date: %d-%d-%d<\n", 
          i, &dble_name2[64 * i], &filter2[16 * i], day2[i], month2[i], year2[i]);
*/
  if(!strcmp(dble_name1, buffer)
     && !strcmp(filter1, &filter2[16 * i])
     && (ABS(day1 - day2[i]) < 2) && (month1 == month2[i]) && (year1 == year2[i])) {
     *eyepiece2 = eyep2[i];
     *ep_bessel2 = EP_Bessel2[i];
     *ep_julian2 = EP_Julian2[i];
     strcpy(quad2, &quadrant2[16 * i]);
     *ifound = 1;
     } 
  }
      if(*ifound != 1) {
        fprintf(stderr, "Error: Star=%s< Filter=%s< date: %d-%d-%d< not found in LogFile\n", dble_name1, filter1, day1, month1, year1);
        status = -1;
        }

return(status);
}
/*********************************************************************
*
**********************************************************************/
static int read_LogFile(char *csv_LogFile, char *dble_name2, 
                        char *filter2, char *quadrant2, double *EP_Bessel2,
                        double *EP_Julian2, int *day2, 
                        int *month2, int *year2, int *eyep2, int *nitems2)
{
char b_in[512], buffer[64], buffer2[64], keywd[64]; 
char dble_kd[64], filter_kd[64], quad_kd[64], eyepiece_kd[64], date_kd[64];
char bessel_kd[64], julian_kd[64];
int nitem2[NMAXI], i, status = 0, it, iline, iw;
double dw;
FILE *fp_in, *fp_out;

if((fp_in = fopen(csv_LogFile, "r")) == NULL) {
  fprintf(stderr, " Fatal error opening LogFile1 %s \n", csv_LogFile);
  return(-1);
  }

if((fp_out = fopen("logfile_tmp.txt", "w")) == NULL) {
  fprintf(stderr, " Fatal error opening out_logfile\n");
  return(-1);
  }

strcpy(date_kd, "Data");
strcpy(dble_kd, "DOPPIA");
strcpy(filter_kd, "Filtro");
strcpy(quad_kd, "Quadrante");
strcpy(eyepiece_kd, "Oculare");
strcpy(bessel_kd, "Epoca Besseliana");
strcpy(julian_kd, "Epoca Juliana");

for(i = 0; i < NMAXI; i++) { nitem2[i] = 0; }

iline = 0;
it = -1;
while(!feof(fp_in))
{
  if(fgets(b_in, 512, fp_in))
  {
  iline++;
#ifdef DEBUG_1
  printf("OKK: iline=%d\n", iline);
  printf("OKK: >%s<\n", b_in);
#endif
// int csv_read_string(char *b_data, int i_column, char *out_string);
  status = csv_read_string(b_in, 1, keywd);
// Date ("Data"):
  if(strncmp(date_kd, keywd, 4) == 0) {
    if(it >= NMAXI) {
      fprintf(stderr, "Fatal error: it=%d >= NMAXI=%d\n", it, NMAXI);
      exit(-1);
      }
    it++;
    status = csv_read_string(b_in, 2, buffer);
    sscanf(buffer, "%d", &iw);
    day2[it] = iw;
    status = csv_read_string(b_in, 3, buffer);
    sscanf(buffer, "%d", &iw);
    month2[it] = iw;
    status = csv_read_string(b_in, 4, buffer);
    sscanf(buffer, "%d", &iw);
    year2[it] = iw;
#ifdef DEBUG
    printf("read_from_LogFile: it=%d date2=%d %d %d\n", it, day2[it], month2[it], year2[it]);
#endif
    nitem2[it]++;
    } 
// Double star discover name ("DOPPIA"):
  else if((it >= 0) && (strncmp(dble_kd, keywd, 6) == 0)) {
    status = csv_read_string(b_in, 2, buffer);
    strcpy(&dble_name2[it * 64], buffer);
    jlp_compact_string(&dble_name2[it * 64], 64);
#ifdef DEBUG
    printf("read_from_LogFile: it=%d dble_name2=%s<\n", it, &dble_name2[it * 64]);
#endif
    nitem2[it]++;
    } 
// Filter ("Filtro"):
  else if((it >= 0) && (strncmp(filter_kd, keywd, 6) == 0)) {
    status = csv_read_string(b_in, 2, buffer);
    strcpy(&filter2[it * 16], buffer);
    jlp_compact_string(&filter2[it * 16], 16);
#ifdef DEBUG
    printf("read_from_LogFile: it=%d filter2=%s\n", it, &filter2[it * 16]);
#endif
    nitem2[it]++;
    } 
// Quadrant ("Quadrante"):
  else if((it >= 0) && (strncmp(quad_kd, keywd, 9) == 0)) {
    status = csv_read_string(b_in, 2, buffer);
    strcpy(&quadrant2[it * 16], buffer);
    jlp_compact_string(&quadrant2[it * 16], 16);
#ifdef DEBUG
    printf("read_from_LogFile: it=%d quadrant2=%s\n", it, &quadrant2[it * 16]);
#endif
    nitem2[it]++;
    } 
// Eyepiece ("Oculare"):
// Oculare (mm),20,,,Epoca Besseliana,,"2017,0430",,
  else if((it >= 0) && (strncmp(eyepiece_kd, keywd, 7) == 0)) {
    status = csv_read_string(b_in, 2, buffer);
    sscanf(buffer, "%d", &iw);
    eyep2[it] = iw;
#ifdef DEBUG
    printf("read_from_LogFile: it=%d eyep2=%d\n", it, eyep2[it]);
#endif
    nitem2[it]++;
// Epoca Besseliana :
    status = csv_read_string(b_in, 5, buffer);
    if(strncmp(bessel_kd, buffer, 16) == 0) {
      status = csv_read_string(b_in, 7, buffer);
      convert_coma_to_dot(buffer, buffer2);
      sscanf(buffer2, "%lf", &dw);
      EP_Bessel2[it] = dw;
#ifdef DEBUG
    printf("read_from_LogFile: buffer=%s ep_bessel2=%f\n", buffer, EP_Bessel2[it]);
#endif
      }
    } 
  } /* EOF if fgets() */
} /* EOF while loop */
*nitems2 = it+1;
printf("read_from_LogFile/Number of lines: nlines=%d\n", iline);
printf("read_from_LogFile/Number of observations: nitems=%d\n", *nitems2);

// Check nitem2 is good for all observations:
for(i = 0; i < *nitems2; i++) { 
   fprintf(fp_out, "i=%d dble_name: >%s<\n", i, &dble_name2[64 * i]);
   fprintf(fp_out, "date %d-%d-%d\n", day2[i], month2[i], year2[i]);
   fprintf(fp_out, "filter: >%s<\n", &filter2[16 * i]);
   fprintf(fp_out, "eyepiece:%d\n", eyep2[i]);
   fprintf(fp_out, "quadrant: >%s<\n\n", &quadrant2[16 * i]);
   if(nitem2[i] != 5) printf("ERROR for i=%d nitems=%d\n", i, nitem2[i]); 
   }

fclose(fp_in);
fclose(fp_out);
return(status);
}

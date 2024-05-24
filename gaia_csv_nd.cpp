/*************************************************************************
* gaia_csv_nd
*
* JLP
* Version 14/01/2023
*************************************************************************/
#include "tex_calib_utils.h" 
#include "csv_utils.h" 
#include "latex_utils.h" 
#include "astrom_utils1.h" 
#include "astrom_utils2.h" 
#include "jlp_string.h"  // jlp_compact_string
#include "csv_utils.h"  // csv_read_string()

#define SQUARE(a) ((a)*(a))

static int gaia_csv_nd0(char *in_csv, char *out_tex);

#define MAX_LENGTH 256
#define NSTARS_MAX 32

#define DEBUG
/*
*/

static int process_gaia_data(char *wds_name0, char *gaia_source_id, 
                             double *ra_deg, double *dec_deg, 
                             double *paral_mas, double *gmag, int nstars,
                             FILE *fp_out);

/*********************************************************************/
int main(int argc, char *argv[])
{
int status;
char out_tex[64], in_csv[64];

/* If command line with "runs" */
if(argc == 7){
 if(*argv[5]) argc = 6;
 else if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
 else argc = 1;
 }

if(argc != 3)
  {
  printf(" Syntax: gaia_csv_nd in_csv_fname out_tex_fname \n");
  printf("argc=%d\n", argc);
  exit(-1);
  }
else
  {
  strcpy(in_csv,argv[1]);
  strcpy(out_tex,argv[2]);
  }

printf(" OK: in_csv=%s out_tex=%s \n", in_csv, out_tex);

// Scan the file and analyse it 
  status = gaia_csv_nd0(in_csv, out_tex);

return(status);
}
/*************************************************************************
* Scan the file and analyse it 
*
*
*************************************************************************/
static int gaia_csv_nd0(char *in_csv, char *out_tex) 
{
char in_line[MAX_LENGTH], wds_name0[64], wds_name1[64]; 
char gaia_source_id[64*NSTARS_MAX];
double ra_deg[NSTARS_MAX], dec_deg[NSTARS_MAX], paral_mas[NSTARS_MAX];
double gmag[NSTARS_MAX];
int iline, icol_wds, icol_ra, icol_dec, icol_paral, icol_gmag;
int jj, nstars, icol_id;
FILE *fp_in, *fp_out;

icol_id = 1;
icol_ra = 2;
icol_dec = 3;
icol_paral = 4;
icol_gmag = 5;
icol_wds = 6;

if((fp_in = fopen(in_csv, "r")) == NULL) {
  fprintf(stderr, "Fatal error opening %s \n", in_csv);
  return(-1);
  }

if((fp_out = fopen(out_tex, "w")) == NULL) {
  fprintf(stderr, "Fatal error opening %s \n", out_tex);
  return(-1);
  }

// int csv_read_string(char *b_data, int i_column, char *out_string);
// int csv_read_dvalue(char *b_data, int i_column, double *dvalue);

iline = 0;
jj = 0;
while(!feof(fp_in)) {
   if(fgets(in_line, MAX_LENGTH, fp_in)) {
     iline++;
printf("UU iline=%d in_line=%s", iline, in_line);
     if(jj == 0) {
        csv_read_string(in_line, icol_wds, wds_name0);
        } else {
        csv_read_string(in_line, icol_wds, wds_name1);
        }
     csv_read_string(in_line, icol_id, &gaia_source_id[jj*64]);
     csv_read_dvalue(in_line, icol_ra, &ra_deg[jj]);
     csv_read_dvalue(in_line, icol_dec, &dec_deg[jj]);
     csv_read_dvalue(in_line, icol_paral, &paral_mas[jj]);
     csv_read_dvalue(in_line, icol_gmag, &gmag[jj]);
     printf("UV wds_name0=%s wds_name1=%s jj=%d\n", wds_name0, wds_name1, jj);
     if(jj > 0) {
       if(strcmp(wds_name0, wds_name1) != 0) {
 printf("\n");
         nstars = jj;
         if(nstars > 1) {
           process_gaia_data(wds_name0, gaia_source_id, ra_deg, dec_deg, 
                             paral_mas, gmag, nstars, fp_out);
           } else {
           fprintf(fp_out, "%s mag=%.2f one star only\n", wds_name0, gmag[jj]); 
           }
// Initialise (values for jj=0) with last values;
         ra_deg[0] = ra_deg[jj];
         dec_deg[0] = dec_deg[jj];
         paral_mas[0] = paral_mas[jj];
         gmag[0] = gmag[jj];
         strcpy(wds_name0, wds_name1);
         strcpy(&gaia_source_id[0], &gaia_source_id[jj*64]);
         jj = 0;
         }
     }
     if(jj < NSTARS_MAX - 1) jj++;
    }
 }


fclose(fp_in);
fclose(fp_out);

return(0);
}
/***************************************************************************
*
***************************************************************************/
static int process_gaia_data(char *wds_name0, char *gaia_source_id, 
                             double *ra_deg, double *dec_deg, 
                             double *paral_mas, double *gmag, int nstars,
                             FILE *fp_out)
{
double gmag0, cos_dec1, dra, ddec;
double rho[NSTARS_MAX], theta[NSTARS_MAX];
int i, i_star0;

  printf("MMwds_name=%s nstars=%d\n", wds_name0, nstars);

// Look for the brightest star:
  gmag0 = 100;
  for(i = 0; i < nstars; i++) {
    if(gmag[i] < gmag0) {
      gmag0 = gmag[i];
      i_star0 = i;
      } 
    }

// Compute rho and theta for all stars:
  for(i = 0; i < nstars; i++) {
     dra = ra_deg[i] - ra_deg[i_star0];
     ddec = dec_deg[i] - dec_deg[i_star0];
     cos_dec1 = cos(dec_deg[i_star0] * (3.14159 / 180.));
     dra *= cos_dec1;
// Conversion to arcseconds:
     dra *= (60. * 60.);
     ddec *= (60. * 60.);

// rho:
     rho[i] = SQUARE(dra) + SQUARE(ddec);
     if(rho[i] > 0) rho[i] = sqrt(rho[i]);
      else rho[i] = 0.; 

// theta is measured from the north direction and positive to the east
     if(rho[i] <= 0.) {
       theta[i] = 0.;
     } else {
       theta[i] = acos((ddec / rho[i]));
// Conversion to degrees:
       theta[i] *= (180. / 3.14159);
       if(dra < 0.) theta[i] += 180.;
     }
//       if(theta[i] < 0.) theta[i] += 180.;
   }

// Comments
  for(i = 0; (i < nstars); i++) {
    if(i == i_star0) {
       printf("%s %d %s mag=%.2f paral=%.2f (primary) \n", 
               wds_name0, i + 1, &gaia_source_id[i*64], gmag[i], paral_mas[i]);
       fprintf(fp_out, "%s %d %s mag=%.2f paral=%.2f (primary) \n", 
               wds_name0, i + 1, &gaia_source_id[i*64], gmag[i], paral_mas[i]);
       } else {
       printf("%s %d %s mag=%.2f rho=%.3f theta=%.1f paral=%.2f\n", 
               wds_name0, i + 1, &gaia_source_id[i*64], gmag[i], rho[i], theta[i], 
               paral_mas[i]);
       fprintf(fp_out, "%s %d %s mag=%.2f rho=%.3f theta=%.1f paral=%.2f\n", 
               wds_name0, i + 1, &gaia_source_id[i*64], gmag[i], rho[i], theta[i], 
               paral_mas[i]);
       }
    }
       fprintf(fp_out, "\n"); 

return(0);
}

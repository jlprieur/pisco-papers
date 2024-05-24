/*************************************************************************
* Program stat_ascii
* To compute statistics on a ASCII list 
*
*
* JLP
* Version 06/09/2020
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h> // isdigit()
#include <math.h>
#include <string.h>
#include "csv_utils.h" // blank_read_dvalue(); 

#define MINI(a,b) ((a) < (b)) ? (a) : (b)
#define MAXI(a,b) ((a) < (b)) ? (b) : (a)

#define DEBUG
/*
*/

static int read_list_from_file(char *filein, int icol_data, int icol_weights,
                               double *dvals, double *dweights, int *nvals, 
                               int nmaxi);
static int read_file_to_list(char *filein, int icol_data, int icol_weights,
                             double *dvals, double *dweights, int *nvals, 
                             int nmaxi);
static int compute_stats_from_list(double *dvals, double *dweights, int nvals,
                                   double lowcut, double highcut, 
                                   double *mean, double *sigma,
                                   double *vmini, double *vmaxi);

int main(int argc, char *argv[])
{
char filein[128];
double mean, sigma, *dvals, *dweights;
double lowcut, highcut, vmini, vmaxi;
int icol_data, icol_weights = 0, nvals, nmaxi;
register int i;

  printf("stat_ascii/ JLP/ Version 29/08/2023\n");

if(argc == 7 && *argv[4]) argc = 5;
if(argc == 7 && *argv[3]) argc = 4;
if(argc == 7 && *argv[2]) argc = 3;
if(argc == 7 && *argv[1]) argc = 2;
if(argc != 3 && argc != 4)
  {
  printf("Error/Bad syntax: argc=%d\n\n", argc);
  printf("Syntax:        stat_ascii in_file icol_data [icol_weights,nmaxi]\n");
  exit(-1);
  }
else
  {
  strcpy(filein,argv[1]);
  sscanf(argv[2],"%d", &icol_data);
   if(argc > 3) {
     sscanf(argv[3],"%d,%d", &icol_weights, &nmaxi);
     } else {
     icol_weights = 0;
     nmaxi = 2048;
     }
  }

printf(" OK: filein=%s icol_data=%d nmaxi=%d\n", filein, icol_data, nmaxi);

dvals = new double[nmaxi];
dweights = new double[nmaxi];

/* Scan the file and build the lists: */
  read_list_from_file(filein, icol_data, icol_weights, dvals, dweights, &nvals, 
                      nmaxi);

  lowcut = -1.e+12;
  highcut = -lowcut;
  compute_stats_from_list(dvals, dweights, nvals, lowcut, highcut, 
                          &mean, &sigma, &vmini, &vmaxi);

/* Series of iterations with rejection of outliers: */
for (i = 0; i < 3; i++) {
  lowcut = mean - (1.8 * sigma);
  highcut = mean + (1.8 * sigma);
  compute_stats_from_list(dvals, dweights, nvals, lowcut, highcut, 
                          &mean, &sigma, &vmini, &vmaxi);
  }

delete dvals;
delete dweights;
return(0);
}
/*************************************************************************
*
* INPUT:
*  icol_data : column number of data (in the range [1, 8])
*  icol_weights : column number of weights (in the range [1, 8])
*************************************************************************/
static int read_list_from_file(char *filein, int icol_data, int icol_weights,
                               double *dvals, double *dweights, int *nvals, 
                               int nmaxi)
{
int i, nn, stat1, stat2;
double val, weight;
char buffer[256];
FILE *fp_in;

if((fp_in = fopen(filein,"r")) == NULL)
{
printf(" Fatal error opening input file %s \n",filein);
return(-1);
}

nn = 0;
i = 0;
while(1) {

/* Read new line */
  if(fgets(buffer,256,fp_in) == NULL) break;

  val = 0.;
  weight = 1.; 
  if(buffer[0] != '#' && buffer[0] != '%') {
    i++;
    stat1 = blank_read_dvalue(buffer, icol_data, &val);
    if(icol_weights != 0) 
      stat2 = blank_read_dvalue(buffer, icol_weights, &weight);
    else 
      stat2 = 0;

    if((stat1 == 0) && (stat2 == 0)) {
      dvals[nn] = val;
      dweights[nn] = weight; 
      nn++;
      if(nn >= nmaxi) {
        fprintf(stderr, "Fatal error nn > nmaxi=%d\n", nmaxi);
        exit(-1);
        }

#ifdef DEBUG
  printf(" Buffer=%s\n", buffer);
  printf(" icol_data=%d icol_weights=%d nn=%d val=%f weight=%f\n",
         icol_data, icol_weights, nn, val, weight);
#endif
      }
  }
}

*nvals = nn;

fclose(fp_in);
return(0);
}
/*************************************************************************
*
* INPUT:
*  icol_data : column number of data (in the range [1, 8])
*  icol_weights : column number of weights (in the range [1, 8])
*************************************************************************/
static int read_file_to_list(char *filein, int icol_data, int icol_weights,
                             double *dvals, double *dweights, int *nvals, 
                             int nmaxi)
{
int ival, iweight, i, nn;
float val, weight;
char buffer[256];
FILE *fp_in;

if((fp_in = fopen(filein,"r")) == NULL)
{
printf(" Fatal error opening input file %s \n",filein);
return(-1);
}

if(icol_data < 1 || icol_data > 13) {
  printf("compute_stats_from_list/Only icol_data = 1 to 13 are allowed yet\n");
  fclose(fp_in);
  return(-1);
  }

if(icol_weights < 0 || icol_weights > 13) {
  printf("compute_stats_from_list/Only icol_weights = 1 to 13 are allowed yet\n");
  fclose(fp_in);
  return(-1);
  }

nn = 0;
i = 0;
while(1) {

/* Read new line */
  if(fgets(buffer,256,fp_in) == NULL) break;

  ival = 0;
  val = 0.;
  weight = 1.; 

// if(buffer[0] != '#' && buffer[0] != '%') printf("ZZZ %s\n", buffer);

if(buffer[0] != '#' && buffer[0] != '%'
   && isdigit(buffer[0])) {
  i++;

  switch(icol_data) { 
    case 1:
      ival = sscanf(buffer, "%f", &val);
      break;
    case 2:
      ival = sscanf(buffer, "%*f %f", &val);
      break;
    case 3:
      ival = sscanf(buffer, "%*f %*f %f", &val);
      break;
    case 4:
      ival = sscanf(buffer, "%*f %*f %*f %f", &val);
      break;
    case 5:
      ival = sscanf(buffer, "%*f %*f %*f %*f %f", &val);
      break;
    case 6:
      ival = sscanf(buffer, "%*f %*f %*f %*f %*f %f", &val);
      break;
    case 7:
      ival = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %f", &val);
      break;
    case 8:
      ival = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %*f %f", &val);
      break;
    case 9:
      ival = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %*f %*f %f", &val);
      break;
    case 10:
      ival = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %*f %*f %*f %f", &val);
      break;
    case 11:
      ival = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f", &val);
      break;
    case 12:
      ival = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f", &val);
      break;
    case 13:
      ival = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f", &val);
      break;
    default:
      fprintf(stderr,"compute_stats_from_list/Fatal error, wrong column (icol_data=%d)\n", icol_data);
    exit(-1);
    }

  switch(icol_weights) { 
    case 1:
      iweight = sscanf(buffer, "%f", &weight);
      break;
    case 2:
      iweight = sscanf(buffer, "%*f %f", &weight);
      break;
    case 3:
      iweight = sscanf(buffer, "%*f %*f %f", &weight);
      break;
    case 4:
      iweight = sscanf(buffer, "%*f %*f %*f %f", &weight);
      break;
    case 5:
      iweight = sscanf(buffer, "%*f %*f %*f %*f %f", &weight);
      break;
    case 6:
      iweight = sscanf(buffer, "%*f %*f %*f %*f %*f %f", &weight);
      break;
    case 7:
      iweight = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %f", &weight);
      break;
    case 8:
      iweight = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %*f %f", &weight);
      break;
    case 9:
      iweight = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %*f %*f %f", &weight);
      break;
    case 10:
      iweight = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %*f %*f %*f %f", &weight);
      break;
    case 11:
      iweight = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f", &weight);
      break;
    case 12:
      iweight = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f", &weight);
      break;
    case 13:
      iweight = sscanf(buffer, "%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f", &weight);
      break;
    default:
      weight = 1.;
      iweight = 1;
      break;
    }

    if((ival != 0) && (iweight != 0)) {
      dvals[nn] = val;
      dweights[nn] = weight; 
      nn++;
      if(nn >= nmaxi) {
        fprintf(stderr, "Fatal error nn > nmaxi=%d\n", nmaxi);
        exit(-1);
        }

#ifdef DEBUG
  printf(" Buffer=%s\n", buffer);
  printf(" icol_data=%d icol_weights=%d ival=%d iweight=%d nn=%d val=%f weight=%f\n",
         icol_data, icol_weights, ival, iweight, nn, val, weight);
#endif
      }
 }
}

*nvals = nn;

fclose(fp_in);
return(0);
}
/*************************************************************************
*
* INPUT:
*  icol_data : column number of data (in the range [1, 8])
*  icol_weights : column number of weights (in the range [1, 8])
*************************************************************************/
static int compute_stats_from_list(double *dvals, double *dweights, int nvals,
                                   double lowcut, double highcut, 
                                   double *mean, double *sigma,
                                   double *vmini, double *vmaxi)
{
double sum, sumsq, sumw, val, weight, bad_value;
int nn, i;

bad_value = -12345.;

#ifdef DEBUG
printf("DEBUG 2020: bad_value=%f nvals=%d\n", bad_value, nvals);
#endif

*vmini = 1.e+12;
*vmaxi = -1.e+12;

sumw = 0.;
sum = 0.; sumsq = 0.;
nn = 0;
for(i = 0; i < nvals; i++) {
   val = dvals[i];
   weight = dweights[i];
    if((val != bad_value) && (val > lowcut) && (val < highcut)) {
     *vmini = MINI(val, *vmini);
     *vmaxi = MAXI(val, *vmaxi);
     sumw += weight;
     sum += val * weight;
     sumsq += val * val * weight;
     nn++;
     }
  }

if(nn > 3) {
  sum = sum / sumw; 
  *mean = sum;
  sumsq = sumsq / sumw - sum * sum;
  *sigma = sqrt(sumsq); 
  printf("compute_stats_from_list/OK: nvalues=%d mean=%f sigma=%f (min=%f max=%f) (lowcut=%f highcut=%f)\n", 
       nn, *mean, *sigma, *vmini, *vmaxi, lowcut, highcut);
  }
else {
  fprintf(stderr," Error: too few points for computing statistics, nvalues= %d\n",
          nn);
  *sigma = 0.;
  *mean = 0.;
  }

return(0);
}

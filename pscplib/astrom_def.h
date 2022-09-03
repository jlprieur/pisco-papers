/*************************************************************************
*
* JLP
* Version 03/03/2020
*************************************************************************/
#ifndef _astrom_def
#define _astrom_def

/* Maximum length for a line will be 1024 characters: */
#define NMAX 1024 

/* Positive number (to be able to compute sqrt in compute_mean) */
#define NO_DATA 10000.
#ifndef ABS   /* Defined in jlp_catalog_utils.h */
#define ABS(x) ((x) >= 0. ? (x) : (-(x)))
#endif
 
/*
#define DEBUG
*/

/* Maximum number of measurements per object: */
#define NMEAS 64
/* Maximum number of objects: */
#define NOBJ_MAX 8000

/* Structure to define a measurement: */
typedef struct {
char comments[40][80]; 	/* Comments (lines starting with %) */
char data[180];         /* Data: line with filename, filter, measurements */
char filename[40];      /* Name of the FITS file used for this measurement */
char date[12];          /* Date, e.g., 29/12/2004 */
double epoch;            /* Julian day as a fraction of year, e;g., 2004.23 */
char filter[10];        /* Filter name: V, R, Sf */
int eyepiece;           /* Focal length of the magnifying eyepiece */
int is_new_double;      /* flag set to one if nd or ND in comments */
int quadrant;           /* Quadrant value */
int dquadrant;          /* Quadrant uncertainty (0 or 1) */
double rho;              /* Angular separation of the binary (arcseconds) */
double drho;             /* Error on rho (arcseconds) */
double theta;            /* Position angle of the companion (relative to North)*/
double dtheta; 		/* Error on theta (degrees) */
double dmag;            /* Delta mag */
double ddmag;           /* Error on Delta mag */
char notes[80];         /* Notes in the last column of the data line */
int from_recorded_file; /* Flag set to 1 if measure from recorded file (in SVHS tape) */
int flagged_out;        /* Flag to cancel output (for publication) */
} MEASURE; 

/* Structure to define an object */
typedef struct {
char wds[40];		/* WDS name */
char discov_name[40];	/* Discoverer's binary name */
char comp_name[40];	/* Companion name */
char ads[40];		/* ADS name */
MEASURE meas[NMEAS];	/* Measurements concerning this object */
MEASURE cmp_meas[NMEAS]; /* Measurements in compared file concerning this object */
char notes[80];		/* Notes which are common to all measurements */
double ra; 		/* Right ascension */
int dec; 		/* Declination */
int nmeas;              /* Nber of measurements for this object*/
double WR;               /* Radius value of last measurement in WDS_CHARA data base */
double WT;               /* Theta value of last measurement in WDS_CHARA data base */
double WY;               /* Year of last measurement in WDS_CHARA data base */
} OBJECT;

#endif

/*************************************************************************
* orbit_plot1
* To plot the curves rho(epoch), theta(epoch) and XY orbits in the sky plane
* measurements only (and/without the orbit).
* 
*
* JLP
* Version 12/01/2012
**************************************************************************/
#include <orbit_plot_utils.h>

/****************** Main program **************************************/
int main(int argc, char *argv[])
{
int status, npts_max, iformat, iplot, smoothed_values, nber_of_orbits;
int resid_vectors;
char measures_infile[64], plotfile[64], comments[81];
char orbit_file[64], plot_title[64];

 printf("Program orbit_plot1 to plot raw data (measurements with/without orbit)\n");
status = -1;
// General case (more than 6 arguments): assume that orbit parameters have been entered:
if (argc >= 6) {
 strcpy(measures_infile, argv[1]);
 strcpy(plotfile, argv[2]);
 sscanf(argv[3], "%d,%d,%d", &iformat, &iplot, &resid_vectors);
 strcpy(orbit_file, argv[4]);
 strcpy(plot_title, argv[5]);
 if (argc == 7) { 
   sscanf(argv[6], "%d", &nber_of_orbits);
   if(nber_of_orbits < 0 || nber_of_orbits > 3) {
    fprintf(stderr, "Fatal error: nber_of_orbits=%d\n", nber_of_orbits);
    return(-1);
    }
 } else {
   nber_of_orbits = 1;
 }
 if((iformat == 3) || (iformat == 4)) status = 0;
// Case of 5 arguments only: assume that no orbit parameters have been entered:
} else if (argc == 5) {
 strcpy(measures_infile, argv[1]);
 strcpy(plotfile, argv[2]);
 sscanf(argv[3], "%d,%d,%d", &iformat, &iplot, &resid_vectors);
 strcpy(plot_title, argv[4]);
 orbit_file[0] = '\0';
 if((iformat > 0 && iformat < 4) || iformat == 5) status = 0;
 }  

printf("OK: argc=%d infile=%s plot_file=%s orbit_file=%s\n", 
       argc, measures_infile, plotfile, orbit_file);
printf("OK: iformat=%d iplot=%d resid_vectors=%d nber_of_orbits=%d\n", 
       iformat, iplot, resid_vectors, nber_of_orbits);
printf("OK: plot_title=%s\n", plot_title);

// Handle errors: 
if(status) {
 printf("Error, syntax is: \norbit_plot1 infile extension_for_pst_files iformat,iplot,resid_vectors  [file_with_orbital_elements, if iformat==3 or 4] title [nber_of_orbits]\n");
 printf("iformat: 1=raw (from WDS)\n");
 printf("         2=corrected measures (from 1BIN.FOR of orbit_weight.c)\n");
 printf("         3=measures and ephemerids (from residuals_2.c)\n");
 printf("         4=raw or corrected measures and orbit\n");
 printf("         5=residuals (dx,dy) (from residuals_2.c)\n");
 printf("iplot: 1=rho,theta 2=dx,dy 3=XY_full_orbit 4=XY_part_of_orbit\n");
 printf("resid_vectors: 0 or 1 if residual vectors to be plotted (for XY orbit only)\n");
 printf("Example: \norbit_plot1 A15971.205 tt 2,1,0 \"ADS 15971\" \n");
 printf("Example: \norbit_plot1 A15971.205 tt 4,3,0 A15971_orbit.txt \"ADS 15971\"\n");
 printf("Example: \norbit_plot1 A15971.205 tt 4,3,0 A15971_orbits.txt  \"ADS 15971\" 2\n");
 printf("Example: \norbit_plot1 A15971.205 tt 3,4,1 A15971_orbit.txt \"ADS 15971\"\n");
 return(-1);
}
 
/* 
* iformat = 1: epoch, rho_O, theta_O, n_nights, author, aperture, instrument
* iformat = 2: epoch, rho_O, theta_O, n_nights, author, aperture, weight 
* iformat = 3: epoch, rho_O, rho_C, Drho_O-C, theta_O, theta_C, Dtheta, author
* iformat = 4: epoch, rho_O, theta_O, n_nights, author, aperture, ...
* iformat = 5: epoch, Dx, Dy, n_nights, author, aperture, weight 
*/
/* 
* iplot=1: plot rho vs epoch and theta vs epoch 
* iplot=2: plot dx vs epoch and dy vs epoch 
* iplot=3: plot XY_full_orbit 
* iplot=4: plot XY_part_of_orbit
*/
/*
* resid_vectors=0 : do not plot residual vectors 
* resid_vectors=1 : plot residual vectors 
*/
strcpy(comments,"Program orbit_plot1.c  -- Version 23/07/2018");

status = orbit_get_npts_max(measures_infile, &npts_max);
if(status) return(-1);

/* Minimum of 1000 points when iformat=4 (since we want a full orbit) */
if((iformat == 3) || (iformat == 4)) npts_max = MAXI(1000,npts_max);

/* strcpy(plotdev,"square/rho.ps");
*/
/* Only useful for Drho, Dtheta, not relevant for Dx,Dy (i.e., when iformat=5)
*/
status = 0;
switch (iplot) {
 case 1: 
   if(iformat != 5) {
   orbit_plot_rho(measures_infile, comments, orbit_file, nber_of_orbits, 
                  npts_max, iformat, plotfile);
   orbit_plot_theta(measures_infile, comments, orbit_file, nber_of_orbits, 
                    npts_max, iformat, plotfile);
   } else {
   fprintf(stderr, "Fatal error: cannot plot rho,theta with dx,dy file !\n");
   status = -1;
   }
   break;
 case 2: 
   smoothed_values = 1;
   orbit_plot_Dx(measures_infile, comments, npts_max, plotfile, smoothed_values);
   orbit_plot_Dy(measures_infile, comments, npts_max, plotfile, smoothed_values);
   break;
/* Plot orbit in the plane of the sky
* Full orbit if iformat==3 or 4, or part of the orbit only in other cases 
*/
 case 3:
 case 4:
   orbit_plot_skyplane(measures_infile, comments, orbit_file, nber_of_orbits, 
                       npts_max, iformat, iplot, resid_vectors, plotfile, 
                       plot_title);
   break;
 default:
   fprintf(stderr, "Fatal error: unkown option: iplot=%d !\n", iplot);
   status = -1;
   break;
} // EOF switch

return(status);
}

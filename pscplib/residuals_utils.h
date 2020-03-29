/************************************************************************
* "residuals_utils.h"
*
* JLP 
* Version 06/07/2009
*************************************************************************/
#ifndef _residuals_utils_h /* BEOF sentry */
#define _residuals_utils_h
#include "jlp_catalog_utils.h"

/* Declaring linkage specification to have "correct names"
* that can be linked with C programs */

#ifdef __cplusplus
extern "C" {
#endif

int compute_ephemerid_of_multiple_system(int nber_of_orbits,
                      double *Omega_node, double *omega_peri, double *i_incl, 
                      double *e_eccent, double *T_periastron, double *Period, 
                      double *a_smaxis, double *mean_motion, double epoch_o, 
                      double *c_tolerance, double *theta_c, double *rho_c);
int compute_ephemerid(double Omega_node, double omega_peri, double i_incl,
                      double e_eccent, double T_periastron, double Period, 
                      double a_smaxis, double mean_motion, double epoch_o, 
                      double c_tolerance, double *theta_c, double *rho_c);
int precession_correction(double *dtheta_precess, double alpha, double delta, 
                          double epoch_o, double orbit_equinox);
int read_orbital_elements_from_file(char *orbit_infile, int iformat,
              double *Omega_node, double *omega_peri, double *i_incl,
              double *e_eccent, double *T_periastron, double *Period,
              double *a_smaxis, double *mean_motion, double *orbit_equinox,
              int nber_of_orbits);
int read_orbital_elements_and_errors_from_file(char *orbit_infile,
              int orbit_format,
              double *Omega_node, double *omega_peri, double *i_incl,
              double *e_eccent, double *T_periastron, double *Period,
              double *a_smaxis, double *mean_motion, double *orbit_equinox,
              double *err_Omega_node, double *err_omega_peri, 
              double *err_i_incl, double *err_e_eccent, 
              double *err_T_periastron, double *err_Period,
              double *err_a_smaxis, double *err_mean_motion, 
              int nber_of_orbits);
int compute_Thiele_elements(double Omega_node, double omega_peri, double i_incl,
              double e_eccent, double T_periastron, double Period,
              double a_smaxis, double mean_motion, double orbit_equinox,
              double *AA, double *BB, double *FF, double *GG);

#ifdef __cplusplus
}
#endif

#endif /* EOF sentry */

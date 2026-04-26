/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_Ham,         // Hamiltonian constraint
    c_Mom1,         // Momentum constraints
    c_Mom2,
    c_Mom3,
    c_Ham_abs_sum, // Sum of absolute value of each term in Ham
    c_Mom_abs_sum1, // Sum of absolute value of each term in Mom
    c_Mom_abs_sum2,
    c_Mom_abs_sum3,
   
    c_sqrt_gamma,  // \sqrt(\gamma) = pow(chi,-3/2) volume factor of spatial
                   // metric
    
    c_rho_scaled,  // \rho * volume factor
    c_S_scaled,    // S * volume factor
    c_K_scaled,    // K * volume factor
    c_rho_phi,
    c_rho_g2,
    c_rho_g3,
    c_rho_GB,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",      "Mom1",     "Mom2",
    "Mom3",
    "Ham_abs",    "Mom_abs1","Mom_abs2","Mom_abs3", "sqrt_gamma",
    "rho_scaled", "S_scaled", "K_scaled",
    

    "rho_phi",  "rho_g2",   "rho_g3",
    "rho_GB"

    };
}

#endif /* DIAGNOSTICVARIABLES_HPP */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(COSMODIAGNOSTICS_HPP_)
#error "This file should only be included through CosmoDiagnostics.hpp"
#endif

#ifndef COSMODIAGNOSTICS_IMPL_HPP_
#define COSMODIAGNOSTICS_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class theory_t>
template <class data_t>
void CosmoDiagnostics<theory_t>::compute(Cell<data_t> current_cell) const
{
    // Load local vars and calculate derivs
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Vars>(current_cell);

    // Inverse metric and Christoffel symbol
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    // Define quantities
    data_t rho;
    data_t sqrt_gamma;
    data_t S;
    data_t rho_scaled;
    data_t S_scaled;
    data_t K_scaled;
    data_t a;
    const data_t ln_chi       = log(vars.chi);
    const data_t ln_chi_bg    = log(m_params.chi_background);  
    const data_t delta_ln_chi = ln_chi - ln_chi_bg;

    const data_t zeta         = -0.5 * delta_ln_chi;

    Tensor<1, data_t> d_ln_gamma;
    FOR(i) { d_ln_gamma[i] = -3.0 * d1.chi[i] / vars.chi; }
    const data_t d0_ln_gamma = -2.0 * vars.K;
    const data_t d0_phi = vars.lapse * vars.Pi;
    Tensor<1, data_t> d_phi;
    FOR(i) { d_phi[i] = d1.phi[i]; }
    data_t R_sq = 0.0;
    Tensor<1, data_t> R_vec;
    FOR(i)
    {
        R_vec[i] = -d_ln_gamma[i] / 6.0
                 + (d0_ln_gamma / (d0_phi + 1e-30)) * d_phi[i] / 6.0;
        R_sq += R_vec[i] * R_vec[i];
    }
    const data_t R_mag = sqrt(R_sq);
    const data_t zeta2 = zeta * zeta;
    const data_t zeta3 = zeta * zeta * zeta;
    const data_t zeta_w  = zeta  * sqrt_gamma;
    const data_t zeta2_w = zeta2 * sqrt_gamma;
    const data_t zeta3_w = zeta3 * sqrt_gamma;
    // Energy Momentum Tensor
    const auto rho_and_Si =
        theory_t.compute_rho_and_Si(vars, d1, d2, coords);
    const auto Sij_TF_and_S =
        theory_t.compute_Sij_TF_and_S(vars, d1, d2, vars,coords);
    const auto all_rhos =
        theory_t.compute_all_rhos(vars, d1, d2, coords);

    const data_t sqrt_gamma = pow(vars.chi, -3. / 2.);
    const data_t rho       = rho_and_Si.rho;
    const data_t S         = Sij_TF_and_S.S;
    const data_t K_scaled = vars.K / pow(vars.chi, 3. / 2.);
    const data_t rho_scaled = rho / pow(vars.chi, 3. / 2.);
    const data_t S_scaled = S / pow(vars.chi, 3. / 2.);
    const data_t rho_phi   = all_rhos.phi;
    const data_t rho_g2    = all_rhos.g2;
    const data_t rho_g3    = all_rhos.g3;
    const data_t rho_phi_scaled = rho_phi / pow(vars.chi, 3. / 2.);
    const data_t rho_g2_scaled  = rho_g2 / pow(vars.chi, 3. / 2.);
    const data_t rho_g3_scaled  = rho_g3 / pow(vars.chi, 3. / 2.);
    // Write the diagnostics into the output vector of cells
    current_cell.store_vars(sqrt_gamma, c_sqrt_gamma);
    current_cell.store_vars(rho, c_rho);
    current_cell.store_vars(rho_scaled, c_rho_scaled);
    current_cell.store_vars(S_scaled, c_S_scaled);
    current_cell.store_vars(K_scaled, c_K_scaled);
    current_cell.store_vars(rho_phi,       c_rho_phi);
    current_cell.store_vars(rho_phi_scaled,c_rho_phi_scaled);
    current_cell.store_vars(rho_g2,        c_rho_g2);
    current_cell.store_vars(rho_g2_scaled, c_rho_g2_scaled);
    current_cell.store_vars(rho_g3,        c_rho_g3);
    current_cell.store_vars(rho_g3_scaled, c_rho_g3_scaled);
    current_cell.store_vars(zeta,            c_zeta);
    current_cell.store_vars(zeta_w,            c_zeta_w);
    current_cell.store_vars(zeta2_w,           c_zeta2_w);
    current_cell.store_vars(zeta3_w,           c_zeta3_w);
    current_cell.store_vars(R_mag,           c_R_mag);
}

#endif /* COSMODIAGNOSTICS_IMPL_HPP_ */

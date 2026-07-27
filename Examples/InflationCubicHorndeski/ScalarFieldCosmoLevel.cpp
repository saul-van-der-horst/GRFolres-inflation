/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "ScalarFieldCosmoLevel.hpp"
#include "AMRReductions.hpp"

#include "BoxLoops.hpp"
#include "CosmoHamTaggingCriterion.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ComputePack.hpp"

#include "ModifiedCCZ4RHS.hpp"
#include "ModifiedGravityConstraints.hpp"
#include "ModifiedGravityWeyl4.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "CosmoDiagnostics.hpp"
#include "RhoDiagnostics.hpp"
#include "CustomExtraction.hpp"
#include "SetValue.hpp"
#include "SixthOrderDerivatives.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"
#include "CubicHorndeski.hpp"
#include "GammaCalculator.hpp"
#include "CouplingAndPotential.hpp"

#include "CosmoModifiedPunctureGauge.hpp"
// Things to do during the advance step after RK4 steps
void CosmoLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck("NaNCheck in specific Advance: "), m_state_new,
                       m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// This initial data uses an approximation for the metric which
// is valid for small boosts
void CosmoLevel::initialData()
{
    CH_TIME("CosmoLevel::initialData");
    if (m_verbosity)
        pout() << "CosmoLevel::initialData " << m_level << endl;

    // First set everything to zero then initial conditions for scalar field -
    // Set initial condition of inflaton, see details in Potential.hpp and
    // InitialScalarData.hpp
   
   
    fillAllGhosts();
    // Note that the GammaCaluculator is not necessary since the data is
    // conformally flat. It is left here for generality.
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);

    
    CouplingAndPotential coupling_and_potential(
            m_p.coupling_and_potential_params);
    CubicHorndeskiWithCouplingAndPotential cubic_horndeski(
            coupling_and_potential);    
    // Calculate constraints and some diagnostics as we need it in tagging
    // criterion
    BoxLoops::loop(
     ModifiedGravityConstraints<CubicHorndeskiWithCouplingAndPotential>(
         cubic_horndeski, m_dx, m_p.center, m_p.G_Newton, c_Ham,
         Interval(c_Mom1, c_Mom3), c_Ham_abs_sum,
         Interval(c_Mom_abs_sum1, c_Mom_abs_sum3)),
     m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    CosmoDiagnostics<CubicHorndeskiWithCouplingAndPotential> cosmo_diagnostics(
        cubic_horndeski, m_dx, m_p.G_Newton);
    BoxLoops::loop(cosmo_diagnostics, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
    
}

// Calculate RHS during RK4 substeps
void CosmoLevel::specificEvalRHS(GRLevelData &a_soln,
                                                  GRLevelData &a_rhs,
                                                  const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate ModifiedCCZ4 right hand side with theory_t = CubicHorndeski
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    CubicHorndeskiWithCouplingAndPotential cubic_horndeski(
        coupling_and_potential);
    CosmoModifiedPunctureGauge cosmo_modified_puncture_gauge(m_p.modified_ccz4_params);
    cosmo_modified_puncture_gauge.set_K_mean(m_cosmo_amr.get_K_mean());
    if (m_p.max_spatial_derivative_order == 4)
    {
        ModifiedCCZ4RHS<CubicHorndeskiWithCouplingAndPotential,
                        CosmoModifiedPunctureGauge, FourthOrderDerivatives>
            my_ccz4_theory(cubic_horndeski, m_p.modified_ccz4_params,
                           cosmo_modified_puncture_gauge, m_dx, m_p.sigma, m_p.center,
                           m_p.G_Newton);

        BoxLoops::loop(my_ccz4_theory, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        ModifiedCCZ4RHS<CubicHorndeskiWithCouplingAndPotential,
                        CosmoModifiedPunctureGauge, SixthOrderDerivatives>
            my_ccz4_theory(cubic_horndeski, m_p.modified_ccz4_params,
                           cosmo_modified_puncture_gauge, m_dx, m_p.sigma, m_p.center,
                           m_p.G_Newton);

        BoxLoops::loop(my_ccz4_theory, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// enforce trace removal during RK4 substeps
void CosmoLevel::specificUpdateODE(GRLevelData &a_soln,
                                                    const GRLevelData &a_rhs,
                                                    Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void CosmoLevel::preTagCells()
{
     fillAllGhosts();
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    CubicHorndeskiWithCouplingAndPotential cubic_horndeski(
        coupling_and_potential);
    BoxLoops::loop(
    ModifiedGravityConstraints<CubicHorndeskiWithCouplingAndPotential>(
        cubic_horndeski, m_dx, m_p.center, m_p.G_Newton, c_Ham,
        Interval(c_Mom1, c_Mom3), c_Ham_abs_sum,
        Interval(c_Mom_abs_sum1, c_Mom_abs_sum3)),
    m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    CosmoDiagnostics<CubicHorndeskiWithCouplingAndPotential> cosmo_diagnostics(
        cubic_horndeski, m_dx, m_p.G_Newton);
    BoxLoops::loop(cosmo_diagnostics, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
}

// specify the cells to tag
void CosmoLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)

{
    double rho_mean = m_cosmo_amr.get_rho_mean();
    BoxLoops::loop(CosmoHamTaggingCriterion(m_dx, m_p.tagging_center,
                                            m_p.tagging_radius, rho_mean),
                   current_state_diagnostics, tagging_criterion);
}

void CosmoLevel::specificPostTimeStep()
{
    int min_level = m_p.extraction_params.min_extraction_level();
    bool calculate_diagnostics = at_level_timestep_multiple(min_level);
    bool first_step = (m_time == 0.);
    CH_TIME("InflationCubicHorndeskiLevel::specificPostTimeStep");

    // No need to evaluate the diagnostics more frequently than every coarse
    // timestep, but must happen on every level (not just level zero or data
    // will not be populated on finer levels)

    if (calculate_diagnostics)
    {
        fillAllGhosts();
        CouplingAndPotential coupling_and_potential(
         m_p.coupling_and_potential_params);
        CubicHorndeskiWithCouplingAndPotential cubic_horndeski(
         coupling_and_potential);
        BoxLoops::loop(ModifiedGravityConstraints<CubicHorndeskiWithCouplingAndPotential>(
            cubic_horndeski, m_dx, m_p.center, m_p.G_Newton, c_Ham,
            Interval(c_Mom1, c_Mom3), c_Ham_abs_sum,
            Interval(c_Mom_abs_sum1, c_Mom_abs_sum3)),
    m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        CosmoDiagnostics<CubicHorndeskiWithCouplingAndPotential> cosmo_diagnostics(
            cubic_horndeski, m_dx, m_p.G_Newton);
        BoxLoops::loop(cosmo_diagnostics, m_state_new, m_state_diagnostics,
                       EXCLUDE_GHOST_CELLS);

        if (m_level == min_level)
        {
            // AMRReductions for diagnostic variables
            AMRReductions<VariableType::diagnostic> amr_reductions_diagnostic(
                m_cosmo_amr);
            double phys_vol = amr_reductions_diagnostic.sum(c_sqrt_gamma);
            double L2_Ham = amr_reductions_diagnostic.norm(c_Ham);
            double L2_Mom = amr_reductions_diagnostic.norm(Interval(c_Mom1, c_Mom3));
            double K_total = amr_reductions_diagnostic.sum(c_K_scaled);
            m_cosmo_amr.set_rho_mean(
                amr_reductions_diagnostic.sum(c_rho_scaled) / phys_vol);
            m_cosmo_amr.set_S_mean(amr_reductions_diagnostic.sum(c_S_scaled) /
                                   phys_vol);
            m_cosmo_amr.set_K_mean(K_total / phys_vol);
            const double R_mag_mean = amr_reductions_diagnostic.sum(c_R_mag_w) / phys_vol;
            // AMRReductions for evolution variables
            AMRReductions<VariableType::evolution> amr_reductions_evolution(
                m_cosmo_amr);

            double chi_mean = amr_reductions_evolution.sum(c_chi) / phys_vol;
            double chi_background = chi_mean;
            m_cosmo_amr.set_chi_mean(chi_background);

            CosmoDiagnostics<CubicHorndeskiWithCouplingAndPotential>
                cosmo_diag_pass2(cubic_horndeski, m_dx, m_p.G_Newton,
                                 chi_background);
            BoxLoops::loop(cosmo_diag_pass2, m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);
            fillAllGhosts();
            AMRReductions<VariableType::diagnostic> amr_diag2(m_cosmo_amr);
            const double zeta_mean  = amr_diag2.sum(c_zeta_w)  / phys_vol;
            const double zeta2_mean = amr_diag2.sum(c_zeta2_w) / phys_vol;
            const double zeta3_mean = amr_diag2.sum(c_zeta3_w) / phys_vol;
            
            const double variance    = zeta2_mean - zeta_mean * zeta_mean;
            const double third_central = zeta3_mean
                                       - 3.0 * zeta_mean * zeta2_mean
                                       + 2.0 * zeta_mean * zeta_mean * zeta_mean;
            const double skewness = (variance > 1e-30)
                                    ? third_central / std::pow(variance, 1.5)
                                    : 0.0;
            const double sigma_zeta = std::sqrt(std::max(variance, 0.0));

            // Write output file
            SmallDataIO data_out_file(m_p.data_path + "data_out", m_dt, m_time,
                                      m_restart_time, SmallDataIO::APPEND,
                                      first_step);
            data_out_file.remove_duplicate_time_data();
            if (first_step)
            {
                data_out_file.write_header_line(
                    {"L^2_Ham", "L^2_Mom", "<chi>", "<rho>", "<K>", "<zeta>", "sigma_zeta","zeta3_mean", "skewness"});
            }
            data_out_file.write_time_data_line({L2_Ham, L2_Mom, chi_mean,
                                                m_cosmo_amr.get_rho_mean(),
                                                m_cosmo_amr.get_K_mean(), zeta_mean, 
                                                sigma_zeta, zeta3_mean, skewness, R_mag_mean});
            
            // Use AMR Interpolator and do lineout data extraction
            // set up an interpolator
            // pass the boundary params so that we can use symmetries if
            // applicable
            AMRInterpolator<Lagrange<4>> interpolator(
                m_cosmo_amr, m_p.origin, m_p.dx, m_p.boundary_params,
                m_p.verbosity);

            // this should fill all ghosts including the boundary ones according
            // to the conditions set in params.txt
            interpolator.refresh();

            // set up the query and execute it
            std::array<double, CH_SPACEDIM> extraction_origin = {
                0., m_p.L / 2, m_p.L / 2}; // specified point {x \in [0,L],y \in
                                           // [0,L], z \in [0,L]}
            // rho lineout
            CustomExtraction rho_extraction(c_rho, m_p.lineout_num_points,
                                            m_p.L, extraction_origin, m_dt,
                                            m_time);
            rho_extraction.execute_query(&interpolator,
                                         m_p.data_path + "rho_lineout");
        }
    }
      if (m_p.activate_extraction == 1 || m_p.activate_scalar_extraction == 1)
    {
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
           CouplingAndPotential coupling_and_potential(
                m_p.coupling_and_potential_params);
           CubicHorndeskiWithCouplingAndPotential cubic_horndeski(
                coupling_and_potential);
            CosmoModifiedPunctureGauge cosmo_modified_puncture_gauge(
                m_p.modified_ccz4_params);
            cosmo_modified_puncture_gauge.set_K_mean(m_cosmo_amr.get_K_mean());
            ModifiedGravityWeyl4< CubicHorndeskiWithCouplingAndPotential,
                                 CosmoModifiedPunctureGauge, FourthOrderDerivatives>
                weyl4(cubic_horndeski, m_p.modified_ccz4_params,
                      cosmo_modified_puncture_gauge,
                      m_p.extraction_params.extraction_center, m_dx, m_p.sigma,
                      CCZ4RHS<>::USE_CCZ4);
            // CCZ4 is required since this code only works in this
            // formulation
            BoxLoops::loop(weyl4, m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);
            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                CH_TIME("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;
                m_cosmo_amr.m_interpolator->refresh(fill_ghosts);
                m_cosmo_amr.fill_multilevel_ghosts(
                    VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
                    min_level);
                if (m_p.activate_extraction)
                {
                    WeylExtraction weyl_extraction(m_p.extraction_params, m_dt,
                                                   m_time, first_step,
                                                   m_restart_time);
                    weyl_extraction.execute_query(m_cosmo_amr.m_interpolator);
                }
                if (m_p.activate_scalar_extraction)
                {
                    ScalarExtraction phi_extraction(
                        m_p.scalar_extraction_params, m_dt, m_time, first_step,
                        m_restart_time);
                    phi_extraction.execute_query(m_cosmo_amr.m_interpolator);
                }
            }
        }
    }
}
   
void CosmoLevel::postRestart()
{
    // only want to do this on the first restart and also every restart
    if (m_time == 0.0)
    {
        fillAllGhosts();
        CouplingAndPotential coupling_and_potential(
            m_p.coupling_and_potential_params);
        CubicHorndeskiWithCouplingAndPotential cubic_horndeski(
            coupling_and_potential)

        // Calculate constraints and some diagnostics as we need it in tagging
        // criterion
        BoxLoops::loop(ModifiedGravityConstraints<CubicHorndeskiWithCouplingAndPotential>(
                cubic_horndeski, m_dx, m_p.center, m_p.G_Newton, c_Ham,
                Interval(c_Mom1, c_Mom3), c_Ham_abs_sum,
                Interval(c_Mom_abs_sum1, c_Mom_abs_sum3)),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        CosmoDiagnostics<CubicHorndeskiWithCouplingAndPotential> cosmo_diagnostics(
            cubic_horndeski, m_dx, m_p.G_Newton);
        BoxLoops::loop(cosmo_diagnostics, m_state_new, m_state_diagnostics,
                       EXCLUDE_GHOST_CELLS);

        pout() << "Setting K mean on restart at t = " << m_time << " on level "
               << m_level << endl;

        // AMRReductions for diagnostic variables
        AMRReductions<VariableType::diagnostic> amr_reductions_diagnostic(
            m_cosmo_amr);
        double phys_vol = amr_reductions_diagnostic.sum(c_sqrt_gamma);
        double K_total = amr_reductions_diagnostic.sum(c_K_scaled);

        // Set rho_mean
        m_cosmo_amr.set_rho_mean(amr_reductions_diagnostic.sum(c_rho_scaled) /
                                 phys_vol);
        pout() << "Set rho_mean = " << m_cosmo_amr.get_rho_mean()
               << " at t = " << m_time << " on initial data at level "
               << m_level << endl;
        m_cosmo_amr.set_S_mean(amr_reductions_diagnostic.sum(c_S_scaled) /
                               phys_vol);
        // Set K_mean
        m_cosmo_amr.set_K_mean(K_total / phys_vol);
        pout() << "Calculated K mean as " << m_cosmo_amr.get_K_mean()
               << " at t = " << m_time << " on restart at level " << m_level
               << endl;

        // Use AMR Interpolator and do lineout data extraction
        // pass the boundary params so that we can use symmetries
        AMRInterpolator<Lagrange<2>> interpolator(m_cosmo_amr, m_p.origin,
                                                  m_p.dx, m_p.boundary_params,
                                                  m_p.verbosity);

        // this should fill all ghosts including the boundary ones according
        // to the conditions set in params.txt
        interpolator.refresh();

        // restart works from level 0 to highest level, so want this to happen
        // last on finest level
        int write_out_level = m_p.max_level;
        if (m_level == write_out_level)
        {
            // AMRReductions for diagnostic variables
            AMRReductions<VariableType::diagnostic> amr_reductions_diagnostic(
                m_gr_amr);
            double L2_Ham = amr_reductions_diagnostic.norm(c_Ham);
            double L2_Mom = amr_reductions_diagnostic.norm(Interval(c_Mom1, c_Mom3));

            // only on rank zero write out the result
            if (procID() == 0)
            {
                pout()
                    << "The initial norm of the constraint vars on restart is "
                    << L2_Ham << " for the Hamiltonian constraint and "
                    << L2_Mom << " for the momentum constraints" << endl;
            }

            // set up the query and execute it
            int num_points = 3 * m_p.ivN[0];
            ConstraintsExtraction constraints_extraction(
                c_Ham, Interval(c_Mom1, c_Mom3), num_points, m_p.L, m_p.center, m_dt, m_time);
            constraints_extraction.execute_query(
                &interpolator, m_p.data_path + "constraints_lineout");
        }
    }
}

#ifdef CH_USE_HDF5
// Things to do before a plot level - need to calculate the Wyl scalars
void CosmoLevel::prePlotLevel()
{
    fillAllGhosts();
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    CubicHorndeskiWithCouplingAndPotential cubic_horndeski(
        coupling_and_potential);
    BoxLoops::loop(
        ModifiedGravityConstraints<CubicHorndeskiWithCouplingAndPotential>(
            cubic_horndeski, m_dx,m_p.center, m_p.G_Newton, c_Ham, Interval(c_Mom1, c_Mom3)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    CosmoDiagnostics<CubicHorndeskiWithCouplingAndPotential> cosmo_diagnostics(
        cubic_horndeski, m_dx, m_p.G_Newton);
    BoxLoops::loop(cosmo_diagnostics, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
}
#endif /* CH_USE_HDF5 */

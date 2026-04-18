/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "ModifiedGravitySimulationParametersBase.hpp"

// Problem specific includes:
#include "InitialScalarData.hpp"
#include "CouplingAndPotential.hpp"
#include "CubicHorndeski.hpp"
#include "CosmoModifiedPunctureGauge.hpp"
class SimulationParameters : public ModifiedGravitySimulationParametersBase<
                                 CubicHorndeski<CouplingAndPotential>>
{
  public:
    SimulationParameters(GRParmParse &pp)
        : ModifiedGravitySimulationParametersBase(pp)
    {
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // Initial scalar field data
        initial_params.center =
            center; // already read in SimulationParametersBase
        pp.load("G_Newton", G_Newton, 1.0);
        pp.load("scalar_amplitude", initial_params.amplitude, 0.1);
        pp.load("scalar_mass", coupling_and_potential_params.scalar_mass);
        pp.load("g3", coupling_and_potential_params.g3);
        pp.load("g2", coupling_and_potential_params.g2);
        // Lineout params
        pp.load("lineout_num_points", lineout_num_points, 10);

        // Tagging params
        pp.load("tagging_center", tagging_center, center);
        pp.load("tagging_radius", tagging_radius, L);
    }

    void check_params()
    {
        warn_parameter("scalar_mass", coupling_and_potential_params.scalar_mass,
                       coupling_and_potential_params.scalar_mass <
                           0.2 / coarsest_dx / dt_multiplier,
                       "oscillations of scalar field do not appear to be "
                       "resolved on coarsest level");
        
       
    }

    // Initial data for matter and potential and BH
    double G_Newton, tagging_radius;
    int lineout_num_points;
    std::array<double, CH_SPACEDIM> tagging_center;
    InitialScalarData::params_t initial_params;
    CouplingAndPotential::params_t coupling_and_potential_params;


};

#endif /* SIMULATIONPARAMETERS_HPP_ */

/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COSMODIAGNOSTICS_HPP_
#define COSMODIAGNOSTICS_HPP_

#include "BSSNVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DiagnosticVariables.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "simd.hpp"
#include <array>

//! Calculates all relevant variables, which
//! are stored as diagnostics

template <class matter_t>
class CosmoDiagnostics // public MatterConstraints<matter_t>
{
  public:
    // Inherit the variable definitions from CCZ4 + matter_t
    template <class data_t>
    using TheoyVars = typename theory_t::template Vars<data_t>;

    template <class data_t>
    using TheoryVars = typename theory_t::template Vars<data_t>;
 
    template <class data_t>
    struct BSSNTheoryVars : public BSSNVars::VarsNoGauge<data_t>,
                            public TheoryVars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            CCZ4Vars<data_t>::enum_mapping(mapping_function);
            MatterVars<data_t>::enum_mapping(mapping_function);
        }
    };

    //! Constructor of class CosmoDiagnostics
    CosmoDiagnostics(const theory_t a_theory, double a_dx,
                     const std::array<double, CH_SPACEDIM> a_center);
    {
    }

    //! The compute member which calculates the diagnostics at each point in the
    //! box
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    theory_t m_ytheory; //!< The matter object, e.g. a scalar field
    const std::array<double, CH_SPACEDIM> m_center;
    const FourthOrderDerivatives m_deriv;
};

#include "CosmoDiagnostics.impl.hpp"

#endif /* COSMODIAGNOSTICS_HPP_ */

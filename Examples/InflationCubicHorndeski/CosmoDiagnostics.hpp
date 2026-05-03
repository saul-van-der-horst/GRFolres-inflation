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

template <class theory_t>
class CosmoDiagnostics
{
  public:
    // Combined BSSN + theory variables loaded at each cell
    template <class data_t>
    using TheoryVars = typename theory_t::template Vars<data_t>;

    template <class data_t>
    struct Vars : public BSSNVars::VarsNoGauge<data_t>,
                  public TheoryVars<data_t>
    {
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            BSSNVars::VarsNoGauge<data_t>::enum_mapping(mapping_function);
            TheoryVars<data_t>::enum_mapping(mapping_function);
        }
    };

    //! Constructor
    CosmoDiagnostics(const theory_t a_theory,
                     double a_dx,
                     double a_G_Newton,
                     double a_chi_background,
                     const std::array<double, CH_SPACEDIM> a_center)
        : m_theory(a_theory),
          m_deriv(a_dx),
          m_G_Newton(a_G_Newton),
          m_chi_background(a_chi_background),
          m_center(a_center)
    {
    }

    //! Compute diagnostics at each cell
    template <class data_t>
    void compute(Cell<data_t> current_cell) const;

  protected:
    theory_t m_theory;
    const FourthOrderDerivatives m_deriv;
    const double m_G_Newton;
    const double m_chi_background;   // a(t)^{-2}, set from spatial mean each step
    const std::array<double, CH_SPACEDIM> m_center;
};

#include "CosmoDiagnostics.impl.hpp"
#endif /* COSMODIAGNOSTICS_HPP_ */

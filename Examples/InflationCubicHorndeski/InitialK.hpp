/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALK_HPP
#define INITIALK_HPP

#include "Cell.hpp"
#include "Coordinates.hpp"            // needed for coords
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.H"
#include "ModifiedCCZ4RHS.hpp"
#include "simd.hpp"



template <class theory_t> class InitialK
{
  public:
    
    template <class data_t>
    using Vars = typename ModifiedCCZ4RHS<theory_t>::template Vars<data_t>;

    /
    template <class data_t>
    using Diff2Vars =
        typename ModifiedCCZ4RHS<theory_t>::template Diff2Vars<data_t>;

    //! Constructor
    InitialK(const theory_t &a_theory, const double dx,
             const std::array<double, CH_SPACEDIM> &a_center,
             double G_Newton = 1.0);

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const theory_t &m_theory;          // consistent naming throughout
    FourthOrderDerivatives m_deriv;
    const std::array<double, CH_SPACEDIM> m_center; // needed for Coordinates
    double m_G_Newton;
};

template <class theory_t>
InitialK<theory_t>::InitialK(const theory_t &a_theory, const double dx,
                              const std::array<double, CH_SPACEDIM> &a_center,
                              double G_Newton)
    : m_theory(a_theory), m_deriv(dx), m_center(a_center),
      m_G_Newton(G_Newton)
{
}

template <class theory_t>
template <class data_t>
void InitialK<theory_t>::compute(Cell<data_t> current_cell) const
{
    
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1   = m_deriv.template diff1<Vars>(current_cell);

    
    const auto d2   = m_deriv.template diff2<Diff2Vars>(current_cell);

    
    Coordinates<data_t> coords(current_cell, m_deriv.m_dx, m_center);

    
    const auto rho_and_Si =
        m_theory.compute_rho_and_Si(vars, d1, d2, coords);

    
    data_t K = -sqrt(24.0 * M_PI * m_G_Newton * rho_and_Si.rho);
    current_cell.store_vars(K, c_K);
}

#endif /* INITIALK_HPP */

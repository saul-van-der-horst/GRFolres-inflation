#ifndef COUPLINGANDPOTENTIAL_HPP_
#define COUPLINGANDPOTENTIAL_HPP_

#include <cmath>
#include <variant>


struct MatterParams
{
    
    double scalar_mass = 0.0;
    double ga3_scalar_mass = 0.0;
    double g2 = 0.0;
    double g3 = 0.0;
    double rbs_g3 = 0.0;
    double usr_y1 = 0.0;
    double y1 = 0.0;
    double y2 = 0.0;
    double v0 = 0.0;
    double usr_v0 = 0.0;
    double dbin_Lambda = 0.0;
    double rbs_Lambda = 0.0;
    double eta = 0.0;
    double f = 0.0;
    double mu = 0.0;
    double p = 0.0;
    double Mpl = 0.0;
    double nu = 0.0;
    double dbin_lambda1 = 0.0;
    double exph_lambda = 0.0;
    double v = 0.0;
    double b = 0.0;
    double gamma = 0.0;
    double eps = 0.0;
    double X_star = 0.0;
};

template <class Derived>
class CouplingAndPotential
{
  public:
    using params_t = MatterParams; 
    MatterParams m_params;

    
    explicit CouplingAndPotential(const MatterParams& a_params) : m_params(a_params) {}

    template <class data_t>
    ALWAYS_INLINE data_t V(const data_t phi, const data_t X) const
    { return static_cast<const Derived*>(this)->V_impl(phi, X); }

    template <class data_t>
    ALWAYS_INLINE data_t dV_dphi(const data_t phi, const data_t X) const
    { return static_cast<const Derived*>(this)->dV_dphi_impl(phi, X); }

    template <class data_t>
    ALWAYS_INLINE data_t G2(const data_t phi, const data_t X) const
    { return static_cast<const Derived*>(this)->G2_impl(phi, X); }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi(const data_t phi, const data_t X) const
    { return static_cast<const Derived*>(this)->dG2_dphi_impl(phi, X); }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dX(const data_t phi, const data_t X) const
    { return static_cast<const Derived*>(this)->dG2_dX_impl(phi, X); }

    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXX(const data_t phi, const data_t X) const
    { return static_cast<const Derived*>(this)->d2G2_dXX_impl(phi, X); }

    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXphi(const data_t phi, const data_t X) const
    { return static_cast<const Derived*>(this)->d2G2_dXphi_impl(phi, X); }

    template <class data_t>
    ALWAYS_INLINE data_t G3(const data_t phi, const data_t X) const
    { return static_cast<const Derived*>(this)->G3_impl(phi, X); }

    template <class data_t>
    ALWAYS_INLINE data_t dG3_dphi(const data_t phi, const data_t X) const
    { return static_cast<const Derived*>(this)->dG3_dphi_impl(phi, X); }

    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX(const data_t phi, const data_t X) const
    { return static_cast<const Derived*>(this)->dG3_dX_impl(phi, X); }

    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXX(const data_t phi, const data_t X) const
    { return static_cast<const Derived*>(this)->d2G3_dXX_impl(phi, X); }

    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXphi(const data_t phi, const data_t X) const
    { return static_cast<const Derived*>(this)->d2G3_dXphi_impl(phi, X); }

    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dphiphi(const data_t phi, const data_t X) const
    { return static_cast<const Derived*>(this)->d2G3_dphiphi_impl(phi, X); }
  protected:
    
    template <class data_t>
    ALWAYS_INLINE data_t V_impl(const data_t phi, const data_t X) const
    { return static_cast<data_t>(0); }

    template <class data_t>
    ALWAYS_INLINE data_t dV_dphi_impl(const data_t phi, const data_t X) const
    { return static_cast<data_t>(0); }

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    { return static_cast<data_t>(0); }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi_impl(const data_t phi, const data_t X) const
    { return static_cast<data_t>(0); }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dX_impl(const data_t phi, const data_t X) const
    { return static_cast<data_t>(0); }

    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXX_impl(const data_t phi, const data_t X) const
    { return static_cast<data_t>(0); }

    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXphi_impl(const data_t phi, const data_t X) const
    { return static_cast<data_t>(0); }

    template <class data_t>
    ALWAYS_INLINE data_t G3_impl(const data_t phi, const data_t X) const
    { return static_cast<data_t>(0); }

    template <class data_t>
    ALWAYS_INLINE data_t dG3_dphi_impl(const data_t phi, const data_t X) const
    { return static_cast<data_t>(0); }

    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX_impl(const data_t phi, const data_t X) const
    { return static_cast<data_t>(0); }

    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXX_impl(const data_t phi, const data_t X) const
    { return static_cast<data_t>(0); }

    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXphi_impl(const data_t phi, const data_t X) const
    { return static_cast<data_t>(0); }

    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dphiphi_impl(const data_t phi, const data_t X) const
    { return static_cast<data_t>(0); }
};
template <class Derived>
class Canonical : public CouplingAndPotential<Derived>
{
  public:
    explicit Canonical(const MatterParams& p) : CouplingAndPotential<Derived>(p) {}
};

template <class Derived>
class KEssence : public CouplingAndPotential<Derived>
{
  public:
    explicit KEssence(const MatterParams& p) : CouplingAndPotential<Derived>(p) {}
};

template <class Derived>
class KGB : public CouplingAndPotential<Derived>
{
  public:
    explicit KGB(const MatterParams& p) : CouplingAndPotential<Derived>(p) {}
};

class KGBCubic_galileon : public KGB<KGBCubic_galileon>
{
  public:
    explicit KGBCubic_galileon(const MatterParams& p) : KGB<KGBCubic_galileon>(p) {}

    
    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
        const auto X2= X * X;
        const double ga3_scalar_mass3= this->m_params.ga3_scalar_mass * this->m_params.ga3_scalar_mass * this->m_params.ga3_scalar_mass;
        return -X + X2/(2. * ga3_scalar_mass3 * this->m_params.mu);
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dX_impl(const data_t phi, const data_t X) const
    { 
      const double ga3_scalar_mass3= this->m_params.ga3_scalar_mass * this->m_params.ga3_scalar_mass * this->m_params.ga3_scalar_mass;
      return -1.+ X/(ga3_scalar_mass3 * this->m_params.mu); 
    }

    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXX_impl(const data_t phi, const data_t X) const
    { 
      const double ga3_scalar_mass3= this->m_params.ga3_scalar_mass * this->m_params.ga3_scalar_mass * this->m_params.ga3_scalar_mass;
      return 1./(ga3_scalar_mass3 * this->m_params.mu);
    }

    template <class data_t>
    ALWAYS_INLINE data_t G3_impl(const data_t phi, const data_t X) const
    { 
      const double ga3_scalar_mass3= this->m_params.ga3_scalar_mass * this->m_params.ga3_scalar_mass * this->m_params.ga3_scalar_mass;
      return X/(ga3_scalar_mass3) ;
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX_impl(const data_t phi, const data_t X) const
    { 
      const double ga3_scalar_mass3= this->m_params.ga3_scalar_mass * this->m_params.ga3_scalar_mass * this->m_params.ga3_scalar_mass;
      return 1./(ga3_scalar_mass3); 
    }
};

class KGBUltra_slow_roll : public KGB<KGBUltra_slow_roll>
{
  public:
    explicit KGBUltra_slow_roll(const MatterParams& p) : KGB<KGBUltra_slow_roll>(p) {}

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
      return - this->m_params.usr_v0 + this->m_params.y2 * phi;
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi_impl(const data_t phi, const data_t X) const
    { 
      return this->m_params.y2; 
    }

    

    template <class data_t>
    ALWAYS_INLINE data_t G3_impl(const data_t phi, const data_t X) const
    { return this->m_params.usr_y1 * X; }

    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX_impl(const data_t phi, const data_t X) const
    { return this->m_params.usr_y1; }
};

class KGBRunning_braiding_starobinsky : public KGB<KGBRunning_braiding_starobinsky>
{
  public:
    explicit KGBRunning_braiding_starobinsky(const MatterParams& p) : KGB<KGBRunning_braiding_starobinsky>(p) {}

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
      const double rbs_Lambda4= this->m_params.rbs_Lambda * this->m_params.rbs_Lambda * this->m_params.rbs_Lambda * this->m_params.rbs_Lambda;
      return -rbs_Lambda4 * (1.-exp(- this->m_params.nu * phi /this->m_params.Mpl)) * (1.-exp(- this->m_params.nu * phi /this->m_params.Mpl)) ;
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi_impl(const data_t phi, const data_t X) const
    { 
      const double rbs_Lambda4= this->m_params.rbs_Lambda * this->m_params.rbs_Lambda * this->m_params.rbs_Lambda * this->m_params.rbs_Lambda;
      return -rbs_Lambda4 * 2. * (this->m_params.nu / this->m_params.Mpl) * exp(- this->m_params.nu * phi /this->m_params.Mpl) *(1. - exp(- this->m_params.nu * phi /this->m_params.Mpl));
    }

    template <class data_t>
    ALWAYS_INLINE data_t G3_impl(const data_t phi, const data_t X) const
    { 
      return this->m_params.rbs_g3 * X * X; 
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX_impl(const data_t phi, const data_t X) const
    { 
      return 2. * this->m_params.rbs_g3 * X; 
    }

    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXX_impl(const data_t phi, const data_t X) const
    { 
      return 2. * this->m_params.rbs_g3; 
    }
};

class KGBExponential_hilltop : public KGB<KGBExponential_hilltop>
{
  public:
    explicit KGBExponential_hilltop(const MatterParams& p) : KGB<KGBExponential_hilltop>(p) {}

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
      const auto phi2= phi * phi;
      const double v2= this->m_params.v * this->m_params.v;
      return (phi2 - v2)*(phi2 - v2) * this->m_params.exph_lambda / 4.;
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi_impl(const data_t phi, const data_t X) const
    { 
      const auto phi2= phi * phi;
      const double v2= this->m_params.v * this->m_params.v;
      return this->m_params.exph_lambda * phi * (phi2 - v2); 
    }

    template <class data_t>
    ALWAYS_INLINE data_t G3_impl(const data_t phi, const data_t X) const
    { 
      return this->m_params.y1 * exp(this->m_params.eta * X /(1. + sqrt( X * X + this->m_params.eps * this->m_params.eps) * abs(this->m_params.eta))); 
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX_impl(const data_t phi, const data_t X) const
    { 
      const double eps2= this->m_params.eps * this->m_params.eps; 
      const auto sqrt_term= sqrt(X * X + eps2);
      const double abs_eta= abs(this->m_params.eta);
      return this->m_params.y1 * exp(this->m_params.eta * X /(1. + sqrt_term * abs_eta))*(eps2 * abs_eta+ sqrt_term)/(sqrt_term*(abs_eta*sqrt_term+1.)*(abs_eta*sqrt_term+1.)); 
    }

    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXX_impl(const data_t phi, const data_t X) const
    { 
      const auto X2= X * X;
      const auto X3= X2 * X;
      const double eps2= this->m_params.eps * this->m_params.eps; 
      const auto sqrt_term= sqrt(X2+ eps2);
      const double eta2= this->m_params.eta * this->m_params.eta;
      const double eta3= eta2 * this->m_params.eta;
      const double abs_eta= abs(this->m_params.eta);
      const double abs_eta3= abs_eta * abs_eta * abs_eta;
      const double eps4= eps2 * eps2;
      return (this->m_params.y1 * this->m_params.eta * exp(this->m_params.eta * X /(1. + sqrt_term * abs_eta))*(-3. * X *eps2* abs_eta3 * (X2 + eps2)+abs_eta*(2. * this->m_params.eta * X2 *eps2+ 2. * this->m_params.eta * eps4 -2. *X3-3. * X * eps2)+ eta2 * sqrt_term *(this->m_params.eta * eps4- 2. * X3-6. *X* eps2 )+ this->m_params.eta * pow(sqrt_term, 3.)))/(pow(sqrt_term, 3.) *(abs_eta* sqrt_term+1.)*(abs_eta* sqrt_term+1.)*(abs_eta* sqrt_term+1.)*(abs_eta* sqrt_term+1.)); 
    }
};

class KGBDefault : public KGB<KGBDefault>
{
  public:
    explicit KGBDefault(const MatterParams& p) : KGB<KGBDefault>(p) {}

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
      const auto phi2= phi * phi;
      const auto X2= X * X;
      const double scalar_mass2= this->m_params.scalar_mass * this->m_params.scalar_mass;
      return this->m_params.g2 * X2 - 0.5 * scalar_mass2 * phi2;
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi_impl(const data_t phi, const data_t X) const
    { 
      const double scalar_mass2= this->m_params.scalar_mass * this->m_params.scalar_mass;
      return -scalar_mass2 * phi; 
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dX_impl(const data_t phi, const data_t X) const
    { 
      return 2. * this->m_params.g2 * X;
    }

    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXX_impl(const data_t phi, const data_t X) const
    { 
      return 2. * this->m_params.g2; 
    }

    template <class data_t>
    ALWAYS_INLINE data_t G3_impl(const data_t phi, const data_t X) const
    { 
      return this->m_params.g3 * X; 
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX_impl(const data_t phi, const data_t X) const
    { 
      return this->m_params.g3; 
    }
};

class KGBDBI_natural : public KGB<KGBDBI_natural>
{
  public:
    explicit KGBDBI_natural(const MatterParams& p) : KGB<KGBDBI_natural>(p) {}

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
      const auto phi4= phi * phi * phi * phi;
      const double Lambda4= this->m_params.dbin_Lambda * this->m_params.dbin_Lambda * this->m_params.dbin_Lambda * this->m_params.dbin_Lambda;
      const double eps2= this->m_params.eps * this->m_params.eps;
      return - sqrt(1. - 2. * X+eps2) * phi4/this->m_params.dbin_lambda1 + phi4/this->m_params.dbin_lambda1 + Lambda4 * (1. + cos(phi/this->m_params.f));
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi_impl(const data_t phi, const data_t X) const
    { 
      const auto phi3= phi * phi * phi;
      const double Lambda4= this->m_params.dbin_Lambda * this->m_params.dbin_Lambda * this->m_params.dbin_Lambda * this->m_params.dbin_Lambda;
      const double eps2= this->m_params.eps * this->m_params.eps;
      return -sqrt(1. - 2. * X+eps2) * 4. * phi3/this->m_params.dbin_lambda1 +4. * phi3/this->m_params.dbin_lambda1 - Lambda4 * sin(phi/this->m_params.f)/this->m_params.f; 
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dX_impl(const data_t phi, const data_t X) const
    { 
      const auto phi4= phi * phi * phi * phi;
      const double eps2= this->m_params.eps * this->m_params.eps;
      return phi4/(this->m_params.dbin_lambda1 *sqrt(1. - 2. * X + eps2) ); 
    }

    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXX_impl(const data_t phi, const data_t X) const
    { 
      const double eps2= this->m_params.eps * this->m_params.eps;
      const double Dbi_term= 1.-2.*X+eps2;
      return (phi * phi * phi * phi)/(this->m_params.dbin_lambda1 * (sqrt(Dbi_term) * (Dbi_term))); 
    }
};

class KGBDBI_power_law : public KGB<KGBDBI_power_law>
{
  public:
    explicit KGBDBI_power_law(const MatterParams& p) : KGB<KGBDBI_power_law>(p) {}

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
      const double eps2= this->m_params.eps * this->m_params.eps;
      const double b2= this->m_params.b * this->m_params.b;
      const double pl_term1= (3. * (this->m_params.gamma +1)) / (4. * (b2 -1.));
      const auto exp_term= exp(2. *this->m_params.b *phi /this->m_params.Mpl);
      return - sqrt(1. - 2. * X + eps2) * this->m_params.v0 * exp_term /((this->m_params.gamma -1.) * (pl_term1)) +this->m_params.v0 * exp_term /((this->m_params.gamma -1.) * (pl_term1))- exp_term;
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi_impl(const data_t phi, const data_t X) const
    { 
      const double eps2= this->m_params.eps * this->m_params.eps;
      const double b2= this->m_params.b * this->m_params.b;
      const double pl_term1= (3. * (this->m_params.gamma +1)) / (4. * (b2 -1.));
      const auto exp_term= exp(2. *this->m_params.b *phi /this->m_params.Mpl);
      return -sqrt(1. - 2. * X + eps2) * this->m_params.v0 * 2. *this->m_params.b * exp_term /(this->m_params.Mpl * (this->m_params.gamma -1.) * (pl_term1))  +this->m_params.v0 * exp_term /(this->m_params.Mpl * (this->m_params.gamma -1.) * (pl_term1))  - 2. * this->m_params.b * exp_term/this->m_params.Mpl ; 
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dX_impl(const data_t phi, const data_t X) const
    { 
      const double eps2= this->m_params.eps * this->m_params.eps;
      const double b2= this->m_params.b * this->m_params.b;
      const double pl_term1= (3. * (this->m_params.gamma +1)) / (4. * (b2 -1.));
      const auto exp_term= exp(2. *this->m_params.b *phi /this->m_params.Mpl);
      return 1. + this->m_params.v0 * exp_term /((this->m_params.gamma -1.) * (pl_term1)) / sqrt(1. - 2. * X + eps2) ; 
    }

    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXX_impl(const data_t phi, const data_t X) const
    { 
      const double eps2= this->m_params.eps * this->m_params.eps;
      const double b2= this->m_params.b * this->m_params.b;
      const double pl_term1= (3. * (this->m_params.gamma +1)) / (4. * (b2 -1.));
      const auto exp_term= exp(2. *this->m_params.b *phi /this->m_params.Mpl);
      return this->m_params.v0 * exp_term /((this->m_params.gamma -1.) * (pl_term1)) /( sqrt(1. - 2. * X+ eps2) * (1. - 2. * X+ eps2)); 
    }
};


using MatterModelVariant = std::variant<
    KGBDefault,
    KGBUltra_slow_roll,
    KGBCubic_galileon,
    KGBRunning_braiding_starobinsky,
    KGBExponential_hilltop,
    KGBDBI_natural,
    KGBDBI_power_law>;

#endif /* COUPLINGANDPOTENTIAL_HPP_ */

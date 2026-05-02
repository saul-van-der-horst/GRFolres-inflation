#ifndef COUPLINGANDPOTENTIAL_HPP_
#define COUPLINGANDPOTENTIAL_HPP_


template <class Derived>
class CouplingAndPotential
{
  public:
    struct params_t
    {
        mutable double scalar_mass = 0.0;
        mutable double ga3_scalar_mass = 0.0;
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
    };
    params_t m_params;

    CouplingAndPotential(params_t a_params) : m_params(a_params) {}


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
    using params_t = typename CouplingAndPotential<Derived>::params_t;
    Canonical(params_t p) : CouplingAndPotential<Derived>(p) {}
};


template <class Derived>
class KEssence : public CouplingAndPotential<Derived>
{
  public:
    using params_t = typename CouplingAndPotential<Derived>::params_t;
    KEssence(params_t p) : CouplingAndPotential<Derived>(p) {}
};


template <class Derived>
class KGB : public CouplingAndPotential<Derived>
{
  public:
    using params_t = typename CouplingAndPotential<Derived>::params_t;
    KGB(params_t p) : CouplingAndPotential<Derived>(p) {}
};


class KGBCubic_galileon : public KGB<KGBCubic_galileon>
{
  public:
    using params_t = typename KGB<KGBCubic_galileon>::params_t;
    KGBCubic_galileon(params_t p) : KGB<KGBCubic_galileon>(p) {}

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
        return -X + pow(X, 2)/(2 * pow(this->m_params.ga3_scalar_mass, 3) * this->m_params.mu);
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dX_impl(const data_t phi, const data_t X) const
    { return -1+ X/(pow(this->m_params.ga3_scalar_mass, 3) * this->m_params.mu); }

    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXX_impl(const data_t phi, const data_t X) const
    { return 1/(pow(this->m_params.ga3_scalar_mass, 3) * this->m_params.mu); }

    template <class data_t>
    ALWAYS_INLINE data_t G3_impl(const data_t phi, const data_t X) const
    { return X/(pow(this->m_params.ga3_scalar_mass, 3)) ; }

    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX_impl(const data_t phi, const data_t X) const
    { return 1/(pow(this->m_params.ga3_scalar_mass, 3)); }
};
class KGBUltra_slow_roll : public KGB<KGBUltra_slow_roll>
{
  public:
    using params_t = typename KGB<KGBUltra_slow_roll>::params_t;
    KGBUltra_slow_roll(params_t p) : KGB<KGBUltra_slow_roll>(p) {}

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
        return -this-> m_params.usr_v0 + this-> m_params.y2 * phi;
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi_impl(const data_t phi, const data_t X) const
    { return this->m_params.y2; }

   

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
    using params_t = typename KGB<KGBRunning_braiding_starobinsky>::params_t;
    KGBRunning_braiding_starobinsky(params_t p) : KGB<KGBRunning_braiding_starobinsky>(p) {}

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
        return -pow(this->m_params.rbs_Lambda, 4) * pow(1-exp(- this->m_params.nu * phi /this->m_params.Mpl), 2) ;
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi_impl(const data_t phi, const data_t X) const
    { return -pow(this->m_params.rbs_Lambda, 4) * 2 * (this->m_params.nu / this->m_params.Mpl) * exp(- this->m_params.nu * phi /this->m_params.Mpl) *
       (1-exp(- this->m_params.nu * phi /this->m_params.Mpl)); }

    template <class data_t>
    ALWAYS_INLINE data_t G3_impl(const data_t phi, const data_t X) const
    { return this->m_params.rbs_g3 * pow(X, 2); }

    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX_impl(const data_t phi, const data_t X) const
    { return 2. * this->m_params.rbs_g3 * X; }

    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXX_impl(const data_t phi, const data_t X) const
    { return 2. * this->m_params.rbs_g3; }
};

class KGBExponential_hilltop : public KGB<KGBExponential_hilltop>
{
  public:
    using params_t = typename KGB<KGBExponential_hilltop>::params_t;
    KGBExponential_hilltop(params_t p) : KGB<KGBExponential_hilltop>(p) {}

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
        return  pow((pow(phi, 2)-pow(this->m_params.v, 2)), 2) * this->m_params.exph_lambda / 4. ;
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi_impl(const data_t phi, const data_t X) const
    { return this->m_params.exph_lambda * phi * (pow(phi, 2)-pow(this->m_params.v, 2)); }

    

    

    template <class data_t>
    ALWAYS_INLINE data_t G3_impl(const data_t phi, const data_t X) const
    { return this->m_params.y1 * exp(this->m_params.eta * X); }

    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX_impl(const data_t phi, const data_t X) const
    { return this->m_params.y1 * this->m_params.eta* exp(this->m_params.eta * X); }

    template <class data_t>
    ALWAYS_INLINE data_t d2G3_dXX_impl(const data_t phi, const data_t X) const
    { return this->m_params.y1 * pow(this->m_params.eta, 2)* exp(this->m_params.eta * X); }
};
class KGBDefault : public KGB<KGBDefault>
{
  public:
    using params_t = typename KGB<KGBDefault>::params_t;
    KGBDefault(params_t p) : KGB<KGBDefault>(p) {}

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
        return this->m_params.g2 * X * X
             - 0.5 * this->m_params.scalar_mass * this->m_params.scalar_mass * phi * phi;
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi_impl(const data_t phi, const data_t X) const
    { return -this->m_params.scalar_mass * this->m_params.scalar_mass * phi; }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dX_impl(const data_t phi, const data_t X) const
    { return 2. * this->m_params.g2 * X; }

    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXX_impl(const data_t phi, const data_t X) const
    { return 2. * this->m_params.g2; }

    template <class data_t>
    ALWAYS_INLINE data_t G3_impl(const data_t phi, const data_t X) const
    { return this->m_params.g3 * X; }

    template <class data_t>
    ALWAYS_INLINE data_t dG3_dX_impl(const data_t phi, const data_t X) const
    { return this->m_params.g3; }
};
class KGBDBI_natural : public KGB<KGBDBI_natural>
{
  public:
    using params_t = typename KGB<KGBDBI_natural>::params_t;
    KGBDBI_natural(params_t p) : KGB<KGBDBI_natural>(p) {}

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
        return -sqrt(1. - 2. * X) *pow(phi, 4)/this->m_params.dbin_lambda1 +pow(phi, 4)/this->m_params.dbin_lambda1 + pow(this->m_params.dbin_Lambda, 4) * (1. +cos(phi/this->m_params.f));
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi_impl(const data_t phi, const data_t X) const
    { return -sqrt(1. - 2. * X) * 4. * pow(phi, 3)/this->m_params.dbin_lambda1 +4. * pow(phi, 3)/this->m_params.dbin_lambda1 -pow(this->m_params.dbin_Lambda, 4) * sin(phi/this->m_params.f); }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dX_impl(const data_t phi, const data_t X) const
    { return pow(phi, 4)/(this->m_params.dbin_lambda1 / sqrt(1. - 2. * X) ); }

    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXX_impl(const data_t phi, const data_t X) const
    { return pow(phi, 4)/(this->m_params.dbin_lambda1 / (sqrt(1. - 2. * X) * (1. - 2. * X))); }

};

class KGBDBI_power_law : public KGB<KGBDBI_power_law>
{
  public:
    using params_t = typename KGB<KGBDBI_power_law>::params_t;
    KGBDBI_power_law(params_t p) : KGB<KGBDBI_power_law>(p) {}

    template <class data_t>
    ALWAYS_INLINE data_t G2_impl(const data_t phi, const data_t X) const
    {
        return -sqrt(1. - 2. * X) * this->m_params.v0 * exp(2 *this->m_params.b *phi /this->m_params.Mpl) 
          /((this->m_params.gamma -1) * ((3 * (this->m_params.gamma +1)) / (4 * (pow(this->m_params.b, 2) -1)))) +this->m_params.v0 * exp(2 *this->m_params.b *phi /this->m_params.Mpl) 
          /((this->m_params.gamma -1) * ((3 * (this->m_params.gamma +1)) / (4 * (pow(this->m_params.b, 2) -1))))-exp(2 *this->m_params.b *phi /this->m_params.Mpl);
    }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dphi_impl(const data_t phi, const data_t X) const
    { return -sqrt(1. - 2. * X) * this->m_params.v0 * 2 *this->m_params.b * exp(2 *this->m_params.b *phi /this->m_params.Mpl) 
          /(this->m_params.Mpl * (this->m_params.gamma -1) * ((3 * (this->m_params.gamma +1)) / (4 * (pow(this->m_params.b, 2) -1))))  +this->m_params.v0 * exp(2 *this->m_params.b *phi /this->m_params.Mpl) 
          /(this->m_params.Mpl * (this->m_params.gamma -1) * ((3 * (this->m_params.gamma +1)) / (4 * (pow(this->m_params.b, 2) -1))))  - 2 * this->m_params.b *exp(2 *this->m_params.b * phi /this->m_params.Mpl)/this->m_params.Mpl ; }

    template <class data_t>
    ALWAYS_INLINE data_t dG2_dX_impl(const data_t phi, const data_t X) const
    { return this->m_params.v0 * exp(2 *this->m_params.b *phi /this->m_params.Mpl) 
          /((this->m_params.gamma -1) * ((3 * (this->m_params.gamma +1)) / (4 * (pow(this->m_params.b, 2) -1)))) / sqrt(1. - 2. * X) ; }

    template <class data_t>
    ALWAYS_INLINE data_t d2G2_dXX_impl(const data_t phi, const data_t X) const
    { return this->m_params.v0 * exp(2 *this->m_params.b *phi /this->m_params.Mpl) 
          /((this->m_params.gamma -1) * ((3 * (this->m_params.gamma +1)) / (4 * (pow(this->m_params.b, 2) -1)))) /( sqrt(1. - 2. * X) * (1. - 2. * X)); }

};

struct ICouplingAndPotential
{
    virtual ~ICouplingAndPotential() = default;

    virtual double V            (double phi, double X) const = 0;
    virtual double dV_dphi      (double phi, double X) const = 0;
    virtual double G2           (double phi, double X) const = 0;
    virtual double dG2_dphi     (double phi, double X) const = 0;
    virtual double dG2_dX       (double phi, double X) const = 0;
    virtual double d2G2_dXX     (double phi, double X) const = 0;
    virtual double d2G2_dXphi   (double phi, double X) const = 0;
    virtual double G3           (double phi, double X) const = 0;
    virtual double dG3_dphi     (double phi, double X) const = 0;
    virtual double dG3_dX       (double phi, double X) const = 0;
    virtual double d2G3_dXX     (double phi, double X) const = 0;
    virtual double d2G3_dXphi   (double phi, double X) const = 0;
    virtual double d2G3_dphiphi (double phi, double X) const = 0;
};

template <class Model>
struct ModelWrapper : public ICouplingAndPotential
{
    Model m_model;
    ModelWrapper(typename Model::params_t p) : m_model(p) {}

    double V            (double phi, double X) const override { return m_model.V(phi, X); }
    double dV_dphi      (double phi, double X) const override { return m_model.dV_dphi(phi, X); }
    double G2           (double phi, double X) const override { return m_model.G2(phi, X); }
    double dG2_dphi     (double phi, double X) const override { return m_model.dG2_dphi(phi, X); }
    double dG2_dX       (double phi, double X) const override { return m_model.dG2_dX(phi, X); }
    double d2G2_dXX     (double phi, double X) const override { return m_model.d2G2_dXX(phi, X); }
    double d2G2_dXphi   (double phi, double X) const override { return m_model.d2G2_dXphi(phi, X); }
    double G3           (double phi, double X) const override { return m_model.G3(phi, X); }
    double dG3_dphi     (double phi, double X) const override { return m_model.dG3_dphi(phi, X); }
    double dG3_dX       (double phi, double X) const override { return m_model.dG3_dX(phi, X); }
    double d2G3_dXX     (double phi, double X) const override { return m_model.d2G3_dXX(phi, X); }
    double d2G3_dXphi   (double phi, double X) const override { return m_model.d2G3_dXphi(phi, X); }
    double d2G3_dphiphi (double phi, double X) const override { return m_model.d2G3_dphiphi(phi, X); }
};

inline std::unique_ptr<ICouplingAndPotential>
makeCouplingAndPotential(const std::string &model_name,
                         CouplingAndPotential<KGBDefault>::params_t params)
{
    if (model_name == "kgb-default")
        return std::make_unique<ModelWrapper<KGBDefault>>(params);
    else if (model_name == "kgb-usr")
        return std::make_unique<ModelWrapper<KGBUltra_slow_roll>>(params);
    else if (model_name == "kgb-ga3")
        return std::make_unique<ModelWrapper<KGBCubic_galileon>>(params);
    else if (model_name == "kgb-rbs")
        return std::make_unique<ModelWrapper<KGBRunning_braiding_starobinsky>>(params);
    else if (model_name == "kgb-exph")
        return std::make_unique<ModelWrapper<KGBExponential_hilltop>>(params);
    else if (model_name == "kgb-dbin")
        return std::make_unique<ModelWrapper<KGBDBI_natural>>(params);
    else if (model_name == "kgb-dbipl")
        return std::make_unique<ModelWrapper<KGBDBI_power_law>>(params);
    else
        throw std::invalid_argument("Unknown model: '" + model_name + "'");
}

#endif /* COUPLINGANDPOTENTIAL_HPP_ */

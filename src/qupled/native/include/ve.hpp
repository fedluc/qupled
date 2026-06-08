#ifndef VE_HPP
#define VE_HPP

#include "input.hpp"
#include "qstls.hpp"
#include <memory>

class VE : public Qstls {

public:

  explicit VE(const std::shared_ptr<const VEInput> &in_);

  const std::vector<double> &getFreeEnergyGrid() const { return freeEnergyGridOut; }
  const std::vector<std::vector<double>> &getFreeEnergyIntegrand() const {
    return freeEnergyIntegrandOut;
  }
  const std::vector<double> &getVCoeff() const { return aCoeff; }
  const std::vector<double> &getA1Coeff() const { return a1Coeff; }
  const std::vector<double> &getMxcL() const { return mxcLDiag; }
  const std::vector<double> &getMatsubaraGrid() const { return matsubaraGrid; }
  double getA0Coeff() const { return a0Coeff; }
  double getKxc0() const { return kxc0; }
  double getMxcInf() const { return mxcInf; }
  double getCxc() const { return cxc; }
  double getOmegaM2() const { return omegaM2; }

protected:

  void init() override;
  void computeLfc() override;

  const VEInput &in() const {
    return *StlsUtil::dynamic_pointer_cast<Input, VEInput>(inPtr);
  }

private:

  Vector3D adrFixedVe;
  std::vector<double> aCoeff;
  std::vector<double> a1Coeff;
  std::vector<double> mxcLDiag;
  std::vector<double> matsubaraGrid;
  std::vector<double> freeEnergyGridOut;
  std::vector<std::vector<double>> freeEnergyIntegrandOut;
  double a0Coeff;
  double kxc0;
  double mxcInf;
  double cxc;
  double omegaM2;

  void computeAdrFixedVe();
  void initializeFreeEnergyIntegrand();
  void updateFreeEnergyIntegrandCurrent();
  double currentFreeEnergyIntegrand() const;
  double exactZeroCouplingIntegrand() const;
  double computeStaticKxc0() const;
  double computeSmallXCoefficient(const Vector2D &psi, size_t l) const;
  std::vector<double> computeA1(const Vector2D &psi1) const;
  double estimateCxcInitial(const std::vector<double> &a1) const;
  double estimateCxcFromMxc(const std::vector<double> &mxcL) const;
};

namespace VEUtil {

  class AdrFixed : public QstlsUtil::AdrFixedBase {

  public:

    AdrFixed(const double &Theta_,
             const double &yMin_,
             const double &yMax_,
             const double &x_,
             const double &mu_,
             const std::vector<double> &itgGrid_,
             std::shared_ptr<Integrator2D> itg_)
        : QstlsUtil::AdrFixedBase(Theta_, yMin_, yMax_, x_, mu_),
          itg(itg_),
          itgGrid(itgGrid_) {}

    void get(const std::vector<double> &wvg, Vector3D &res) const;

  private:

    double integrand1(const double &y, const double &l) const;
    double integrand2(const double &s, const double &w, const double &l) const;
    double q(const double &s) const;

    const std::shared_ptr<Integrator2D> itg;
    const std::vector<double> &itgGrid;
  };

} // namespace VEUtil

#endif

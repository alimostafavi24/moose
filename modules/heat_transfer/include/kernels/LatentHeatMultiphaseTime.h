
#pragma once

#include "TimeDerivative.h"
#include "Material.h"

class LatentHeatMultiphaseTime : public TimeDerivative
{
public:
  LatentHeatMultiphaseTime(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;


  /// Multiplier for the coupled force term

  

  const VariableValue & _c;  // w is the chemical potential
  const MaterialProperty<Real> & _dfldT;

  const Real & _rho_l;
  const Real & _latent_heat;



};


#pragma once

#include "Kernel.h"

class ColumbForce : public Kernel
{
public:
  ColumbForce(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
 // virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Multiplier for the coupled force term

  
  const Real _coef;
  const VariableGradient & _grad_phit; // c is the order parameter 
  const VariableValue & _rhoet;  // w is the chemical potential


  unsigned int _component;

};

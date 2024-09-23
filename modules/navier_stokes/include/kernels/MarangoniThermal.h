
#pragma once

#include "Kernel.h"

class MarangoniThermal : public Kernel
{
public:
  MarangoniThermal(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const Real _coef;
  const Real _gamma_T;



  /// Multiplier for the coupled force term

  

  const VariableGradient & _grad_c; // c is the order parameter 
  const VariableGradient & _grad_T;  // w is the chemical potential


  unsigned int _component;
  unsigned int _c_var;
  unsigned int _T_var;

};

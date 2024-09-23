
#pragma once

#include "Kernel.h"

class Abels : public Kernel
{
public:
  Abels(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const Real _rho2;
  const Real _rho1;
  const Real _mobility;


  /// Multiplier for the coupled force term

//  const VariableValue & _u;  // w is the chemical potential
//  const VariableGradient & _grad_u;  // w is the chemical potential
  const VariableGradient & _grad_w;  // w is the chemical potential


  //unsigned int _component;
 // unsigned int _vel_x_var;
  unsigned int _w_var;

};

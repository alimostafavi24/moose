
#pragma once

#include "Kernel.h"

class ElectricForce22 : public Kernel
{
public:
  ElectricForce22(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const Real _epsilon_zero;
  const Real _epsilon_water;
  const Real _epsilon_air;


  /// Multiplier for the coupled force term

  
  const VariableGradient & _grad_c; // c is the order parameter
  const VariableValue & _phie;  // w is the chemical potential 
  const VariableGradient & _grad_phie; // c is the order parameter 


  unsigned int _component;
  unsigned int _c_var;
  unsigned int _phie_var;

};

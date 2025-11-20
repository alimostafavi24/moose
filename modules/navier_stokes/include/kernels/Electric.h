
#pragma once

#include "Kernel.h"

class Electric : public Kernel
{
public:
  Electric(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

 // virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Multiplier for the coupled force term
  const MaterialProperty<Real> & _permittivity;
  const MaterialProperty<Real> & _dpermittivitydc;
  const MaterialProperty<Real> & _d2permittivitydc2;

  const VariableGradient & _grad_phit; // phit is the electric potential 
  const VariableValue & _rhoe;  // rhoe is the charge density
  const VariableGradient & _grad_c;
   
  unsigned int _component;
  unsigned int _phit_var;
  unsigned int _rhoe_var;
  unsigned int _c_var;

};

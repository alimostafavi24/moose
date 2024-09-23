
#pragma once

#include "Kernel.h"

class SurfaceTensionThermal : public Kernel
{
public:
  SurfaceTensionThermal(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  
  const MaterialProperty<Real> & _f;
  const Real _coef;



  /// Multiplier for the coupled force term

  

  const VariableGradient & _grad_c; // c is the order parameter 
  const VariableValue & _w;  // w is the chemical potential


  unsigned int _component;
  unsigned int _c_var;
  unsigned int _w_var;

};

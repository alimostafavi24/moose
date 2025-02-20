//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSBase.h"

// Forward Declarations

/**
 * This class computes the mass equation residual and Jacobian
 * contributions for the incompressible Navier-Stokes momentum
 * equation.
 */
class INSMass : public INSBase
{
public:
  static InputParameters validParams();

  INSMass(const InputParameters & parameters);

  virtual ~INSMass() {}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  virtual Real computeQpPGResidual();
  virtual Real computeQpPGJacobian();
  virtual Real computeQpPGOffDiagJacobian(unsigned comp);

  bool _pspg;
  bool _multiphase; // New parameter for multiphase consideration
  const Function & _x_ffn;
  const Function & _y_ffn;
  const Function & _z_ffn;
  const VariableGradient & _grad_c; // c is the order parameter 
  const VariableValue & _w;  // w is the chemical potential
};

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.KEViscousWater1RZ, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.KEViscousWater1RZ.html

#pragma once

#include "ElementIntegralPostprocessor.h"
#include "MooseVariableInterface.h"

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class KEViscousWater1RZ : public ElementIntegralPostprocessor,
                                             public MooseVariableInterface<Real>
{
public:
  static InputParameters validParams();

  KEViscousWater1RZ(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  /// Holds the solution at current quadrature points
  const VariableValue & _u;
  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;

  const VariableValue & _vel_x;
  const VariableGradient & _grad_vel_x;
  const VariableSecond & _second_vel_x;

    const VariableValue & _vel_y;
  const VariableGradient & _grad_vel_y;
  const VariableSecond & _second_vel_y;

  const VariableValue & _c;
  const VariableGradient & _grad_c;
  
  //const Real & _mu;
  const MaterialProperty<Real> & _mu;

};

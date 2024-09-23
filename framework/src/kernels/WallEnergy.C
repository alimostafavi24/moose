//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WallEnergy.h"
// MOOSE includes


registerMooseObject("MooseApp", WallEnergy);

InputParameters
WallEnergy::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredCoupledVar("variable", "The name of the variable that this object operates on");

  params.addRequiredCoupledVar("c", "c");
  params.addParam<Real>("gamma", 1.0, "gamma");
  params.addParam<Real>("cos", 1.0, "cos");
  params.addClassDescription("Computes a volume integral of the specified variable");
  return params;
}

WallEnergy::WallEnergy(
    const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this,
                                 false,
                                 "variable",
                                 Moose::VarKindType::VAR_ANY,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD),
    _u(coupledValue("variable")),
    _grad_u(coupledGradient("variable")),
    _c(coupledValue("c")),
    _gamma(getParam<Real>("gamma")),
    _cos(getParam<Real>("cos"))


{
  addMooseVariableDependency(&mooseVariableField());
}

Real
WallEnergy::computeQpIntegral()
{
  if (_q_point[_qp](1) == 0.0)
  return -1.0 * _gamma * _cos * 0.25 * _c[_qp] * (3.0-_c[_qp]*_c[_qp]);
  else
  return 0;


 }
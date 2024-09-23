//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TotalDiffusion.h"
// MOOSE includes


registerMooseObject("MooseApp", TotalDiffusion);

InputParameters
TotalDiffusion::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredCoupledVar("variable", "The name of the variable that this object operates on");
  params.addRequiredCoupledVar("c", "c");

  params.addClassDescription("Computes a volume integral of the specified variable");
  params.addParam<MaterialPropertyName>("mobility", "mobility", "The name of the density");
   params.addParam<Real>("coef", 1.0, "Coefficient");

  return params;
}

TotalDiffusion::TotalDiffusion(
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
    _mobility(getMaterialProperty<Real>("mobility")),
    _coef(getParam<Real>("coef"))

{
  addMooseVariableDependency(&mooseVariableField());
}

Real
TotalDiffusion::computeQpIntegral()
{
 // if (_c[_qp] < -0.0)
 // return 0;

  return 1 * _mobility[_qp] * _coef * _coef *
  ( _grad_u[_qp](0) *  _grad_u[_qp](0) +
    _grad_u[_qp](1) *  _grad_u[_qp](1));
 }
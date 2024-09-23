

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MarangoniThermal.h"
#include "NS.h"
#include "MooseVariableFE.h"
#include "SystemBase.h"

registerMooseObject("NavierStokesApp", MarangoniThermal);

InputParameters
MarangoniThermal::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Adds the surface tension force to the NS equation.");  
 params.addParam<Real>("coef", 1.0, "Coefficient");
  params.addParam<Real>("gamma_T", 1.0, "Marangoni coef");
  params.addRequiredParam<unsigned int>("component", "number of component (0 = x, 1 = y, 2 = z)");
  params.addRequiredCoupledVar("c", "The order parameter");
  params.addRequiredCoupledVar("T", "Temperature");


  return params;
}

MarangoniThermal::MarangoniThermal(const InputParameters & parameters)
  : Kernel(parameters),
    _coef(getParam<Real>("coef")),
    _gamma_T(getParam<Real>("gamma_T")),
    _grad_c(coupledGradient("c")),
    _grad_T(coupledGradient("T")),
    _component(getParam<unsigned int>("component")), 

    _c_var(coupled("c")),
    _T_var(coupled("T"))
{
}

Real
MarangoniThermal::computeQpResidual()
{
 //return -_c[_qp] * _w[_qp] * _grad_test[_i][_qp](_component);
  Real A = (_grad_c[_qp](0) * _grad_c[_qp](0)) + (_grad_c[_qp](1) * _grad_c[_qp](1));

  Real B = ( _grad_T[_qp](0) *  _grad_c[_qp](0)) + ( _grad_T[_qp](1) *  _grad_c[_qp](1));


 return _coef * ((_gamma_T * A * _grad_T[_qp](_component))-(_gamma_T * B * _grad_c[_qp](_component)));

//   Real A = _gamma_T * (_grad_c[_qp](0) * _grad_c[_qp](0)) + (_grad_c[_qp](1) * _grad_c[_qp](1));

  //Real B = _gamma_T * ( _grad_T[_qp](0) *  _grad_c[_qp](0)) + ( _grad_T[_qp](1) *  _grad_c[_qp](1));


 //return _coef * (-1 * A * _grad_T[_qp](_component) + B * _grad_c[_qp](_component));

}

Real
MarangoniThermal::computeQpJacobian()
{
  return 0.;
}

Real
MarangoniThermal::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _c_var)
  {
    return 0;
  }
  else if (jvar == _T_var)
  {
    return 0;
  }
    return 0.;
}
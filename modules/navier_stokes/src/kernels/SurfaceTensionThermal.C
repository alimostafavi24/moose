

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SurfaceTensionThermal.h"
#include "NS.h"
#include "MooseVariableFE.h"
#include "SystemBase.h"

registerMooseObject("NavierStokesApp", SurfaceTensionThermal);

InputParameters
SurfaceTensionThermal::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Adds the surface tension force to the NS equation."); 
  params.addParam<MaterialPropertyName>("f_name", "f", "f function"); 
  params.addParam<Real>("coef", 1.0, "Coefficient");
  params.addRequiredParam<unsigned int>("component", "number of component (0 = x, 1 = y, 2 = z)");
  params.addRequiredCoupledVar("c", "The order parameter");
  params.addRequiredCoupledVar("w", "The chemical potential");


  return params;
}

SurfaceTensionThermal::SurfaceTensionThermal(const InputParameters & parameters)
  : Kernel(parameters),
    _f(getMaterialProperty<Real>("f_name")),
    _coef(getParam<Real>("coef")),
    _grad_c(coupledGradient("c")),
    _w(coupledValue("w")),
    _component(getParam<unsigned int>("component")), 

    _c_var(coupled("c")),
    _w_var(coupled("w"))
{
}

Real
SurfaceTensionThermal::computeQpResidual()
{
 //return -_c[_qp] * _w[_qp] * _grad_test[_i][_qp](_component);
 return _coef * _w[_qp] *_f[_qp]*  _grad_c[_qp](_component)* _test[_i][_qp];

}

Real
SurfaceTensionThermal::computeQpJacobian()
{
  return 0.;
}

Real
SurfaceTensionThermal::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _c_var)
  {
    return _coef * _w[_qp] *_f[_qp]*  _grad_phi[_j][_qp](_component)* _test[_i][_qp];;
  }
  else if (jvar == _w_var)
  {
    return _coef * _phi[_j][_qp] *_f[_qp]*  _grad_c[_qp](_component)* _test[_i][_qp];;
  }
    return 0.;
}
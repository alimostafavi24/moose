

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SurfaceTension3Multi.h"
#include "NS.h"
#include "MooseVariableFE.h"
#include "SystemBase.h"

registerMooseObject("NavierStokesApp", SurfaceTension3Multi);

InputParameters
SurfaceTension3Multi::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Adds the surface tension force to the NS equation.");  
 params.addParam<Real>("coef", 1.0, "Coefficient");
  params.addRequiredParam<unsigned int>("component", "number of component (0 = x, 1 = y, 2 = z)");
  params.addRequiredCoupledVar("c", "The order parameter");
  params.addRequiredCoupledVar("w", "The chemical potential");


  return params;
}

SurfaceTension3Multi::SurfaceTension3Multi(const InputParameters & parameters)
  : Kernel(parameters),
    _coef(getParam<Real>("coef")),
    _grad_c(coupledGradient("c")),
    _w(coupledValue("w")),
    _component(getParam<unsigned int>("component")) 

 //   _c_var(coupled("c")),
//    _w_var(coupled("w"))
{
}

Real
SurfaceTension3Multi::computeQpResidual()
{
 //return -_c[_qp] * _w[_qp] * _grad_test[_i][_qp](_component);
 return _coef * _w[_qp] * _grad_c[_qp](_component)* _test[_i][_qp];

}

Real
SurfaceTension3Multi::computeQpJacobian()
{
  return 0.;
}

Real
SurfaceTension3Multi::computeQpOffDiagJacobian(unsigned jvar)
{

    return 0.;
}
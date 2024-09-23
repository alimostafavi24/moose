

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ColumbForce.h"
#include "NS.h"
#include "MooseVariableFE.h"
#include "SystemBase.h"

registerMooseObject("NavierStokesApp", ColumbForce);

InputParameters
ColumbForce::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addParam<Real>("coef", 1.0, "coef");
  params.addClassDescription("Adds the Columb force to the NS equation.");  
  params.addRequiredParam<unsigned int>("component", "number of component (0 = x, 1 = y, 2 = z)");
  params.addRequiredCoupledVar("phit", "Potential");
  params.addRequiredCoupledVar("rhoet", "Charge density");


  return params;
}

ColumbForce::ColumbForce(const InputParameters & parameters)
  : Kernel(parameters),
    _coef(getParam<Real>("coef")),
    _grad_phit(coupledGradient("phit")),
    _rhoet(coupledValue("rhoet")),
    _component(getParam<unsigned int>("component"))
{
}

Real
ColumbForce::computeQpResidual()
{
 //return -_c[_qp] * _w[_qp] * _grad_test[_i][_qp](_component);
 return  _coef * _rhoet[_qp] * _grad_phit[_qp](_component)* _test[_i][_qp];

}

Real
ColumbForce::computeQpJacobian()
{
  return 0.;
}

//Real
//ColumbForce::computeQpOffDiagJacobian(unsigned jvar)
//{
//    return 0.;
//}
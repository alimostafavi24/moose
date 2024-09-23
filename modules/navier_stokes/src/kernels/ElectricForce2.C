

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElectricForce2.h"
#include "NS.h"
#include "MooseVariableFE.h"
#include "SystemBase.h"

registerMooseObject("NavierStokesApp", ElectricForce2);

InputParameters
ElectricForce2::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Adds the Electric force to the NS equation.");  
  params.addParam<Real>("epsilon_zero", 1.0, "epsilon_zero");
  params.addParam<Real>("epsilon_water", 1.0, "epsilon_water");
  params.addParam<Real>("epsilon_air", 1.0, "epsilon_air");
  params.addRequiredParam<unsigned int>("component", "number of component (0 = x, 1 = y, 2 = z)");
  params.addRequiredCoupledVar("phie", "phie");
  params.addRequiredCoupledVar("c", "c");

  return params;
}

ElectricForce2::ElectricForce2(const InputParameters & parameters)
  : Kernel(parameters),
    _epsilon_zero(getParam<Real>("epsilon_zero")),
    _epsilon_water(getParam<Real>("epsilon_water")),
    _epsilon_air(getParam<Real>("epsilon_air")),
    _grad_c(coupledGradient("c")), 
    _phie(coupledValue("phie")),
    _grad_phie(coupledGradient("phie")),
    _component(getParam<unsigned int>("component")), 
    _c_var(coupled("c"))
 //   _phie_var(coupled("phie"))

{
}

Real
ElectricForce2::computeQpResidual()
{
 //return -_c[_qp] * _w[_qp] * _grad_test[_i][_qp](_component);
 return 0.5 * _epsilon_zero * (_grad_phie[_qp](0)*_grad_phie[_qp](0)+_grad_phie[_qp](1)*_grad_phie[_qp](1)) *
 0.5 * (_epsilon_water - _epsilon_air) * (_grad_c[_qp](_component)) * _test[_i][_qp];

}

Real
ElectricForce2::computeQpJacobian()
{
  return 0.;
}

Real
ElectricForce2::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _c_var)
  {
 return 0.5 * _epsilon_zero * (_grad_phie[_qp](0)*_grad_phie[_qp](0)+_grad_phie[_qp](1)*_grad_phie[_qp](1)) *
        0.5 * (_epsilon_water - _epsilon_air) * (_grad_phi[_j][_qp](_component)) * _test[_i][_qp];
   }
  
//  if (jvar == _phie_var)
//  {
// return 0.5 * _epsilon_zero * (_grad_phi[_j][_qp](0)*_grad_phi[_j][_qp](0)+_grad_phi[_j][_qp](1)*_grad_phi[_j][_qp](1)) *
//        0.5 * (_epsilon_water - _epsilon_air) * (_grad_c[_qp](_component)) * _test[_i][_qp];
 //  }

    return 0.;
}
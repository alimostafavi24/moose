

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Abels.h"
#include "NS.h"
#include "MooseVariableFE.h"
#include "SystemBase.h"

registerMooseObject("NavierStokesApp", Abels);

InputParameters
Abels::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Adds the surface tension force to the NS equation.");  
 params.addParam<Real>("rho2", 1.0, "rho2");
  params.addParam<Real>("rho1", 1.0, "rho1");
  params.addParam<Real>("mobility", 1.0, "mobility");

//  params.addRequiredCoupledVar("variable", "The order parameter");
  params.addRequiredCoupledVar("w", "The chemical potential");
  //params.addRequiredParam<unsigned int>("component", "number of component (0 = x, 1 = y, 2 = z)");

  return params;
}

Abels::Abels(const InputParameters & parameters)
  : Kernel(parameters),
    _rho2(getParam<Real>("rho2")),
    _rho1(getParam<Real>("rho1")),
    _mobility(getParam<Real>("mobility")),        
  //  _u(coupledValue("variable")),
   // _grad_u(coupledGradient("variable")),
    _grad_w(coupledGradient("w")),


    //_component(getParam<unsigned int>("component")), 
   // _vel_x_var(coupled("c")),
    _w_var(coupled("w"))
{
}

Real
Abels::computeQpResidual()
{
 //return -_c[_qp] * _w[_qp] * _grad_test[_i][_qp](_component);
  return (_rho2 - _rho1) / 2 * (_mobility) * (_grad_w[_qp](0) * _grad_u[_qp](0) + _grad_w[_qp](1) * _grad_u[_qp](1)) ;


}

Real
Abels::computeQpJacobian()
{
  return (_rho2 - _rho1) / 2 * (_mobility) * (_grad_w[_qp](0) * _grad_phi[_j][_qp](0) + _grad_w[_qp](1) * _grad_phi[_j][_qp](1)) ;
}

Real
Abels::computeQpOffDiagJacobian(unsigned jvar)
{

   if (jvar == _w_var)
  {
    return (_rho2 - _rho1) / 2 * (_mobility) * (_grad_phi[_j][_qp](0) * _grad_u[_qp](0) + _grad_phi[_j][_qp](1) * _grad_u[_qp](1)) ;
  }
    return 0.;
}
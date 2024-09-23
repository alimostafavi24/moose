//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WettingBCEnergy.h"
#include "TimeDerivative.h"

registerMooseObject("MooseApp", WettingBCEnergy);

InputParameters
WettingBCEnergy::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription(
      "Computes a boundary residual contribution consistent with the Diffusion Kernel. "
      "Does not impose a boundary condition; instead computes the boundary "
      "contribution corresponding to the current value of grad(u) and accumulates "
      "it in the residual vector.");
  //params.addRequiredCoupledVar("variable", "The name of the variable that this object operates on");
  //  params.addRequiredCoupledVar("c", "c");
  params.addRequiredCoupledVar("vel_x", "vel_x");
  params.addRequiredCoupledVar("vel_y", "vel_y");
   params.addParam<Real>("cos", 0.0, "cos");
   params.addParam<Real>("sigma", 1.0, "sigma");
   params.addParam<Real>("gamma", 1.0, "gamma");
   params.addParam<Real>("eta", 1.0, "eta");

  return params;
}

WettingBCEnergy::WettingBCEnergy(const InputParameters & parameters) : IntegratedBC(parameters),
   // _u(coupledValue("variable")),
   // _grad_u(coupledGradient("variable")),
     _u_dot(_var.uDot()),
    _du_dot_du(_var.duDotDu()),
    _vel_x(coupledValue("vel_x")),
    _vel_y(coupledValue("vel_x")),    
    _cos(getParam<Real>("cos")),
    _sigma(getParam<Real>("sigma")),
    _gamma(getParam<Real>("gamma")),
    _eta(getParam<Real>("eta"))
 {
}

Real
WettingBCEnergy::computeQpResidual()
{
  return (_test[_i][_qp] * _u_dot[_qp]) + 

       ((_vel_x [_qp] *  _grad_u[_qp](0)) + (_vel_y[_qp] * _grad_u[_qp](1) )) * 
       _test[_i][_qp] + 

       _eta * ( _sigma * _grad_u[_qp] * _normals[_qp] -
        _cos * 0.75 * _gamma * (1 - _u[_qp] * _u[_qp])) * _test[_i][_qp];

}

Real
WettingBCEnergy::computeQpJacobian()
{
return (_test[_i][_qp] * _phi[_j][_qp] * _du_dot_du[_qp]) + 

       (_vel_x [_qp] * _grad_phi[_j][_qp](0) + _vel_y [_qp] * _grad_phi[_j][_qp](0)) *
       _test[_i][_qp] +

        _eta * ( _sigma * _grad_phi[_j][_qp] * _normals[_qp] -
        _cos * 0.75 * _gamma * (1 - 2 * _u[_qp] * _phi[_j][_qp])) * _test[_i][_qp] ;

       ;

}
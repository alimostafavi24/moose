//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WettingBC3D.h"

registerMooseObject("MooseApp", WettingBC3D);

InputParameters
WettingBC3D::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription(
      "Computes a boundary residual contribution consistent with the Diffusion Kernel. "
      "Does not impose a boundary condition; instead computes the boundary "
      "contribution corresponding to the current value of grad(u) and accumulates "
      "it in the residual vector.");
   params.addParam<Real>("cos", 0.0, "cos");


  return params;
}

WettingBC3D::WettingBC3D(const InputParameters & parameters) : IntegratedBC(parameters),
       _cos(getParam<Real>("cos"))
 {
}

Real
WettingBC3D::computeQpResidual()
{
  return (_grad_u[_qp] * _normals[_qp]) * _test[_i][_qp] - sqrt((_grad_u[_qp](0)*_grad_u[_qp](0)+_grad_u[_qp](1)*_grad_u[_qp](1)+_grad_u[_qp](2)*_grad_u[_qp](2)))*_cos* _test[_i][_qp];
}

Real
WettingBC3D::computeQpJacobian()
{
  return (_grad_phi[_j][_qp] * _normals[_qp]) * _test[_i][_qp] - sqrt((_grad_phi[_j][_qp](0)*_grad_phi[_j][_qp](0)+_grad_phi[_j][_qp](1)*_grad_phi[_j][_qp](1)+_grad_phi[_j][_qp](2)*_grad_phi[_j][_qp](2)))*_cos* _test[_i][_qp];
}
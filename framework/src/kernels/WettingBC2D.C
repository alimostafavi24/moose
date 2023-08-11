//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WettingBC2D.h"

registerMooseObject("MooseApp", WettingBC2D);

InputParameters
WettingBC2D::validParams()
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

WettingBC2D::WettingBC2D(const InputParameters & parameters) : IntegratedBC(parameters),
       _cos(getParam<Real>("cos"))
 {
}

Real
WettingBC2D::computeQpResidual()
{
  return (_grad_u[_qp] * _normals[_qp]) * _test[_i][_qp] - sqrt((_grad_u[_qp](0)*_grad_u[_qp](0)+_grad_u[_qp](1)*_grad_u[_qp](1)))*_cos* _test[_i][_qp];
}

Real
WettingBC2D::computeQpJacobian()
{
  return (_grad_phi[_j][_qp] * _normals[_qp]) * _test[_i][_qp] - sqrt((_grad_phi[_j][_qp](0)*_grad_phi[_j][_qp](0)+_grad_phi[_j][_qp](1)*_grad_phi[_j][_qp](1)))*_cos* _test[_i][_qp];
}
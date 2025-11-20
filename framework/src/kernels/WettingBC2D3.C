
#include "WettingBC2D3.h"

registerMooseObject("MooseApp", WettingBC2D3);

InputParameters
WettingBC2D3::validParams()
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

WettingBC2D3::WettingBC2D3(const InputParameters & parameters) : IntegratedBC(parameters),
       _cos(getParam<Real>("cos"))
 {
}

Real
WettingBC2D3::computeQpResidual()
{
  return (_grad_u[_qp] * _normals[_qp]) * _test[_i][_qp] 
  - sqrt((_grad_u[_qp]*_grad_u[_qp]))*_cos* _test[_i][_qp];
}

Real
WettingBC2D3::computeQpJacobian()
{
  //return (_grad_phi[_j][_qp] * _normals[_qp]) * _test[_i][_qp] - sqrt((_grad_phi[_j][_qp](0)*_grad_phi[_j][_qp](0)+_grad_phi[_j][_qp](1)*_grad_phi[_j][_qp](1)))*_cos* _test[_i][_qp];

   return (_grad_phi[_j][_qp] * _normals[_qp]) * _test[_i][_qp] - 
1/( sqrt(_grad_u[_qp]*_grad_u[_qp])+std::numeric_limits<double>::epsilon())*_grad_u[_qp] * _grad_phi[_j][_qp] * _cos* _test[_i][_qp];
}





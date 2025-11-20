

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LatentHeatMultiphaseTime.h"


registerMooseObject("HeatTransferApp", LatentHeatMultiphaseTime);

InputParameters
LatentHeatMultiphaseTime::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addClassDescription("Adds transient latent heat.");
  params.addRequiredCoupledVar("c", "The order parameter");  
  params.addRequiredParam<MaterialPropertyName>("dfldT", "dfldT");
  params.addRequiredParam<Real>("rho_l", "rho_l");
  params.addRequiredParam<Real>("latent_heat", "latent_heat");



  return params;
}

LatentHeatMultiphaseTime::LatentHeatMultiphaseTime(const InputParameters & parameters)
  : TimeDerivative(parameters),
  _c(coupledValue("c")),
 
    _dfldT(getMaterialProperty<Real>("dfldT")),

    _rho_l(getParam<Real>("rho_l")),
    _latent_heat(getParam<Real>("latent_heat"))

{
}

Real
LatentHeatMultiphaseTime::computeQpResidual()
{
 //return -_c[_qp] * _w[_qp] * _grad_test[_i][_qp](_component);
 return _dfldT[_qp] * _rho_l * _latent_heat * TimeDerivative::computeQpResidual() * 0.5 * (1 + _c[_qp]);

}

Real
LatentHeatMultiphaseTime::computeQpJacobian()
{
  return _dfldT[_qp] * _rho_l * _latent_heat * TimeDerivative::computeQpJacobian() * 0.5 * (1 + _c[_qp]);
}


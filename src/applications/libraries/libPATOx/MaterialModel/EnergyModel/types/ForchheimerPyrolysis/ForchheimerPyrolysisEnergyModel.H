/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If Pyrolysist, see <http://www.gnu.org/licenses/>.

Class
    Foam::ForchheimerPyrolysisEnergyModel

Description

SourceFiles
    ForchheimerPyrolysisEnergyModel.H
    ForchheimerPyrolysisEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef ForchheimerPyrolysisEnergyModel_H
#define ForchheimerPyrolysisEnergyModel_H

#include "simpleEnergyModel.H"
#include "simpleMaterialPropertiesModel.H"
#include "simpleGasPropertiesModel.H"
#include "simplePyrolysisModel.H"
#include "simpleMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ForchheimerPyrolysisEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class ForchheimerPyrolysisEnergyModel
:
public simpleEnergyModel
{

 protected:

  //- Heat capacity [J/kg/K]
  const volScalarField& cp;

  //- Solid thermal conductivity [W/m/K]
  const volTensorField& k;

  //- Pyrolysis flux [W/m3] - updated in Material Properties Model
  const volScalarField& pyrolysisFlux_;

  //- Solid density [kg/m3]
  const volScalarField& rho_s;

  //- Temperature [K]
  volScalarField& T;

  //- Universal gas constant [J/mol/K]
  const dimensionedScalar R;

  //- Mass Model
  const simpleMassModel& massModel_;

  //- Forchheimer's coefficient [1/m]
  const volTensorField& Beta;

  //- Gravity field [m/s2]
  const uniformDimensionedVectorField& g;

  //- Permeability [m2]
  const volTensorField& K;

  //- Pressure [Pa]
  const volScalarField& p;

  //- Gaseous velocity times gas volume fraction [m/s]
  const volVectorField& U;

  //- Gas Properties Model
  const simpleGasPropertiesModel& gasPropertiesModel_;

  //- Gas volume fraction (1 - porosity) [-]
  const volScalarField& eps_g;

  //- Gas enthalpy [J/kg]
  const volScalarField& h_g;

  //- Gas density [kg/m3]
  const volScalarField& rho_g;

  //- Molar Mass [kg/kmol]
  const volScalarField& M;

  //- Viscosity [kg/m/s]
  const volScalarField& mu;

  //- Gas volume fraction [-] * Gas density [kg/m3] * Total gas energy [J/kg]
  volScalarField& epsgRhogEg;

  //- Convection term / grad(pressure) - gravitational term
  volTensorField& GammaHg;

  //- Convection term (gravitational part)
  volTensorField& Gamma_GHg;

  //- Interpolation on the face of the convection term coming from gravity field [kg/m3/s]
  surfaceScalarField& phiHg;

 public:

  //- Runtime type information
  TypeName("ForchheimerPyrolysis");

  // Constructors

  //- Construct from mesh and phase name
  ForchheimerPyrolysisEnergyModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~ForchheimerPyrolysisEnergyModel();

  //- Member Functions

  //- Update the energy equation fields
  virtual void update();

  //- Before to solve the energy equation
  void beforeSolve();

  //- After to solve the energy equation
  void afterSolve();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

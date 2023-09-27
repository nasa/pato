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
    Foam::PyrolysisEnergyModel

Description

SourceFiles
    PyrolysisEnergyModel.H
    PyrolysisEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef PyrolysisEnergyModel_H
#define PyrolysisEnergyModel_H

#include "simpleEnergyModel.H"
#include "simpleGasPropertiesModel.H"
#include "simpleMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class PyrolysisEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class PyrolysisEnergyModel
:
public simpleEnergyModel
{

 protected:

  //- Reference to the mesh database
  const fvMesh& mesh_;

  //- Reference to the region name
  const word& dictName_;

  //- Temperature [K]
  volScalarField& T;

  //- Solid density [kg/m3]
  const volScalarField& rho_s;

  //- Solid thermal conductivity [W/m/K]
  const volTensorField& k;

  //- Heat capacity [J/kg/K]
  const volScalarField& cp;

  //- Pyrolysis flux [W/m3] - updated in Material Properties Model
  const volScalarField& pyrolysisFlux_;

  //- Gas Properties Model
  const simpleGasPropertiesModel& gasPropertiesModel_;

  //- Gas enthalpy [J/kg]
  const volScalarField& h_g;

  //- Molar Mass [kg/kmol]
  const volScalarField& M;

  //- Viscosity [kg/m/s]
  const volScalarField& mu;

  //- Gas volume fraction (1 - porosity) [-]
  const volScalarField& eps_g;

  //- Gas density [kg/m3]
  const volScalarField& rho_g;

  //- Mass Properties Model
  const simpleMassModel& massModel_;

  //- Pressure [Pa]
  const volScalarField& p;

  //- Gaseous velocity times gas volume fraction [m/s]
  const volVectorField& U;

  //- Permeability [m2]
  const volTensorField& K;

  //- Gravity field [m/s2]
  const uniformDimensionedVectorField& g;

  //- Universal gas constant [J/mol/K]
  const dimensionedScalar R;

  //- Gas volume fraction [-] * Gas density [kg/m3] * Total gas energy [J/kg]
  volScalarField& epsgRhogEg;

  //- Convection term / grad(pressure)
  volTensorField& GammaHg;

  //- Convection term (gravitational part)
  volTensorField& Gamma_GHg;

  //- Interpolation on the face of the convection term coming from gravity field [kg/m3/s]
  surfaceScalarField& phiHg;

 public:

  //- Runtime type information
  TypeName("Pyrolysis");

  // Constructors

  //- Construct from mesh and phase name
  PyrolysisEnergyModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~PyrolysisEnergyModel();

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

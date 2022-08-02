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
#include "simpleMaterialPropertiesModel.H"
#include "simpleGasPropertiesModel.H"
#include "simplePyrolysisModel.H"

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

  //- Material Properties Model
  simpleMaterialPropertiesModel& materialPropertiesModel_;

  //- Gas Properties Model
  simpleGasPropertiesModel& gasPropertiesModel_;

  //- Pyrolysis Model
  simplePyrolysisModel& pyrolysisModel_;

  //- Gas enthalpy [J/kg]
  volScalarField& h_g;

  //- Pressure [Pa]
  volScalarField& p;

  //- Molar Mass [kg/kmol]
  volScalarField& M;

  //- Viscosity [kg/m/s]
  volScalarField& mu;

  //- Universal gas constant [J/mol/K]
  const dimensionedScalar R;

  //- Temperature [K]
  volScalarField& T;

  //- Permeability [m2]
  const volTensorField& K;

  //- Gas volume fraction (1 - porosity) [-]
  const volScalarField& eps_g;

  //- Gas density [kg/m3]
  const volScalarField& rho_g;

  //- Solid density [kg/m3]
  volScalarField& rho_s;

  //- Heat capacity [J/kg/K]
  volScalarField& cp;

  //- Pyrolysis flux [W/m3] - updated in Material Properties Model
  volScalarField& pyrolysisFlux_;

  //- Gas volume fraction [-] * Gas density [kg/m3] * Total gas energy [J/kg]
  volScalarField epsgRhogEg;

  //- Convection term / grad(pressure)
  volTensorField GammaHg;

  //- Solid thermal conductivity [W/m/K]
  volTensorField& k;

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

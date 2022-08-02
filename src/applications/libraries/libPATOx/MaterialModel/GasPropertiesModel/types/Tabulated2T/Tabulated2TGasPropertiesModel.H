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
    along with OpenFOAM.  If Tabulated2Tt, see <http://www.gnu.org/licenses/>.

Class
    Foam::Tabulated2TGasPropertiesModel

Description

SourceFiles
    Tabulated2TGasPropertiesModel.H
    Tabulated2TGasPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef Tabulated2TGasPropertiesModel_H
#define Tabulated2TGasPropertiesModel_H

#include "simpleEnergyModel.H"
#include "simpleGasPropertiesModel.H"
#include "Tabulated2TGasPropertiesObject.H"
#include "simplePyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class Tabulated2TGasPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class Tabulated2TGasPropertiesModel
:
public simpleGasPropertiesModel
{

 protected:

  //- EnergyModel
  simpleEnergyModel& energyModel_;

  //- Average molar mass of the gas phase [kg/mol]
  volScalarField& M;

  //- Gas heat capacity [J/kg/K]
  volScalarField& cp_g;

  //- Average viscosity of the gas phase [kg/m/s]
  volScalarField& mu;

  //- Average enthalpy of the gas phase [J/kg]
  volScalarField& h_g;

  //- Average volume fraction of the gas phase [kg/m3]
  volScalarField& eps_g;

  //- Average density of the gas phase [kg/m3]
  volScalarField& rho_g;

  //- Average thermal conductivity of the gas phase [kg/m3]
  volScalarField& k_g;

  //- Average temperature of the gas phase [K]
  volScalarField& Tg;

  //- Average pressure of the gas phase [Pa]
  volScalarField& p;

  //- Tabulated2TGasPropertiesSolverObject which handles the gas property lookup and interpolations when called
  autoPtr<Tabulated2TGasPropertiesObject> gasPT_ptr;

  //- Table directory
  const fileName gasPropertiesFile_;

  //- Material properties directory
  const fileName materialPropertiesDirectory;

  //- Constant properties dictionary
  IOdictionary constantPropertiesDictionary;

  //- Charred gaseous volume fraction [-]
  dimensionedScalar eps_g_c_;

  //- Virgin gaseous volume fraction [-]
  dimensionedScalar eps_g_v_;

  //- Pyrolysis Model
  simplePyrolysisModel& pyrolysisModel_;

  //- Advancement of the pyrolysis
  volScalarField& tau_;

 public:

  //- Runtime type information
  TypeName("Tabulated2T");

  // Constructors

  //- Construct from mesh and phase name
  Tabulated2TGasPropertiesModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~Tabulated2TGasPropertiesModel();

  //- Member Functions

  //- Update
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef Tabulated2TRepository
#include "Tabulated2TGasPropertiesModel.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

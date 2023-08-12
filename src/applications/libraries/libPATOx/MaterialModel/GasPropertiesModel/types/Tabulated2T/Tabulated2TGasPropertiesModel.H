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


#include "Tabulated2TGasPropertiesObject.H"
#include "simpleEnergyModel.H"
#include "simpleMassModel.H"
#include "simplePyrolysisModel.H"
#include "simpleGasPropertiesModel.H"

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

  //- Gas heat capacity [J/kg/K]
  volScalarField& cp_g;

  //- Average volume fraction of the gas phase [kg/m3]
  volScalarField& eps_g;

  //- Average enthalpy of the gas phase [J/kg]
  volScalarField& h_g;

  //- Average thermal conductivity of the gas phase [kg/m3]
  volScalarField& k_g;

  //- Average molar mass of the gas phase [kg/mol]
  volScalarField& M;

  //- Average viscosity of the gas phase [kg/m/s]
  volScalarField& mu;

  //- Average density of the gas phase [kg/m3]
  volScalarField& rho_g;

  //- Charred gaseous volume fraction [-]
  const dimensionedScalar& eps_g_c_;

  //- Virgin gaseous volume fraction [-]
  const dimensionedScalar& eps_g_v_;

  //- Energy Model
  const simpleEnergyModel& energyModel_;

  //- Average temperature of the gas phase [K]
  const volScalarField& Tg;

  //- Mass Model
  const simpleMassModel& massModel_;

  //- Average pressure of the gas phase [Pa]
  const volScalarField& p;

  //- Pyrolysis Model
  simplePyrolysisModel& pyrolysisModel_;

  //- Advancement of the pyrolysis
  volScalarField& tau_;

  //- Table directory
  const fileName gasPropertiesFile_;

  //- Tabulated2TGasPropertiesSolverObject which handles the gas property lookup and interpolations when called
  autoPtr<Tabulated2TGasPropertiesObject> gasPT_ptr;

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

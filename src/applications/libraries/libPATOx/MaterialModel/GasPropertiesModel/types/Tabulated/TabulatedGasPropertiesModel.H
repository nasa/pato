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
    along with OpenFOAM.  If Tabulatedt, see <http://www.gnu.org/licenses/>.

Class
    Foam::TabulatedGasPropertiesModel

Description

SourceFiles
    TabulatedGasPropertiesModel.H
    TabulatedGasPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef TabulatedGasPropertiesModel_H
#define TabulatedGasPropertiesModel_H

#include "simpleGasPropertiesModel.H"
#include "TabulatedGasPropertiesObject.H"
#include "simplePyrolysisModel.H"
#include "simpleEnergyModel.H"
#include "simpleMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class TabulatedGasPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class TabulatedGasPropertiesModel
:
public simpleGasPropertiesModel
{

 protected:

  // Created fields

  //- Average molar mass of the gas phase [kg/mol]
  volScalarField& M;

  //- Average viscosity of the gas phase [kg/m/s]
  volScalarField& mu;

  //- Average enthalpy of the gas phase [J/kg]
  volScalarField& h_g;

  //- Average volume fraction of the gas phase [kg/m3]
  volScalarField& eps_g;

  //- Average density of the gas phase [kg/m3]
  volScalarField& rho_g;

  // Reference fields

  //- Pyrolysis Model
  const simplePyrolysisModel& pyrolysisModel;

  //- Advancement of the pyrolysis
  const volScalarField& tau_;

  //- Energy Model
  const simpleEnergyModel& energyModel;

  //- Average temperature of the gas phase [K]
  const volScalarField& Tg;

  //- Mass Model
  const simpleMassModel& massModel;

  //- Average pressure of the gas phase [Pa]
  const volScalarField& p;

  //- TabulatedGasPropertiesSolverObject which handles the gas property lookup and interpolations when called
  autoPtr<TabulatedGasPropertiesObject> gasPT_ptr;

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

 public:

  //- Runtime type information
  TypeName("Tabulated");

  // Constructors

  //- Construct from mesh and phase name
  TabulatedGasPropertiesModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~TabulatedGasPropertiesModel();

  //- Member Functions

  //- Update
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef TabulatedRepository
#include "TabulatedGasPropertiesModel.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

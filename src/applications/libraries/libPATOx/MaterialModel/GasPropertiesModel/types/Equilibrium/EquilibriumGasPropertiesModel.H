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
    along with OpenFOAM.  If Equilibriumt, see <http://www.gnu.org/licenses/>.

Class
    Foam::EquilibriumGasPropertiesModel

Description

SourceFiles
    EquilibriumGasPropertiesModel.H
    EquilibriumGasPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef EquilibriumGasPropertiesModel_H
#define EquilibriumGasPropertiesModel_H

#include "simpleGasPropertiesModel.H"
#include "simplePyrolysisModel.H"
#include "simpleMaterialChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class EquilibriumGasPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class EquilibriumGasPropertiesModel
:
public simpleGasPropertiesModel
{

 protected:

  //- MaterialChemistry Model correct flag
  Switch verifyMaterialChemistry_;

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

  //- Average temperature of the gas phase [K]
  volScalarField& Tg;

  //- Average pressure of the gas phase [Pa]
  volScalarField& p;

  //- Previous average temperature of the gas phase [K]
  volScalarField Tg_old;

  //- Previous average pressure of the gas phase [Pa]
  volScalarField p_old;

  //- Temperature threshold to update
  scalar T_threshold_;

  //- Pressure threshold to update
  scalar p_threshold_;

  //- Material properties directory
  const fileName materialPropertiesDirectory;

  //- Material properties dictionary
  IOdictionary constantPropertiesDictionary;

  //- Charred volume fraction [-]
  dimensionedScalar eps_g_c_;

  //- Virgin volume fraction [-]
  dimensionedScalar eps_g_v_;

  //- Pyrolysis Model
  simplePyrolysisModel& pyrolysisModel_;

  //- Advancement of the pyrolysis
  volScalarField& tau_;

  //- MaterialChemistry Model
  simpleMaterialChemistryModel& MaterialChemistryModel_;

  //- mixture
  autoPtr<Mutation::Mixture>& mix_;

  //- pressure and temperature for Mutation++
  double* pTp;

  // mass fractions of the species/elements
  PtrList<volScalarField>& massFractions_;

  // mole fractions of the species/elements
  PtrList<volScalarField>& moleFractions_;

  //- mole fractions of the species/elements for Mutation++
  double* p_Zx;

 public:

  //- Runtime type information
  TypeName("Equilibrium");

  // Constructors

  //- Construct from mesh and phase name
  EquilibriumGasPropertiesModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~EquilibriumGasPropertiesModel();

  //- Member Functions

  //- Update
  virtual void update();

  //- Verify the material chemistry type
  inline Switch verifyMaterialChemistry();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

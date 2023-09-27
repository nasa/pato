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


#include "simpleMaterialChemistryModel.H"
#include "simpleEnergyModel.H"
#include "simpleMassModel.H"
#include "simpleGasPropertiesModel.H"
#include "simplePyrolysisModel.H"
#include <fcntl.h>
#include <unistd.h>

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

  //- Average volume fraction of the gas phase [kg/m3]
  volScalarField& eps_g;

  //- Average enthalpy of the gas phase [J/kg]
  volScalarField& h_g;

  //- Average molar mass of the gas phase [kg/mol]
  volScalarField& M;

  //- Average viscosity of the gas phase [kg/m/s]
  volScalarField& mu;

  //- Average density of the gas phase [kg/m3]
  volScalarField& rho_g;

  //- Charred volume fraction [-]
  const dimensionedScalar& eps_g_c_;

  //- Virgin volume fraction [-]
  const dimensionedScalar& eps_g_v_;

  //- Pressure threshold to update
  const scalar& p_threshold_;

  //- Temperature threshold to update
  const scalar& T_threshold_;

  //- Energy Model
  const simpleEnergyModel& energyModel_;

  //- Average temperature of the gas phase [K]
  const volScalarField& Tg;

  //- Mass Model
  const simpleMassModel& massModel_;

  //- Average pressure of the gas phase [Pa]
  const volScalarField& p;

  //- Pyrolysis Model
  const simplePyrolysisModel& pyrolysisModel_;

  //- Advancement of the pyrolysis
  const volScalarField& tau_;

  //- MaterialChemistry Model
  const simpleMaterialChemistryModel& MaterialChemistryModel_;

  //- mixture
  autoPtr<Mutation::Mixture>& mix_;

  //- pressure and temperature for Mutation++
  double* pTp;

  // mass fractions of the species/elements
  const PList<volScalarField>& massFractions_;

  // mole fractions of the species/elements
  const PList<volScalarField>& moleFractions_;

  //- mole fractions of the species/elements for Mutation++
  double* p_Zx;

  //- Previous average temperature of the gas phase [K]
  volScalarField& Tg_old;

  //- Previous average pressure of the gas phase [Pa]
  volScalarField& p_old;

  //- Initialized flag
  Switch initialized_;

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

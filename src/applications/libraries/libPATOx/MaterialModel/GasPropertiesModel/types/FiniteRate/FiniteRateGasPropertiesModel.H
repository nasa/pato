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
    along with OpenFOAM.  If FiniteRatet, see <http://www.gnu.org/licenses/>.

Class
    Foam::FiniteRateGasPropertiesModel

Description

SourceFiles
    FiniteRateGasPropertiesModel.H
    FiniteRateGasPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef FiniteRateGasPropertiesModel_H
#define FiniteRateGasPropertiesModel_H


#include "simpleMaterialChemistryModel.H"
#include "simpleEnergyModel.H"
#include "simpleMassModel.H"
#include "simpleGasPropertiesModel.H"
#include "simplePyrolysisModel.H"
#include "IOFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class FiniteRateGasPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class FiniteRateGasPropertiesModel
:
public simpleGasPropertiesModel
{

 protected:



  //- MaterialChemistry Model correct flag
  Switch verifyMaterialChemistry_;

  //- Diffusion coefficient [m2/s]
  PList<volScalarField> Dm_;

  //- Species enthalpy [J/kg]
  PList<volScalarField> hSpeciesField_;

  //- Flux transported by diffusion of the molecules [W/m3]
  volScalarField& Ediff;

  //- Average volume fraction of the gas phase [kg/m3]
  volScalarField& eps_g;

  //- Average molar mass of the gas phase [kg/mol]
  volScalarField& M;

  //- Average viscosity of the gas phase [kg/m/s]
  volScalarField& mu;

  //- Average enthalpy of the gas phase [J/kg]
  volScalarField& h_g;

  //- Average density of the gas phase [kg/m3]
  volScalarField& rho_g;

  //- Char porosity [-]
  const dimensionedScalar& eps_g_c_;

  //- Virgin porosity [-]
  const dimensionedScalar& eps_g_v_;

  //- Tortuosity [-]
  const dimensionedScalar& eta0_;

  //- Pressure threshold to update
  const scalar& p_threshold_;

  //- Temperature threshold to update
  const scalar& T_threshold_;

  //- Ideal gas constant [J/mol/K]
  const dimensionedScalar R;

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
  const simpleMaterialChemistryModel& materialChemistryModel_;

  //- mixture
  autoPtr<Mutation::Mixture>& mix_;

  //- Number of species
  int ns_mix;

  //- Pressure and temperature for Mutation++
  double* pTp;

  //- Mass fractions of the species [-]
  const PList<volScalarField>& massFractions_;

  //- Name of the species
  const wordList& speciesNames_;

  //- Index species in Mutation++
  const labelList& speciesIndexMutation_;

  //- Mass fractions of the species for Mutation++ [-]
  double* p_Y;

  //- Diffusion coefficient of the species for Mutation++ [m2/s]
  double* p_Dm;

  //- Enthalpy of the species for Mutation++ [J/kg]
  double* p_hSpecies;

  //- Previous average temperature of the gas phase [K]
  volScalarField& Tg_old;

  //- Previous average pressure of the gas phase [Pa]
  volScalarField& p_old;

  //- Initialized flag
  Switch initialized_;

 public:

  //- Runtime type information
  TypeName("FiniteRate");

  // Constructors

  //- Construct from mesh and phase name
  FiniteRateGasPropertiesModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~FiniteRateGasPropertiesModel();

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

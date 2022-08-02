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

#include "simpleGasPropertiesModel.H"
#include "simplePyrolysisModel.H"
#include "simpleMaterialChemistryModel.H"
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

  //- Ideal gas constant [J/mol/K]
  const dimensionedScalar R;

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

  //- Flux transported by diffusion of the molecules [W/m3]
  volScalarField& Ediff;

  //- Diffusion coefficient [m2/s]
  PtrList<volScalarField>& Dm_;

  //- Species enthalpy [J/kg]
  PtrList<volScalarField> hSpeciesField_;

  //- Previous average temperature of the gas phase [K]
  volScalarField Tg_old;

  //- Previous average pressure of the gas phase [Pa]
  volScalarField p_old;

  //- Temperature threshold to update
  scalar T_threshold_;

  //- Pressure threshold to update
  scalar p_threshold_;

  //- Material Properties Directory
  const fileName materialPropertiesDirectory;

  //- Constant Properties Dictionary
  IOdictionary constantPropertiesDictionary;

  //- Tortuosity [-]
  dimensionedScalar eta0_;

  //- Char porosity [-]
  dimensionedScalar eps_g_c_;

  //- Virgin porosity [-]
  dimensionedScalar eps_g_v_;

  //- Pyrolysis Model
  simplePyrolysisModel& pyrolysisModel_;

  //- Advancement of the pyrolysis
  volScalarField& tau_;

  //- MaterialChemistry Model
  simpleMaterialChemistryModel& MaterialChemistryModel_;

  //- mixture
  autoPtr<Mutation::Mixture>& mix_;

  //- Number of species
  int ns_mix;

  //- Pressure and temperature for Mutation++
  double* pTp;

  //- Mass fractions of the species [-]
  PtrList<volScalarField>& massFractions_;

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

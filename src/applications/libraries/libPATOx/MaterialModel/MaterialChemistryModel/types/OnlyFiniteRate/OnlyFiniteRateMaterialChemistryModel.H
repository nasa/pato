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
    along with OpenFOAM.  If OnlyFiniteRatet, see <http://www.gnu.org/licenses/>.

Class
    Foam::OnlyFiniteRateMaterialChemistryModel

Description

SourceFiles
    OnlyFiniteRateMaterialChemistryModel.H
    OnlyFiniteRateMaterialChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef OnlyFiniteRateMaterialChemistryModel_H
#define OnlyFiniteRateMaterialChemistryModel_H

#include "simpleMaterialChemistryModel.H"
#include "mathFunctions.H"
#include "IOFunctions.H"
#include "simpleGasPropertiesModel.H"
#include "psiReactionThermo.H"
#include "BasicFiniteRateChemistryModel.H"
#include "OFstream.H"
#include "thermoPhysicsTypes.H"
#include "basicMultiComponentMixture.H"
#include "coupledMixedFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class OnlyFiniteRateMaterialChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class OnlyFiniteRateMaterialChemistryModel
:
public simpleMaterialChemistryModel
{

 protected:

  //- Mutation++ Mixture name
  word mixtureMutation;

  //- Mutation++ Mixture Options
  autoPtr<Mutation::MixtureOptions> optEq_;

  //- Mutation++ Mixture
  autoPtr<Mutation::Mixture>& mix_;

  //- Number of species in the mixture
  int ns_mix;

  //- List of the element names
  wordList& elementNames_;

  //- List of the species names
  wordList& speciesNames_;

  //- Species mass fractions
  PtrList<volScalarField>& massFractions_;

  //- Species mole fractions
  PtrList<volScalarField>& moleFractions_;

  //- Temperature (K)
  volScalarField& T_;

  //- Pressure (Pa)
  volScalarField& p_;

  //- Default species mass/mole fractions
  volScalarField Ydefault_;

  //- Surface Interpolation Scheme table
  multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

  //- Psi reaction thermo model
  autoPtr<psiReactionThermo> pThermo_;

  //- Psi reaction thermo model reference
  psiReactionThermo& thermo_;

  //- Chemistry model coupled to psithermo model
  autoPtr<BasicFiniteRateChemistryModel<psiReactionThermo>> pChemistry_;

  //- Referemce of the chemistry model coupled to psithermo model
  BasicFiniteRateChemistryModel<psiReactionThermo>& chemistry_;

  //- Latest estimation of MaterialChemistry integration
  scalar dtChem_;

  //- CHEMKIN composition
  basicSpecieMixture& composition_;

  //- CHEMKIN species mass fractions
  PtrList<volScalarField>& Y_;

  //- Index species in Mutation++
  labelList& speciesIndexMutation_;

  //- Init flag for mixture
  Switch init_;

  //- Gas Properties model
  simpleGasPropertiesModel& gasPropertiesModel_;

  //- Gas density [kg/m3]
  volScalarField& rho_g_;

  //- Gas volume fraction [-]
  volScalarField& eps_g_;

  //- Species diffusion coefficient
  PtrList<volScalarField>& Dm_;

  //- Material properties directory
  const fileName materialPropertiesDirectory;

  //- ConstantProperties dictionary
  IOdictionary constantPropertiesDictionary;

  //- Tortuosity
  dimensionedScalar eta0_;

  //- Diffusive part []
  PtrList<volScalarField> diffY_;

 public:

  //- Runtime type information
  TypeName("OnlyFiniteRate");

  // Constructors

  //- Construct from mesh and phase name
  OnlyFiniteRateMaterialChemistryModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~OnlyFiniteRateMaterialChemistryModel();

  //- Member Functions

  //- Update the finite rate fields
  virtual void update();

  //- Initialize the mixture
  inline Switch init();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

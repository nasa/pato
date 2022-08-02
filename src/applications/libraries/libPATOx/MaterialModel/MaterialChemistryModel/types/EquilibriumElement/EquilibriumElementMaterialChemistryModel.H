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
    along with OpenFOAM.  If EquilibriumElementt, see <http://www.gnu.org/licenses/>.

Class
    Foam::EquilibriumElementMaterialChemistryModel

Description

SourceFiles
    EquilibriumElementMaterialChemistryModel.H
    EquilibriumElementMaterialChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef EquilibriumElementMaterialChemistryModel_H
#define EquilibriumElementMaterialChemistryModel_H

#include "simpleMaterialChemistryModel.H"
#include "mathFunctions.H"
#include "BprimeFvPatchScalarField.H"
#include "IOFunctions.H"
#include "simpleGasPropertiesModel.H"
#include "simpleMassModel.H"
#include "simplePyrolysisModel.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class EquilibriumElementMaterialChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class EquilibriumElementMaterialChemistryModel
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

  //- Number of elements in the mixture
  int ne_mix;

  //- Initial mass fraction of elements
  scalarList p_Z;

  //- Initial reference mass fraction of elements
  scalarList p_Z0;

  //- Initial mole fraction of elements
  scalarList p_Zx;     // mole fractions of elements

  //- Minimum tolerence on the mass and mole fractions
  const scalar tolZ;

  // General

  //- Material properties directory
  const fileName materialPropertiesDirectory;

  //- ConstantProperties dictionary
  IOdictionary constantPropertiesDictionary;

  //- List of the element names
  wordList& elementNames_;

  //- List of the species names
  wordList& speciesNames_;

  //- Elemental mass fractions
  PtrList<volScalarField>& massFractions_;

  //- Elemental mole fractions
  PtrList<volScalarField>& moleFractions_;

  //- Mole fractions in the material from the user if not found in constantProperties
  List<Tuple2<word,scalar> > moleFractionGasInMaterial_;

  //- Elemental mass flux
  PtrList<surfaceVectorField> ZmassFlux_;

  //- Diffusion coefficient in binary diffusion mode
  dimensionedScalar D0_;

  //- Tortuosity
  dimensionedScalar eta0_;

  //- Surface Mass and Energy Balance (SMEB) patches - BoundaryConditions Model
  labelList BprimePatches_;

  //- Bprime BC for the SMEB patches
  List<BprimeBoundaryConditions*> BprimeBC_;

  //- Temperature (K)
  volScalarField& T_;

  //- Default elemental mass/mole fractions
  volScalarField Zdefault_;

  // Surface Interpolation Scheme table
  multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

  //- Init flag for mixture
  Switch init_;

  //- Gas Properties model
  simpleGasPropertiesModel& gasPropertiesModel_;

  //- Gas density [kg/m3]
  volScalarField& rho_g_;

  //- Gas volume fraction [-]
  volScalarField& eps_g_;

  //- Elemental diffusion coefficient [m2/s]
  PtrList<volScalarField>& Dm_;

  //- Mass model
  simpleMassModel& massModel_;

  //- Pyrolysis gas mass flow rate [kg/m3/s]
  const volVectorField& mDotG_;

  //- Pyrolysis gas mass flux flow rate interpolated on the faces [kg/m3/s] - updated in Mass
  const surfaceVectorField& mDotGFace_;

  //- Pyrolysis gas mass flow rate at the walls [kg/m2/s] - updated in BprimeBoundaryConditions
  const volScalarField& mDotGw_;

  //- Pyrolysis model
  simplePyrolysisModel& pyrolysisModel_;

  //- Pyrolysis elemental reaction rates [kg/m3/s]
  const PtrList<volScalarField>& pi_;

  //- FV options
  fv::options fvOptions_;

 public:

  //- Runtime type information
  TypeName("EquilibriumElement");

  // Constructors

  //- Construct from mesh and phase name
  EquilibriumElementMaterialChemistryModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~EquilibriumElementMaterialChemistryModel();

  // Member Functions

  //- Solve the element conservation equations
  virtual void update();

  //- Initialize the mixture
  inline Switch init();

  //- Update the BC for the element conservation equations
  inline void updateBC();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

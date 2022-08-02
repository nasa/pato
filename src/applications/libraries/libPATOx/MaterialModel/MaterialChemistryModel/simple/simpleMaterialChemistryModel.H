/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::simpleMaterialChemistryModel

Description

SourceFiles
    simpleMaterialChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleMaterialChemistryModel_H
#define simpleMaterialChemistryModel_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#undef Log // conflict between OpenFoam and Mutation++ on the alias "Log" -> undefined here from OpenFoam.
#include <mutation++/mutation++.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleMaterialChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class simpleMaterialChemistryModel
:
public IOdictionary
{

 private:

  // Private Member Functions

 protected:

  // Protected data

  //- Mesh database
  const fvMesh& mesh_;

  //- Material dictionary
  IOdictionary materialDict_;

  //- Debug Switch
  const Switch debug_;

  //- dynamic mesh flag
  const bool dynamicMesh_;

  //- Element names
  wordList elementNames_;

  //- Species names
  wordList speciesNames_;

  //- Mass fractions of the species/elements
  PtrList<volScalarField> massFractions_;

  //- Mass fractions of the species/elements from previous time step
  PtrList<volScalarField> oldMassFractions_;

  //- Mole fractions of the species/elements
  PtrList<volScalarField> moleFractions_;

  //- Mixture from Mutation++
  autoPtr<Mutation::Mixture> mixture_;

  //- Index species in Mutation++
  labelList speciesIndexMutation_;

  //- Latest estimation of MaterialChemistry integration
  scalar dtChem_;

 public:

  //- Runtime type information
  TypeName("simpleMaterialChemistryModel");


  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleMaterialChemistryModel,
      fvMesh,
      (const fvMesh& mesh,const word& dictName),
      (mesh, dictName)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleMaterialChemistryModel
  (
      const fvMesh& mesh,
      const word& dictName
  );


  //- Selector
  static autoPtr<simpleMaterialChemistryModel> New
  (
      const fvMesh& mesh,
      const word& dictName
  );

  //- Destructor
  virtual ~simpleMaterialChemistryModel();

  //- Update
  virtual void update() = 0;

  //- Return the elements names
  inline wordList& elementNames() const;

  //- Return the species names
  inline wordList& speciesNames() const;

  //- Return the elemental mass fractions
  inline PtrList<volScalarField>& massFractions() const;

  //- Return the elemental mass fractions of the previous time step
  inline PtrList<volScalarField>& oldMassFractions() const;

  //- Return the elemental mole fractions
  inline PtrList<volScalarField>& moleFractions() const;

  //- Return the Mutation++ mixture
  inline autoPtr<Mutation::Mixture>& mixture();

  //- Return the species index from Mutation++ compared to OpenFOAM
  inline labelList& speciesIndexMutation() const;

  //- Return the delta time for the chemistry model
  inline scalar& dtChem() const;

};

#include "simpleMaterialChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

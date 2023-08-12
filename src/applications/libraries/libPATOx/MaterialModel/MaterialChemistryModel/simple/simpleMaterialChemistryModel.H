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

#include "simpleModel.H"
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
public simpleModel
{

 private:

  // Private Member Functions

 protected:

  // Protected data

  //- Element names
  wordList elementNames_;

  //- Species names
  wordList speciesNames_;

  //- Mass fractions of the species/elements
  PList<volScalarField> massFractions_;

  //- Mass fractions of the species/elements from previous time step
  PList<volScalarField> oldMassFractions_;

  //- Mole fractions of the species/elements
  PList<volScalarField> moleFractions_;

  //- Mixture from Mutation++
  autoPtr<Mutation::Mixture> mixture_;

  //- Index species in Mutation++
  labelList speciesIndexMutation_;

  //- Latest estimation of MaterialChemistry integration
  scalar dtChem_;

 public:

  //- Runtime type information
  TypeName("simpleMaterialChemistryModel");

  //- Model name
  static const word modelName;

  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleMaterialChemistryModel,
      fvMesh,
      (const fvMesh& mesh,const word& regionName),
      (mesh, regionName)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleMaterialChemistryModel
  (
      const fvMesh& mesh,
      const word& regionName
  );


  //- Selector
  static autoPtr<simpleMaterialChemistryModel> New
  (
      const fvMesh& mesh,
      const word& regionName
  );

  //- Destructor
  virtual ~simpleMaterialChemistryModel();

  //- Update
  virtual void update() = 0;

  //- Return the elements names
  inline const wordList& elementNames() const;

  //- Return the species names
  inline const wordList& speciesNames() const;

  //- Return the elemental mass fractions
  inline const PList<volScalarField>& massFractions() const;

  //- Return the elemental mass fractions of the previous time step
  inline const PList<volScalarField>& oldMassFractions() const;

  //- Return the elemental mole fractions
  inline const PList<volScalarField>& moleFractions() const;

  //- Return the Mutation++ mixture
  inline const autoPtr<Mutation::Mixture>& mixture() const;

  //- Return the species index from Mutation++ compared to OpenFOAM
  inline const labelList& speciesIndexMutation() const;

  //- Return the delta time for the chemistry model
  inline const scalar& dtChem() const;

};

#include "simpleMaterialChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

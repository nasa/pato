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
    along with OpenFOAM.  If ConstantFiniteRatet, see <http://www.gnu.org/licenses/>.

Class
    Foam::ConstantFiniteRateMaterialChemistryModel

Description

SourceFiles
    ConstantFiniteRateMaterialChemistryModel.H
    ConstantFiniteRateMaterialChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef ConstantFiniteRateMaterialChemistryModel_H
#define ConstantFiniteRateMaterialChemistryModel_H

#include "simpleMaterialChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ConstantFiniteRateMaterialChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class ConstantFiniteRateMaterialChemistryModel
:
public simpleMaterialChemistryModel
{

 protected:

  //- Mutation++ mixture name
  word mixtureMutationFiniteRate;

  //- Mutation++ mixture
  autoPtr<Mutation::Mixture>& mixFiniteRate;

  //- Mutation++ options
  autoPtr<Mutation::MixtureOptions> optsFiniteRate;

  //- Number of species
  int ns_mix;

  //- Elements names
  wordList& elementNames_;

  //- Species names
  wordList& speciesNames_;

  //- Species mass fractions [-]
  PtrList<volScalarField>& massFractions_;

  //- Species mole fractions [-]
  PtrList<volScalarField>& moleFractions_;

 public:

  //- Runtime type information
  TypeName("ConstantFiniteRate");

  // Constructors

  //- Construct from mesh and phase name
  ConstantFiniteRateMaterialChemistryModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~ConstantFiniteRateMaterialChemistryModel();

  //- Member Functions

  //- Update
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

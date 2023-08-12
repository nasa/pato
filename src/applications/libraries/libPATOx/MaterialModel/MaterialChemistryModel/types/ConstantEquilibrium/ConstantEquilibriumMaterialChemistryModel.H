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
    along with OpenFOAM.  If ConstantEquilibriumt, see <http://www.gnu.org/licenses/>.

Class
    Foam::ConstantEquilibriumMaterialChemistryModel

Description

SourceFiles
    ConstantEquilibriumMaterialChemistryModel.H
    ConstantEquilibriumMaterialChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef ConstantEquilibriumMaterialChemistryModel_H
#define ConstantEquilibriumMaterialChemistryModel_H

#include "simpleMaterialChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ConstantEquilibriumMaterialChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class ConstantEquilibriumMaterialChemistryModel
:
public simpleMaterialChemistryModel
{

 protected:

  //- Mutation++ mixture name
  word mixtureMutation;

  //- Mutation++ mixture options
  autoPtr<Mutation::MixtureOptions> optEq_;

  //- Mutation++ mixture
  autoPtr<Mutation::Mixture>& mix_;

  //- Number of elements in the mixture
  int ne_mix;

  //- Mass fractions of elements
  scalarList p_Z;

  //- Reference of mass fractions of elements
  scalarList p_Z0;

  //- Mass fractions of elements in the subsurface cell
  scalarList p_Zsub;

  //- Mole fractions of elements
  scalarList p_Zx;

  //- Elements names
  wordList& elementNames_;

  //- Species names
  wordList& speciesNames_;

  //- Elemental mass fractions
  PList<volScalarField>& massFractions_;

  //- Elemental mole fractions
  PList<volScalarField>& moleFractions_;

  //- Gaseous mole fraction in the material
  List<Tuple2<word,scalar> > moleFractionGasInMaterial_;

 public:

  //- Runtime type information
  TypeName("ConstantEquilibrium");

  // Constructors

  //- Construct from mesh and phase name
  ConstantEquilibriumMaterialChemistryModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~ConstantEquilibriumMaterialChemistryModel();

  //- Member Functions

  //- Update
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

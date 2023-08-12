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

\*---------------------------------------------------------------------------*/

#include "simpleMaterialChemistryModel.H"


/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simpleMaterialChemistryModel, 0);
defineRunTimeSelectionTable(simpleMaterialChemistryModel, fvMesh);
}

/* * * * * * * * * * * * * * * public static members * * * * * * * * * * * * * */

const Foam::word Foam::simpleMaterialChemistryModel::modelName="MaterialChemistry";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleMaterialChemistryModel::simpleMaterialChemistryModel
(
    const fvMesh& mesh,
    const word& regionName
):
simpleModel(mesh,regionName,modelName)
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simpleMaterialChemistryModel> Foam::simpleMaterialChemistryModel::New
(
    const fvMesh& mesh,
    const word& regionName
)
{
  IOdictionary MaterialDict
  (
      IOobject
      (
          regionName+"Properties",
          mesh.time().constant(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE,
          false
      )
  );

  word MaterialChemistryModelTypeName= "no";

  if (MaterialDict.isDict("MaterialChemistry")) {
    dictionary MaterialChemistryModelDict = MaterialDict.subDict("MaterialChemistry");

    MaterialChemistryModelTypeName =
        word(MaterialChemistryModelDict.lookupOrDefault<word>("MaterialChemistryType","no"));

  }

  typename simpleMaterialChemistryModel::fvMeshConstructorTable::iterator cstrIter =
      simpleMaterialChemistryModel::fvMeshConstructorTablePtr_->find(MaterialChemistryModelTypeName);

  if (cstrIter == simpleMaterialChemistryModel::fvMeshConstructorTablePtr_->end()) {

    wordList list_available =  simpleMaterialChemistryModel::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simpleMaterialChemistryModel::typeName << " type "
        << MaterialChemistryModelTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<simpleMaterialChemistryModel>(cstrIter()(mesh, regionName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleMaterialChemistryModel::~simpleMaterialChemistryModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //


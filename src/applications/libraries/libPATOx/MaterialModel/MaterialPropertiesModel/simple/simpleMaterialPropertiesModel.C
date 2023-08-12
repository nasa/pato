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

#include "simpleMaterialPropertiesModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simpleMaterialPropertiesModel, 0);
defineRunTimeSelectionTable(simpleMaterialPropertiesModel, fvMesh);
}

/* * * * * * * * * * * * * * * public static members * * * * * * * * * * * * * */

const Foam::word Foam::simpleMaterialPropertiesModel::modelName="MaterialProperties";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleMaterialPropertiesModel::simpleMaterialPropertiesModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleModel(mesh,regionName,modelName),
I_(1,0,0,0,1,0,0,0,1)
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simpleMaterialPropertiesModel> Foam::simpleMaterialPropertiesModel::New
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


  word MaterialPropertiesModelTypeName= "no";

  if (MaterialDict.isDict("MaterialProperties")) {
    dictionary MaterialPropertiesModelDict =
        MaterialDict.subDict("MaterialProperties");

    MaterialPropertiesModelTypeName =
        word(MaterialPropertiesModelDict.lookupOrDefault<word>("MaterialPropertiesType","no"));
  }

  typename simpleMaterialPropertiesModel::fvMeshConstructorTable::iterator
  cstrIter =
      simpleMaterialPropertiesModel::fvMeshConstructorTablePtr_->find
      (
          MaterialPropertiesModelTypeName
      );

  if
  (
      cstrIter == simpleMaterialPropertiesModel::fvMeshConstructorTablePtr_->end()
  ) {
    wordList list_available =
        simpleMaterialPropertiesModel::fvMeshConstructorTablePtr_
        ->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simpleMaterialPropertiesModel::typeName << " type "
        << MaterialPropertiesModelTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<simpleMaterialPropertiesModel>(cstrIter()(mesh, regionName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleMaterialPropertiesModel::~simpleMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //

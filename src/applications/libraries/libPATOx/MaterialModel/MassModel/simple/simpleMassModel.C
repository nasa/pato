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

#include "simpleMassModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simpleMassModel, 0);
defineRunTimeSelectionTable(simpleMassModel, fvMesh);
}

/* * * * * * * * * * * * * * * public static members * * * * * * * * * * * * * */

const Foam::word Foam::simpleMassModel::modelName="Mass";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleMassModel::simpleMassModel
(
    const fvMesh& mesh,
    const word& regionName
):
simpleModel(mesh,regionName,modelName)
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simpleMassModel> Foam::simpleMassModel::New
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

  word MassModelTypeName= "no";

  if (MaterialDict.isDict("Mass")) {
    dictionary MassModelDict = MaterialDict.subDict("Mass");

    MassModelTypeName =
        word(MassModelDict.lookupOrDefault<word>("MassType","no"));

  }

  typename simpleMassModel::fvMeshConstructorTable::iterator cstrIter =
      simpleMassModel::fvMeshConstructorTablePtr_->find(MassModelTypeName);

  if (cstrIter == simpleMassModel::fvMeshConstructorTablePtr_->end()) {

    wordList list_available =  simpleMassModel::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simpleMassModel::typeName << " type "
        << MassModelTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<simpleMassModel>(cstrIter()(mesh, regionName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleMassModel::~simpleMassModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //


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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleMassModel::simpleMassModel
(
    const fvMesh& mesh,
    const word& dictName
):
IOdictionary
(
    IOobject
    (
        dictName+"MassModel",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
mesh_(mesh),
materialDict_(
    IOobject
    (
        dictName+"Properties",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    )
),
debug_(materialDict_.lookupOrDefault<Switch>("debug", false)),
dynamicMesh_(isA<dynamicFvMesh>(mesh)),
mDotG_(meshLookupOrConstructVector(mesh,"mDotG", dimensionedVector("zero", dimMass/(dimLength*dimLength*dimTime), vector(0, 0, 0)))),
vG_(meshLookupOrConstructVector(mesh,"vG",dimensionedVector("zero", dimLength/dimTime, vector(0, 0, 0)))),
mDotGFace_(meshLookupOrConstructSurfaceVector(mesh, "mDotGface", linearInterpolate(mDotG_)))
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simpleMassModel> Foam::simpleMassModel::New
(
    const fvMesh& mesh,
    const word& dictName
)
{
  IOdictionary MaterialDict
  (
      IOobject
      (
          dictName+"Properties",
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
        word(MassModelDict.lookup("MassType"));

  }

  Info<< "Selecting MassModel type \"" << MassModelTypeName << "\"" << endl;

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

  return autoPtr<simpleMassModel>(cstrIter()(mesh, dictName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleMassModel::~simpleMassModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //


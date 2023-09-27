/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "simpleOptionModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simpleOptionModel, 0);
defineRunTimeSelectionTable(simpleOptionModel, fvMesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleOptionModel::simpleOptionModel(const fvMesh& mesh, const dictionary& dict, const label& patchID, const word& optionTypeName)
  :
IOdictionary
(
    IOobject
    (
        "simpleOptionModel",
        mesh.time().constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
),
mesh_(mesh),
dict_(dict),
patchID_(patchID),
optionTypeName_(optionTypeName)
{
}

// ************************************************************************* //

Foam::autoPtr<Foam::simpleOptionModel> Foam::simpleOptionModel::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const label& patchID,
    const word& optionTypeName
)
{
  Info<< "Selecting OptionModel type \"" << optionTypeName << "\"" << endl;

  typename simpleOptionModel::fvMeshConstructorTable::iterator cstrIter =
      simpleOptionModel::fvMeshConstructorTablePtr_->find(optionTypeName);

  if (cstrIter == simpleOptionModel::fvMeshConstructorTablePtr_->end()) {

    wordList list_available =  simpleOptionModel::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simpleOptionModel::typeName << " type "
        << optionTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<simpleOptionModel>(cstrIter()(mesh, dict, patchID, optionTypeName));
}

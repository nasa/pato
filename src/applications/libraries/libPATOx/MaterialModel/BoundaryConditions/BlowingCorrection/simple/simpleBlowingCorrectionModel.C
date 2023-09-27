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

#include "simpleBlowingCorrectionModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simpleBlowingCorrectionModel, 0);
defineRunTimeSelectionTable(simpleBlowingCorrectionModel, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleBlowingCorrectionModel::simpleBlowingCorrectionModel
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& regionName,
    const label& currentPatchID,
    const dictionary dict
):
IOdictionary
(
    IOobject
    (
        "BlowingCorrectionModel",
        mesh.time().constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
),
mesh_(mesh),
phaseName_(phaseName),
regionName_(regionName),
currentPatchID_(currentPatchID),
dict_(dict),
initialized_(false),
blowingCorrectionType_(dict.lookupOrDefault<word>("blowingCorrectionType","constantLambda")),
debug_(dict_.lookupOrDefault<Switch>("debug","no"))
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simpleBlowingCorrectionModel> Foam::simpleBlowingCorrectionModel::New
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& regionName,
    const label& currentPatchID,
    const dictionary dict
)
{

  word BlowingCorrectionModelTypeName =
      word(dict.lookupOrDefault<word>("blowingCorrectionType","constantLambda"));

  Info<< simpleModel::getTabLevel("| ") << "Selecting BlowingCorrectionModel type \"" << BlowingCorrectionModelTypeName << "\"" << endl;

  typename simpleBlowingCorrectionModel::fvMeshConstructorTable::iterator cstrIter =
      simpleBlowingCorrectionModel::fvMeshConstructorTablePtr_->find(BlowingCorrectionModelTypeName);

  if (cstrIter == simpleBlowingCorrectionModel::fvMeshConstructorTablePtr_->end()) {

    wordList list_available =  simpleBlowingCorrectionModel::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simpleBlowingCorrectionModel::typeName << " type "
        << BlowingCorrectionModelTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<simpleBlowingCorrectionModel>(cstrIter()(mesh, phaseName, regionName, currentPatchID, dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleBlowingCorrectionModel::~simpleBlowingCorrectionModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //


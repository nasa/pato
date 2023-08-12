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

#include "simpleBoundaryMappingModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simpleBoundaryMappingModel, 0);
defineRunTimeSelectionTable(simpleBoundaryMappingModel, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleBoundaryMappingModel::simpleBoundaryMappingModel
(
    const fvMesh& mesh,
    const wordList& neededFields,
    const dictionary dict
):
IOdictionary
(
    IOobject
    (
        "BoundaryMappingModel",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
mesh_(mesh),
neededFields_(neededFields),
dict_(dict),
debug_(dict_.lookupOrDefault<Switch>("debug","no")),
mappingFileName_(dict_.lookup("mappingFileName")),
mappingFields_(dict_.lookup("mappingFields")),
currentTimePatchesDataFields_(mesh.boundaryMesh().size())
{
  foundFieldsInMesh(mesh_,neededFields_);

  forAll(mappingFields_, fieldI) { // e.g. mappingFields ((p  3) (k 2));
    mappingFieldsName_.append(mappingFields_[fieldI].first()); // mappingFieldsName_ = [p k]
    word num_ = mappingFields_[fieldI].second();
    if (!isNumber(num_)) {
      FatalErrorInFunction << num_ << " from mappingFields_ is not a number." << exit(FatalError);
    }
    mappingFieldsColumn_.append(stof(num_));  // mappingFieldsColumn_ = [3 2]
  }

  foundFieldsInMesh(mesh_,mappingFieldsName_);
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simpleBoundaryMappingModel> Foam::simpleBoundaryMappingModel::New
(
    const fvMesh& mesh,
    const wordList& neededFields,
    const dictionary dict
)
{

  word BoundaryMappingModelTypeName =
      word(dict.lookup("mappingType"));

  Info<< simpleModel::getTabLevel() << "Selecting BoundaryMappingModel type \"" << BoundaryMappingModelTypeName << "\"" << endl;

  typename simpleBoundaryMappingModel::fvMeshConstructorTable::iterator cstrIter =
      simpleBoundaryMappingModel::fvMeshConstructorTablePtr_->find(BoundaryMappingModelTypeName);

  if (cstrIter == simpleBoundaryMappingModel::fvMeshConstructorTablePtr_->end()) {

    wordList list_available =  simpleBoundaryMappingModel::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simpleBoundaryMappingModel::typeName << " type "
        << BoundaryMappingModelTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<simpleBoundaryMappingModel>(cstrIter()(mesh, neededFields, dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleBoundaryMappingModel::~simpleBoundaryMappingModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //


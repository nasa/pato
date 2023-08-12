/*---------------------------------------------------------------------------* \
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

#include "simpleMaterialFailureCriteriaModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simpleMaterialFailureCriteriaModel, 0);
defineRunTimeSelectionTable(simpleMaterialFailureCriteriaModel, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleMaterialFailureCriteriaModel::simpleMaterialFailureCriteriaModel
(
    const fvMesh& mesh,
    const word& dictName
)  :
IOdictionary
(
    IOobject
    (
        dictName+"MaterialFailureCriteriaModel",
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
dynamicMesh_(isA<dynamicFvMesh>(mesh))
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simpleMaterialFailureCriteriaModel> Foam::simpleMaterialFailureCriteriaModel::New
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

  word MaterialFailureModelTypeName= "no";

  if (MaterialDict.isDict("SolidMechanics")) {
    dictionary SolidMechanicsModelDict = MaterialDict.subDict("SolidMechanics");
    if (SolidMechanicsModelDict.isDict("MaterialFailure")) {
      dictionary FailureModelDict = SolidMechanicsModelDict.subDict("MaterialFailure");
      MaterialFailureModelTypeName = word(FailureModelDict.lookup("MaterialFailureCriteriaType"));
    }
  }

  Info<< "Selecting MaterialFailureCriteriaModel type \"" << MaterialFailureModelTypeName << "\"" << endl;

  typename simpleMaterialFailureCriteriaModel::fvMeshConstructorTable::iterator cstrIter =
      simpleMaterialFailureCriteriaModel::fvMeshConstructorTablePtr_->find(MaterialFailureModelTypeName);

  if (cstrIter == simpleMaterialFailureCriteriaModel::fvMeshConstructorTablePtr_->end()) {

    wordList list_available =  simpleMaterialFailureCriteriaModel::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simpleMaterialFailureCriteriaModel::typeName << " type "
        << MaterialFailureModelTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }
  return autoPtr<simpleMaterialFailureCriteriaModel>(cstrIter()(mesh, dictName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleMaterialFailureCriteriaModel::~simpleMaterialFailureCriteriaModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //


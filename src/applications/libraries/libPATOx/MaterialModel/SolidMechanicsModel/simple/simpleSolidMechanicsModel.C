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

#include "simpleSolidMechanicsModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simpleSolidMechanicsModel, 0);
defineRunTimeSelectionTable(simpleSolidMechanicsModel, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleSolidMechanicsModel::simpleSolidMechanicsModel
(
    const fvMesh& mesh,
    const word& dictName
):
IOdictionary
(
    IOobject
    (
        dictName+"SolidMechanicsModel",
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
nu_(meshLookupOrConstructScalar(mesh,"nu",dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.3))),
E_(meshLookupOrConstructScalar(mesh,"E",dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0, 0, 0), 24e06))),// 12e9/500 Keeping OpenFoam notation for now E = E/rho (needs to be changed at some point)
alpha_(meshLookupOrConstructScalar(mesh,"alpha",dimensionedScalar("zero", dimensionSet(0, 0, 0, -1, 0, 0, 0), 5.4e-05))),
rho_(meshLookupOrConstructScalar(mesh,"rho",dimensionedScalar("zero", dimensionSet(1, -3, 0, 0, 0, 0, 0), 500)))
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simpleSolidMechanicsModel> Foam::simpleSolidMechanicsModel::New
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


  word SolidMechanicsModelTypeName= "no";

  if (MaterialDict.isDict("SolidMechanics")) {
    dictionary SolidMechanicsModelDict = MaterialDict.subDict("SolidMechanics");

    SolidMechanicsModelTypeName =
        word(SolidMechanicsModelDict.lookup("SolidMechanicsType"));

  }

  Info<< "Selecting SolidMechanicsModel type \"" << SolidMechanicsModelTypeName << "\"" << endl;

  typename simpleSolidMechanicsModel::fvMeshConstructorTable::iterator cstrIter =
      simpleSolidMechanicsModel::fvMeshConstructorTablePtr_->find(SolidMechanicsModelTypeName);

  if (cstrIter == simpleSolidMechanicsModel::fvMeshConstructorTablePtr_->end()) {

    wordList list_available =  simpleSolidMechanicsModel::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simpleSolidMechanicsModel::typeName << " type "
        << SolidMechanicsModelTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<simpleSolidMechanicsModel>(cstrIter()(mesh, dictName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleSolidMechanicsModel::~simpleSolidMechanicsModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //


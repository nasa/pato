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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleMaterialPropertiesModel::simpleMaterialPropertiesModel
(
    const fvMesh& mesh,
    const word& dictName
):
IOdictionary
(
    IOobject
    (
        dictName+"MaterialPropertiesModel",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
mesh_(mesh),
dictName_(dictName),
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
cp_(meshLookupOrConstructScalar(mesh,"cp", dimensionedScalar("0", dimensionSet(0,2,-2,-1,0,0,0), scalar(0)))),
k_(meshLookupOrConstructTensor(mesh,"k",  dimensionedTensor("0", dimensionSet(1, 1, -3, -1, 0, 0, 0), tensor(1,0,0,0,1,0,0,0,1)))),
K_(meshLookupOrConstructTensor(mesh,"K",  dimensionedTensor("0", dimLength*dimLength, tensor(0,0,0,0,0,0,0,0,0)) )),
rho_s_(meshLookupOrConstructScalar(mesh, "rho_s")),
h_bar_(meshLookupOrConstructScalar(mesh,"h_bar", dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0)))),
emissivity_(meshLookupOrConstructScalar(mesh,"emissivity", dimensionedScalar("0", dimless, scalar(0)))),
absorptivity_(meshLookupOrConstructScalar(mesh,"absorptivity", dimensionedScalar("0", dimless, scalar(0)))),
h_c_(meshLookupOrConstructScalar(mesh,"h_c", dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0)))),
pyrolysisFlux_(meshLookupOrConstructScalar(mesh, "pyrolysisFlux", dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), scalar(0)))),
nSolidPhases_(0)
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simpleMaterialPropertiesModel> Foam::simpleMaterialPropertiesModel::New
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


  word MaterialPropertiesModelTypeName= "no";

  if (MaterialDict.isDict("MaterialProperties")) {
    dictionary MaterialPropertiesModelDict = MaterialDict.subDict("MaterialProperties");

    MaterialPropertiesModelTypeName =
        word(MaterialPropertiesModelDict.lookup("MaterialPropertiesType"));

  }

  Info<< "Selecting MaterialPropertiesModel type \"" << MaterialPropertiesModelTypeName << "\"" << endl;

  typename simpleMaterialPropertiesModel::fvMeshConstructorTable::iterator cstrIter =
      simpleMaterialPropertiesModel::fvMeshConstructorTablePtr_->find(MaterialPropertiesModelTypeName);

  if (cstrIter == simpleMaterialPropertiesModel::fvMeshConstructorTablePtr_->end()) {

    wordList list_available =  simpleMaterialPropertiesModel::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simpleMaterialPropertiesModel::typeName << " type "
        << MaterialPropertiesModelTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<simpleMaterialPropertiesModel>(cstrIter()(mesh, dictName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleMaterialPropertiesModel::~simpleMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //


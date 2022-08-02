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

#include "simpleGasPropertiesModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simpleGasPropertiesModel, 0);
defineRunTimeSelectionTable(simpleGasPropertiesModel, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleGasPropertiesModel::simpleGasPropertiesModel
(
    const fvMesh& mesh,
    const word& dictName
):
IOdictionary
(
    IOobject
    (
        dictName+"GasPropertiesModel",
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
M_(meshLookupOrConstructScalar(mesh,"M_g",dimensionedScalar("zero", dimensionSet(1, 0, 0, 0, -1, 0, 0), 0.02))),
cp_g_(meshLookupOrConstructScalar(mesh,"cp_g",dimensionedScalar("zero", dimensionSet(0, 2, -2, -1, 0, 0, 0), 0.0) )),
mu_(meshLookupOrConstructScalar(mesh,"mu_g",dimensionedScalar("mu", dimensionSet(1, -1, -1, 0, 0, 0, 0), 2e-5))),
h_g_(meshLookupOrConstructScalar(mesh,"h_g",dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0) )),
rho_g_(meshLookupOrConstructScalar(mesh,"rho_g",dimensionedScalar("zero", dimensionSet(1, -3, 0, 0, 0, 0, 0), 0.0) )),
k_g_(meshLookupOrConstructScalar(mesh,"k_g",dimensionedScalar("zero", dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0) )),
eps_g_(meshLookupOrConstructScalar(mesh,"eps_g",dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0))),
Ediff_(meshLookupOrConstructScalar(mesh,"Ediff",dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)))
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simpleGasPropertiesModel> Foam::simpleGasPropertiesModel::New
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


  word GasPropertiesModelTypeName= "no";

  if (MaterialDict.isDict("GasProperties")) {
    dictionary GasPropertiesModelDict = MaterialDict.subDict("GasProperties");

    GasPropertiesModelTypeName =
        word(GasPropertiesModelDict.lookup("GasPropertiesType"));

  }

  Info<< "Selecting GasPropertiesModel type \"" << GasPropertiesModelTypeName << "\"" << endl;

  typename simpleGasPropertiesModel::fvMeshConstructorTable::iterator cstrIter =
      simpleGasPropertiesModel::fvMeshConstructorTablePtr_->find(GasPropertiesModelTypeName);

  if (cstrIter == simpleGasPropertiesModel::fvMeshConstructorTablePtr_->end()) {

    wordList list_available =  simpleGasPropertiesModel::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simpleGasPropertiesModel::typeName << " type "
        << GasPropertiesModelTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<simpleGasPropertiesModel>(cstrIter()(mesh, dictName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleGasPropertiesModel::~simpleGasPropertiesModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //


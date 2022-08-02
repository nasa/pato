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

#include "simplePyrolysisModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simplePyrolysisModel, 0);
defineRunTimeSelectionTable(simplePyrolysisModel, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simplePyrolysisModel::simplePyrolysisModel
(
    const fvMesh& mesh,
    const word& dictName
):
IOdictionary
(
    IOobject
    (
        dictName+"PyrolysisModel",
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
piTotal_(meshLookupOrConstructScalar(mesh, "piTotal",  dimensionedScalar("piTotal0", dimensionSet(1,-3,-1,0,0,0,0), 0))),
tau_(meshLookupOrConstructScalar(mesh, "tau",  dimensionedScalar("tau", dimless, 1.0), "zeroGradient")),
Ta_(meshLookupOrConstructScalar(mesh, "Ta")),
rho_s_(meshLookupOrConstructScalar(mesh, "rho_s")),
rho_v_(meshLookupOrConstructScalar(mesh, "rho_v",  dimensionedScalar("rho_v0", dimMass/dimVolume, 0))),
rho_c_(meshLookupOrConstructScalar(mesh, "rho_c",  dimensionedScalar("rho_c0", dimMass/dimVolume, 0)))
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simplePyrolysisModel> Foam::simplePyrolysisModel::New
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


  word PyrolysisModelTypeName= "no";

  if (MaterialDict.isDict("Pyrolysis")) {
    dictionary PyrolysisModelDict = MaterialDict.subDict("Pyrolysis");

    PyrolysisModelTypeName =
        word(PyrolysisModelDict.lookup("PyrolysisType"));

  }

  Info<< "Selecting PyrolysisModel type \"" << PyrolysisModelTypeName << "\"" << endl;

  typename simplePyrolysisModel::fvMeshConstructorTable::iterator cstrIter =
      simplePyrolysisModel::fvMeshConstructorTablePtr_->find(PyrolysisModelTypeName);

  if (cstrIter == simplePyrolysisModel::fvMeshConstructorTablePtr_->end()) {

    wordList list_available =  simplePyrolysisModel::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simplePyrolysisModel::typeName << " type "
        << PyrolysisModelTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<simplePyrolysisModel>(cstrIter()(mesh, dictName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simplePyrolysisModel::~simplePyrolysisModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //


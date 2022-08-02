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

#include "simpleIOModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(simpleIOModel, 0);
defineRunTimeSelectionTable(simpleIOModel, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleIOModel::simpleIOModel
(
    const fvMesh& mesh,
    const word& dictName
):
IOdictionary
(
    IOobject
    (
        dictName+"IOModel",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
mesh_(mesh),
phaseName_(""),
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
infoDebug_(materialDict_.isDict("IO")?materialDict_.subDict("IO").lookupOrDefault<label>("infoDebug", lduMatrix::debug): lduMatrix::debug),
materialDictPath_(mesh.time().constant()+"/"+dictName+"/" + IOobject::groupName(dictName+"Properties", phaseName_)),
readFilesList_(materialDict_.isDict("IO")?materialDict_.subDict("IO").lookupOrDefault<List<fileName> >("readFiles",nullListFileName):nullListFileName)
{
  if (infoDebug_ != lduMatrix::debug) {
    // IO of the solvers
    lduMatrix::debug=infoDebug_;
    solverPerformance::debug=infoDebug_;
    Info << "lduMatrix::debug=" << infoDebug_ << endl;
    Info << "solverPerformance::debug=" << infoDebug_ << endl;
  }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simpleIOModel> Foam::simpleIOModel::New
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


  word IOModelTypeName= "no";

  if (MaterialDict.isDict("IO")) {
    dictionary IOModelDict = MaterialDict.subDict("IO");

    IOModelTypeName =
        word(IOModelDict.lookupOrDefault<word>("IOType","no"));

  }

  Info<< "Selecting IOModel type \"" << IOModelTypeName << "\"" << endl;

  typename simpleIOModel::fvMeshConstructorTable::iterator cstrIter =
      simpleIOModel::fvMeshConstructorTablePtr_->find(IOModelTypeName);

  if (cstrIter == simpleIOModel::fvMeshConstructorTablePtr_->end()) {

    wordList list_available =  simpleIOModel::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simpleIOModel::typeName << " type "
        << IOModelTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<simpleIOModel>(cstrIter()(mesh, dictName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleIOModel::~simpleIOModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //


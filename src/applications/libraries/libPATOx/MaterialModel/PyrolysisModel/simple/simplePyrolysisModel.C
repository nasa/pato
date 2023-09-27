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

/* * * * * * * * * * * * * * * public static members * * * * * * * * * * * * * */

const Foam::word Foam::simplePyrolysisModel::modelName="Pyrolysis";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simplePyrolysisModel::simplePyrolysisModel
(
    const fvMesh& mesh,
    const word& regionName
):
simpleModel(mesh,regionName,modelName)
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::simplePyrolysisModel> Foam::simplePyrolysisModel::New
(
    const fvMesh& mesh,
    const word& regionName
)
{
  IOdictionary MaterialDict
  (
      IOobject
      (
          regionName+"Properties",
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
        word(PyrolysisModelDict.lookupOrDefault<word>("PyrolysisType","no"));

  }

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

  return autoPtr<simplePyrolysisModel>(cstrIter()(mesh, regionName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simplePyrolysisModel::~simplePyrolysisModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //


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

/* * * * * * * * * * * * * * * public static members * * * * * * * * * * * * * */

const Foam::word Foam::simpleSolidMechanicsModel::modelName="SolidMechanics";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleSolidMechanicsModel::simpleSolidMechanicsModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleModel(mesh,regionName,modelName)
{
  if (materialDict_.isDict("SolidMechanics")) {
    planeStress_ =
        materialDict_.subDict("SolidMechanics").lookupOrDefault<Switch>
        (
            "planeStress",
            false
        );
  }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::simpleSolidMechanicsModel> Foam::simpleSolidMechanicsModel::New
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


  word SolidMechanicsModelTypeName= "no";

  if (MaterialDict.isDict("SolidMechanics")) {
    dictionary SolidMechanicsModelDict =
        MaterialDict.subDict("SolidMechanics");

    SolidMechanicsModelTypeName =
        word(SolidMechanicsModelDict.lookupOrDefault<word>("SolidMechanicsType","no"));
  }

  typename simpleSolidMechanicsModel::fvMeshConstructorTable::iterator
  cstrIter =
      simpleSolidMechanicsModel::fvMeshConstructorTablePtr_->find
      (
          SolidMechanicsModelTypeName
      );

  if
  (
      cstrIter == simpleSolidMechanicsModel::fvMeshConstructorTablePtr_->end()
  ) {
    wordList list_available =
        simpleSolidMechanicsModel::fvMeshConstructorTablePtr_->sortedToc();

    FatalErrorInFunction
        << "Unknown " << simpleSolidMechanicsModel::typeName << " type "
        << SolidMechanicsModelTypeName << nl << nl
        << "Valid types are:"
        << list_available
        << exit(FatalError);
  }

  return autoPtr<simpleSolidMechanicsModel>(cstrIter()(mesh, regionName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleSolidMechanicsModel::~simpleSolidMechanicsModel()
{}


// * * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

// ************************************************************************* //


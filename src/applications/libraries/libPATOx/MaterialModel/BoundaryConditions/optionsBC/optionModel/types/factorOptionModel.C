/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "factorOptionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::factorOptionModel::factorOptionModel(const fvMesh& mesh, const dictionary& dict, const label& patchID, const word& optionTypeName)
  :
simpleOptionModel(mesh, dict, patchID, optionTypeName),
factorList(dict.lookup("factorList"))
{
  print(); // print the fields and factors from factorList
  // check if the field names exist in the mesh (volScalarField)
  forAll(factorList, i) {
    word name = factorList[i].first();
    if (!mesh.objectRegistry::template foundObject<volScalarField>(name)) {
      FatalError << "factorList::fieldName(" << name << ") is not found in mesh (volScalarField)." << exit(FatalError);
    }
  }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::factorOptionModel::update()
{
  if (mesh_.time().value()==currentTime_) {
    return;
  }
  currentTime_=mesh_.time().value();
  print(); // print the fields and factors from factorList
  // update the fields in factorList
  forAll(factorList, i) {
    volScalarField& field = const_cast<volScalarField&>(mesh_.objectRegistry::template lookupObject<volScalarField>(factorList[i].first()));
    scalarList& bf_field = field.boundaryFieldRef()[patchID_];
    forAll(bf_field, faceI) {
      bf_field[faceI] *= factorList[i].second();
    }
  }
}

void Foam::factorOptionModel::print()
{
  Info << "factorOption multiplies the values of the patch " << mesh_.boundaryMesh()[patchID_].name() << " for the fields [";
  forAll(factorList, i) {
    Info << factorList[i].first();
    if (i < factorList.size()-1) {
      Info << " ";
    }
  }
  Info << "] by the factors [";
  forAll(factorList, i) {
    Info << factorList[i].second();
    if (i < factorList.size()-1) {
      Info << " ";
    }
  }
  Info << "] respectively." << endl;
}

// ************************************************************************* //

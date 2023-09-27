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

#include "timeFactorOptionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeFactorOptionModel::timeFactorOptionModel(const fvMesh& mesh, const dictionary& dict, const label& patchID, const word& optionTypeName)
  :
simpleOptionModel(mesh, dict, patchID, optionTypeName),
timeFactorFile_(dict.lookup("timeFactorFile")),
timeFactorList_(dict.lookup("timeFactorList")),
timeFactorData_(readFileData(timeFactorFile_))
{
  print(); // print the fields and time factor indexes from timeFactorList
  // check if the field names exist in the mesh (volScalarField)
  forAll(timeFactorList_, i) {
    word name = timeFactorList_[i].first();
    if (!mesh.objectRegistry::template foundObject<volScalarField>(name)) {
      FatalError << "factorList::fieldName(" << name << ") is not found in mesh (volScalarField)." << exit(FatalError);
    }
  }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeFactorOptionModel::update(label faceI)
{
  print(); // print the fields and time factor indexes from timeFactorList
  // update the fields in timeFactorList
  forAll(timeFactorList_, i) {
    volScalarField& field = const_cast<volScalarField&>(mesh_.objectRegistry::template lookupObject<volScalarField>(timeFactorList_[i].first()));
    scalarList& bf_field = field.boundaryFieldRef()[patchID_];
    scalar factor = linearInterpolation(timeFactorData_[0],timeFactorData_[timeFactorList_[i].second()], mesh_.time().value());
    bf_field[faceI] += factor;
  }
}

void Foam::timeFactorOptionModel::print()
{
  Info << "timeFactorOption add the values of the patch " << mesh_.boundaryMesh()[patchID_].name() << " for the fields [";
  forAll(timeFactorList_, i) {
    Info << timeFactorList_[i].first();
    if (i < timeFactorList_.size()-1) {
      Info << " ";
    }
  }
  Info << "] using a temporal linear interpolation from the file " << timeFactorFile_;
  Info << " with the factor indexes [";
  forAll(timeFactorList_, i) {
    Info << timeFactorList_[i].second();
    if (i < timeFactorList_.size()-1) {
      Info << " ";
    }
  }
  Info << "] respectively." << endl;
}

// ************************************************************************* //

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

#include "simpleBoundaryMappingModel.H"

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline void Foam::simpleBoundaryMappingModel::init()
{
  // Take constant value from BoundaryConditions for the needed fields not found in mappingFields
  forAll(neededFields_, fieldI) {

    if(!foundInList(neededFields_[fieldI],mappingFieldsName_ )) {
      if (mesh_.objectRegistry::foundObject<volScalarField>(neededFields_[fieldI])) {
        if(debug_) {
          Info << "--- neededFields_[fieldI] = " << neededFields_[fieldI] << " --- Foam::BoundaryMapping::BoundaryMapping " << endl;
        }
        volScalarField& field_ = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(neededFields_[fieldI]));

        scalar fieldConstantValue_(readScalar(dict_.lookup(neededFields_[fieldI])));
        constantFields_.append(Tuple2<word,scalar>(neededFields_[fieldI],fieldConstantValue_));
        if(debug_) {
          Info << "--- fieldConstantValue_ = " << fieldConstantValue_ << " --- Foam::BoundaryMapping::BoundaryMapping " << endl;
        }
        forAll(field_, i) {
          field_[i]= fieldConstantValue_;
        }
        if(debug_) {
          Info << "--- boundaryField neededFields_[fieldI] = " << neededFields_[fieldI] << " --- Foam::BoundaryMapping::BoundaryMapping " << endl;
          Info << "--- fieldConstantValue_  = " << fieldConstantValue_ << " --- Foam::BoundaryMapping::BoundaryMapping " << endl;
          Info << "--- field_ = " << field_ <<  " --- Foam::BoundaryMapping::BoundaryMapping " << endl;
        }
        forAll(mesh_.boundaryMesh(),patchI) {
          forAll(mesh_.boundaryMesh()[patchI], faceI) {
            field_.boundaryFieldRef()[patchI][faceI] = fieldConstantValue_;
          }
        }
        if(debug_) {
          Info << "--- end neededFields_[fieldI] = " << neededFields_[fieldI] << " --- Foam::BoundaryMapping::BoundaryMapping " << endl;
        }
      }

      if (mesh_.objectRegistry::foundObject<volVectorField>(neededFields_[fieldI])) {

        volVectorField& field_ = const_cast<volVectorField&>(mesh_.objectRegistry::lookupObject<volVectorField>(mappingFieldsName_[fieldI]));
        vector fieldConstantValue_(dict_.lookup(neededFields_[fieldI]));

        forAll(field_, filedI) {
          field_[fieldI]= fieldConstantValue_;
        }
        forAll(field_.boundaryField(),patchI) {
          forAll(field_.boundaryField()[patchI], faceI) {
            field_.boundaryFieldRef()[patchI][faceI] = fieldConstantValue_;
          }
        }

      }

      if (mesh_.objectRegistry::foundObject<volTensorField>(neededFields_[fieldI])) {
        volTensorField& field_ = const_cast<volTensorField&>(mesh_.objectRegistry::lookupObject<volTensorField>(mappingFieldsName_[fieldI]));
        tensor fieldConstantValue_(dict_.lookup(neededFields_[fieldI]));

        forAll(field_, filedI) {
          field_[fieldI]= fieldConstantValue_;
        }
        forAll(field_.boundaryField(),patchI) {
          forAll(field_.boundaryField()[patchI], faceI) {
            field_.boundaryFieldRef()[patchI][faceI] = fieldConstantValue_;
          }
        }
      }
    }
  }

  forAll(currentTimePatchesDataFields_, patchI) {
    forAll(mappingFields_, fieldI) {
      currentTimePatchesDataFields_[patchI].append(-1);
    }
  }
}

inline wordList& Foam::simpleBoundaryMappingModel::mappingFieldsName() const
{
  return const_cast<wordList&>(mappingFieldsName_);
}

// ************************************************************************* //


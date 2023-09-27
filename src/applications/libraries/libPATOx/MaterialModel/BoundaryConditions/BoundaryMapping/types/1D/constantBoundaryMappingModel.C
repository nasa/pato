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

#include "constantBoundaryMappingModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantBoundaryMappingModel::constantBoundaryMappingModel
(
    const fvMesh& mesh,
    const wordList& neededFields,
    const dictionary dict
)
  :
simpleBoundaryMappingModel(mesh, neededFields, dict)
{
  listMappingFileData_.append(readFileData(mappingFileName_));
  Info << simpleModel::getTabLevel() << "constant Boundary Mapping: fields = ( ";
  forAll(mappingFieldsName_, fieldI) {
    Info <<   mappingFieldsName_[fieldI] << " ";
  }
  Info << ")" << endl;

  forAll(listMappingFileData_, fileI) {
    List<scalarList>& mappingFileData_ = listMappingFileData_[fileI];
    if (mappingFileData_.size()< 2) {
      FatalErrorInFunction << mappingFileName_ << " column size is less than 2.\n File format should be \"//time(s) field1 field2 ...\"" << exit(FatalError) ;
    }
  }

  simpleBoundaryMappingModel::init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantBoundaryMappingModel::~constantBoundaryMappingModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constantBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)
{
  if(foundInList(fieldName,mappingFieldsName_ )) {
    int indexField_=-1;
    forAll(mappingFieldsName_, fieldI) {
      if(fieldName == mappingFieldsName_[fieldI] ) {
        indexField_=fieldI;
        break;
      }
    }


    if (indexField_ < 0) {
      FatalErrorInFunction << fieldName << " not found in mappingFieldsName." << exit(FatalError);
    }

    label fieldI = indexField_;

    if (patchID<0) {
      FatalErrorInFunction << "patchID not correct." << exit(FatalError);
    }
    if (currentTimePatchesDataFields_[patchID][fieldI]==timeValue) {
      return;
    }
    Info << "update " << mappingFieldsName_[fieldI] << " from " << mappingFileName_  << endl;
    currentTimePatchesDataFields_[patchID][fieldI]=timeValue;

    List<scalarList>& mappingFileData_ = listMappingFileData_[0]; // only one file in constant type
    scalarList data_ = mappingFileData_[mappingFieldsColumn_[fieldI]];
    scalarList times_ = mappingFileData_[0];

    if (mesh_.objectRegistry::foundObject<volScalarField>(mappingFieldsName_[fieldI])) {
      volScalarField& field_ = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        field_.boundaryFieldRef()[patchID][faceI] = linearInterpolation(times_, data_, timeValue);
      }
    }

    if (mesh_.objectRegistry::foundObject<volVectorField>(mappingFieldsName_[fieldI])) {
      volVectorField& field_ = const_cast<volVectorField&>(mesh_.objectRegistry::lookupObject<volVectorField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        scalar value_ = linearInterpolation(times_, data_, timeValue);
        field_.boundaryFieldRef()[patchID][faceI].x() = value_;
        field_.boundaryFieldRef()[patchID][faceI].y() = value_;
        field_.boundaryFieldRef()[patchID][faceI].z() = value_;
      }
    }

    if (mesh_.objectRegistry::foundObject<volTensorField>(mappingFieldsName_[fieldI])) {
      volTensorField& field_ = const_cast<volTensorField&>(mesh_.objectRegistry::lookupObject<volTensorField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        scalar value_ = linearInterpolation(times_, data_, timeValue);
        field_.boundaryFieldRef()[patchID][faceI].xx() = value_;
        field_.boundaryFieldRef()[patchID][faceI].xy() = 0;
        field_.boundaryFieldRef()[patchID][faceI].xz() = 0;
        field_.boundaryFieldRef()[patchID][faceI].yx() = 0;
        field_.boundaryFieldRef()[patchID][faceI].yy() = value_;
        field_.boundaryFieldRef()[patchID][faceI].yz() = 0;
        field_.boundaryFieldRef()[patchID][faceI].zx() = 0;
        field_.boundaryFieldRef()[patchID][faceI].zy() = 0;
        field_.boundaryFieldRef()[patchID][faceI].zz() = value_;
      }
    }
  }
}

void Foam::constantBoundaryMappingModel::write(Ostream& os) const
{
  os.writeKeyword("mappingType") << "\"constant\"" << token::END_STATEMENT << nl;
  os.writeKeyword("mappingFileName") << mappingFileName_ << token::END_STATEMENT << nl;
  if (mappingFields_.size()>0) {
    os.writeKeyword("mappingFields")  << "(";
    forAll(mappingFieldsName_, fieldI) {
      os << "(" << mappingFieldsName_[fieldI] << " \"" << mappingFieldsColumn_[fieldI] << "\")";
    }
    os  << ")" << token::END_STATEMENT << nl;
    forAll(simpleBoundaryMappingModel::constantFields_, fieldI) {
      os.writeKeyword(simpleBoundaryMappingModel::constantFields_[fieldI].first())  << simpleBoundaryMappingModel::constantFields_[fieldI].second() <<  token::END_STATEMENT << nl;
    }
  }
}

// ************************************************************************* //

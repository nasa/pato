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

#include "sigmoidBoundaryMappingModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigmoidBoundaryMappingModel::sigmoidBoundaryMappingModel
(
    const fvMesh& mesh,
    const wordList& neededFields,
    const dictionary dict
)
  :
simpleBoundaryMappingModel(mesh, neededFields, dict)
{

  fileName mappingDictName_ = changeEnviVar(mappingFileName_); // change the environment variable
  IOdictionary mappingDict_
  (
      IOobject
      (
          mappingDictName_,
          mesh_,
          IOobject::MUST_READ,
          IOobject::NO_WRITE,
          false
      )
  );

  Info << "gaussian Boundary Mapping:" << endl;
  fieldCoeffs_.resize(mappingFieldsName_.size());
  wordList coefNames_; // Sigmoid coeffient names
  coefNames_.append("a");
  coefNames_.append("b");
  coefNames_.append("c");
  coefNames_.append("d");
  // f(t) = a + d / (1 + exp(b(x+c)))

  forAll(mappingFieldsName_, fieldI) {

    fieldCoeffs_[fieldI].resize(coefNames_.size());
    Info << mappingFieldsName_[fieldI] << "(t) = ";
    forAll(fieldCoeffs_[fieldI], coeffI) {
      fieldCoeffs_[fieldI][coeffI] = readScalar(mappingDict_.lookup(mappingFieldsName_[fieldI]+"["+coefNames_[coeffI]+"]"));
    }

    Info  << fieldCoeffs_[fieldI][0] << " + " << fieldCoeffs_[fieldI][3] << " / ( 1 + exp [ " << fieldCoeffs_[fieldI][1] << " ( t + " << fieldCoeffs_[fieldI][2] << " ) ] )" << endl;

  }

  // verify mappingFieldsColumn == 0 for type 0: e.g.  f(t) = a + d / (1 + exp(b(x+c)))
  forAll(mappingFieldsColumn_, columnI) {
    if(mappingFieldsColumn_[columnI] != 0) {
      FatalErrorInFunction << "Second value of the \"mappingFields\" defines the type of the gaussian. Only type 0 is available."
                           << nl << "Type 0: f(t) = a + d / (1 + exp(b(x+c)))" << exit(FatalError);
    }
  }
  simpleBoundaryMappingModel::init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigmoidBoundaryMappingModel::~sigmoidBoundaryMappingModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigmoidBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)
{
  if(foundInList(fieldName,mappingFieldsName_ )) {
    int indexField_=-1;
    forAll(mappingFieldsName_, fieldI) {
      if(fieldName == mappingFieldsName_[fieldI] ) {
        indexField_=fieldI;
        break;
      }
    }

    const label fieldI = indexField_;

    if (indexField_ < 0) {
      FatalErrorInFunction << fieldName << " not found in mappingFieldsName." << exit(FatalError);
    }
    if (patchID<0) {
      FatalErrorInFunction << "patchID not correct." << exit(FatalError);
    }
    if (currentTimePatchesDataFields_[patchID][fieldI]==timeValue) {
      return;
    }


    Info << "update " << mappingFieldsName_[fieldI] << " from " << mappingFileName_  << endl;
    currentTimePatchesDataFields_[patchID][fieldI]=timeValue;

    if (mesh_.objectRegistry::foundObject<volScalarField>(mappingFieldsName_[fieldI])) {
      volScalarField& field_ = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        field_.boundaryFieldRef()[patchID][faceI] = fieldCoeffs_[fieldI][0] + fieldCoeffs_[fieldI][3] / (1 + exp( fieldCoeffs_[fieldI][1] * (timeValue + fieldCoeffs_[fieldI][2])));
      }
    }

    if (mesh_.objectRegistry::foundObject<volVectorField>(mappingFieldsName_[fieldI])) {
      volVectorField& field_ = const_cast<volVectorField&>(mesh_.objectRegistry::lookupObject<volVectorField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        scalar value_ =  fieldCoeffs_[fieldI][0] + fieldCoeffs_[fieldI][3] / (1 + exp( fieldCoeffs_[fieldI][1] * (timeValue + fieldCoeffs_[fieldI][2])));
        field_.boundaryFieldRef()[patchID][faceI].x() = value_;
        field_.boundaryFieldRef()[patchID][faceI].y() = value_;
        field_.boundaryFieldRef()[patchID][faceI].z() = value_;
      }
    }

    if (mesh_.objectRegistry::foundObject<volTensorField>(mappingFieldsName_[fieldI])) {
      volTensorField& field_ = const_cast<volTensorField&>(mesh_.objectRegistry::lookupObject<volTensorField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        scalar value_ = fieldCoeffs_[fieldI][0] + fieldCoeffs_[fieldI][3] / (1 + exp( fieldCoeffs_[fieldI][1] * (timeValue + fieldCoeffs_[fieldI][2])));
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

void Foam::sigmoidBoundaryMappingModel::write(Ostream& os) const
{
  os.writeKeyword("mappingType") << "\"gaussian\"" << token::END_STATEMENT << nl;
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

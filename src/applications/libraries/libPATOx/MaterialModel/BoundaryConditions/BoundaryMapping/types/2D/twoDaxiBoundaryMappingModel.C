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

#include "twoDaxiBoundaryMappingModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoDaxiBoundaryMappingModel::twoDaxiBoundaryMappingModel
(
    const fvMesh& mesh,
    const wordList& neededFields,
    const dictionary dict
)
  :
simpleBoundaryMappingModel(mesh, neededFields, dict),
mappingAxis_(dict_.lookupOrDefault<word>("mappingAxis","none")),
mappingSymmetry_(dict_.lookupOrDefault<word>("mappingSymmetry","none"))
{
  Info << simpleModel::getTabLevel() << "2D-axi Boundary Mapping:" << endl;
  List<fileName> listFiles_ = filesInFolder(mappingFileName_.path());
  forAll(listFiles_, listI) {
    word name_= listFiles_[listI].name();
    int loc = name_.find(mappingFileName_.name()+"_");
    if (loc>=0) {
      word value_ = name_.replace(mappingFileName_.name()+"_","");
      if(isNumber(value_)) {
        timesMappingFileData_.append(stof(value_));
        listMappingFileData_.append(readFileData(listFiles_[listI]));
      }
    }
  }
  // Sort times and data
  labelList visitOrder;
  sortedOrder(timesMappingFileData_, visitOrder);
  timesMappingFileData_ = scalarList(timesMappingFileData_, visitOrder);
  listMappingFileData_ = List<List<scalarList> >(listMappingFileData_, visitOrder);

  wordList availableAxes_;
  availableAxes_.append("x");
  availableAxes_.append("y");
  availableAxes_.append("z");
  if (!foundInList(mappingAxis_,availableAxes_)) {
    FatalErrorInFunction << "mappingAxis from BoundaryModel is not valid. Please use:\n("<< endl ;

    forAll(availableAxes_, wordI) {
      FatalErrorInFunction << "  " << availableAxes_[wordI] << endl;
    }
    FatalErrorInFunction  << ")"     << exit(FatalError);
  }

  forAll(listMappingFileData_, fileI) {
    List<scalarList>& mappingFileData_ = listMappingFileData_[fileI];
    if (mappingFileData_.size()< 3) {
      FatalErrorInFunction << mappingFileName_ << " column size is less than 3.\n File format should be \"r z field1 field2 ... \"" << exit(FatalError) ;
    }
  }

  simpleBoundaryMappingModel::init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoDaxiBoundaryMappingModel::~twoDaxiBoundaryMappingModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoDaxiBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)
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
    Info << "update " << mappingFieldsName_[fieldI]<< " from \"" << (word) mappingFileName_ << "_*\""  << endl;
    currentTimePatchesDataFields_[patchID][fieldI]=timeValue;
    int index_time=-1;
    forAll(timesMappingFileData_, timeI) {
      if(timesMappingFileData_[timeI] >= timeValue) {
        index_time=timeI - 1;
        break;
      }
      if (timeI == timesMappingFileData_.size()-1) { // time value > last timesMappingFileData
        index_time = timesMappingFileData_.size()- 2;
      }
    }
    if (index_time<0) { // first
      index_time=0;
    }
    if (index_time == timesMappingFileData_.size()-1) { // last
      index_time = timesMappingFileData_.size()- 2; // last - 1
    }

    List<scalarList>& mappingFileDataI_ = listMappingFileData_[index_time];
    scalarList& radiusDataI_ = mappingFileDataI_[0];
    scalarList& zDataI_ = mappingFileDataI_[1];

    List<scalarList>& mappingFileDataIPlus1_ = listMappingFileData_[index_time+1];
    scalarList& radiusDataIPlus1_ = mappingFileDataIPlus1_[0];
    scalarList& zDataIPlus1_ = mappingFileDataIPlus1_[1];

    scalarList times_;
    times_.append(timesMappingFileData_[index_time]);
    times_.append(timesMappingFileData_[index_time+1]);
    if(debug_) {
      Info << "--- times_" << times_ << "---  Foam::BoundaryMapping::update2Daxi(scalar timeValue, label patchID, label fieldI)" << endl;
      Info << "--- radiusDataI_" << radiusDataI_ << "---  Foam::BoundaryMapping::update2Daxi(scalar timeValue, label patchID, label fieldI)" << endl;
      Info << "--- zDataI_" << zDataI_ << "---  Foam::BoundaryMapping::update2Daxi(scalar timeValue, label patchID, label fieldI)" << endl;
      Info << "--- radiusDataIPlus1_" << radiusDataIPlus1_ << "---  Foam::BoundaryMapping::update2Daxi(scalar timeValue, label patchID, label fieldI)" << endl;
      Info << "--- zDataIPlus1_" << zDataIPlus1_ << "---  Foam::BoundaryMapping::update2Daxi(scalar timeValue, label patchID, label fieldI)" << endl;
    }

    if (mesh_.objectRegistry::foundObject<volScalarField>(mappingFieldsName_[fieldI])) {
      volScalarField& field_ = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        scalar x = mesh_.boundaryMesh()[patchID].faceCentres()[faceI].x();
        scalar y = mesh_.boundaryMesh()[patchID].faceCentres()[faceI].y();
        scalar z = mesh_.boundaryMesh()[patchID].faceCentres()[faceI].z();
        scalar radius_ = -1;
        scalar z_ = -1;
        if (mappingAxis_=="x") {
          radius_ = sqrt(pow(y,2)+pow(z,2));
          z_ = x;
        }

        if (mappingAxis_=="y") {
          radius_ = sqrt(pow(x,2)+pow(z,2));
          z_ = y;
        }

        if (mappingAxis_=="z") {
          radius_ = sqrt(pow(x,2)+pow(y,2));
          z_ = z;
        }

        if (radius_ < 0 || z_ < 0) {
          FatalErrorInFunction << "mappingAxis from BoundaryModel is not valid." << exit(FatalError);
        }

        scalarList dataI_ = mappingFileDataI_[mappingFieldsColumn_[fieldI]];
        scalar interpolatedDataI_ = bilinearInterpolation(radiusDataI_, zDataI_, dataI_, radius_, z_);

        scalarList dataIPlus1_ = mappingFileDataIPlus1_[mappingFieldsColumn_[fieldI]];
        scalar interpolatedDataIPlus1_ = bilinearInterpolation(radiusDataIPlus1_, zDataIPlus1_, dataIPlus1_, radius_, z_);

        scalarList fieldsRZ_;
        fieldsRZ_.append(interpolatedDataI_);
        fieldsRZ_.append(interpolatedDataIPlus1_);
        field_.boundaryFieldRef()[patchID][faceI] = linearInterpolation(times_, fieldsRZ_, timeValue);

        if(debug_) {
          Info << "--- radius_ " << radius_ << "---  Foam::BoundaryMapping::update2Daxi(scalar timeValue, label patchID, label fieldI)" << endl;
          Info << "--- z_ " << z_ << "---  Foam::BoundaryMapping::update2Daxi(scalar timeValue, label patchID, label fieldI)" << endl;
          Info << "--- dataI_ " << dataI_ << "---  Foam::BoundaryMapping::update2Daxi(scalar timeValue, label patchID, label fieldI)" << endl;
          Info << "--- interpolatedDataI_ " << interpolatedDataI_ << "---  Foam::BoundaryMapping::update2Daxi(scalar timeValue, label patchID, label fieldI)" << endl;
          Info << "--- dataIPlus1_ " << dataIPlus1_ << "---  Foam::BoundaryMapping::update2Daxi(scalar timeValue, label patchID, label fieldI)" << endl;
          Info << "--- interpolatedDataIPlus1_ " << interpolatedDataIPlus1_ << "---  Foam::BoundaryMapping::update2Daxi(scalar timeValue, label patchID, label fieldI)" << endl;
          Info << "--- field_.boundaryFieldRef()[patchID][faceI]" << field_.boundaryFieldRef()[patchID][faceI] << "---  Foam::BoundaryMapping::update2Daxi(scalar timeValue, label patchID, label fieldI)" << endl;
        }

      }
    }

    if (mesh_.objectRegistry::foundObject<volVectorField>(mappingFieldsName_[fieldI])) {
      volVectorField& field_ = const_cast<volVectorField&>(mesh_.objectRegistry::lookupObject<volVectorField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        scalar x = mesh_.boundaryMesh()[patchID].faceCentres()[faceI].x();
        scalar y = mesh_.boundaryMesh()[patchID].faceCentres()[faceI].y();
        scalar z = mesh_.boundaryMesh()[patchID].faceCentres()[faceI].z();
        scalar radius_ = -1;
        scalar z_ = -1;
        if (mappingAxis_=="x") {
          radius_ = sqrt(pow(y,2)+pow(z,2));
          z_ = x;
        }

        if (mappingAxis_=="y") {
          radius_ = sqrt(pow(x,2)+pow(z,2));
          z_ = y;
        }

        if (mappingAxis_=="z") {
          radius_ = sqrt(pow(x,2)+pow(y,2));
          z_ = z;
        }

        if (radius_ < 0 || z_ < 0) {
          FatalErrorInFunction << "mappingAxis from BoundaryModel is not valid." << exit(FatalError);
        }

        scalarList dataI_ = mappingFileDataI_[mappingFieldsColumn_[fieldI]];
        scalar interpolatedDataI_ = bilinearInterpolation(radiusDataI_, zDataI_, dataI_, radius_, z_);

        scalarList dataIPlus1_ = mappingFileDataIPlus1_[mappingFieldsColumn_[fieldI]];
        scalar interpolatedDataIPlus1_ = bilinearInterpolation(radiusDataIPlus1_, zDataIPlus1_, dataIPlus1_, radius_, z_);

        scalarList fieldsRZ_;
        fieldsRZ_.append(interpolatedDataI_);
        fieldsRZ_.append(interpolatedDataIPlus1_);

        field_.boundaryFieldRef()[patchID][faceI].x() = linearInterpolation(times_, fieldsRZ_, timeValue);
        field_.boundaryFieldRef()[patchID][faceI].y() = field_.boundaryFieldRef()[patchID][faceI].x();
        field_.boundaryFieldRef()[patchID][faceI].z() = field_.boundaryFieldRef()[patchID][faceI].x();

      }
    }

    if (mesh_.objectRegistry::foundObject<volTensorField>(mappingFieldsName_[fieldI])) {
      volTensorField& field_ = const_cast<volTensorField&>(mesh_.objectRegistry::lookupObject<volTensorField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        scalar x = mesh_.boundaryMesh()[patchID].faceCentres()[faceI].x();
        scalar y = mesh_.boundaryMesh()[patchID].faceCentres()[faceI].y();
        scalar z = mesh_.boundaryMesh()[patchID].faceCentres()[faceI].z();
        scalar radius_ = -1;
        scalar z_ = -1;
        if (mappingAxis_=="x") {
          radius_ = sqrt(pow(y,2)+pow(z,2));
          z_ = x;
        }

        if (mappingAxis_=="y") {
          radius_ = sqrt(pow(x,2)+pow(z,2));
          z_ = y;
        }

        if (mappingAxis_=="z") {
          radius_ = sqrt(pow(x,2)+pow(y,2));
          z_ = z;
        }

        if (radius_ < 0 || z_ < 0) {
          FatalErrorInFunction << "mappingAxis from BoundaryModel is not valid." << exit(FatalError);
        }

        scalarList dataI_ = mappingFileDataI_[mappingFieldsColumn_[fieldI]];
        scalar interpolatedDataI_ = bilinearInterpolation(radiusDataI_, zDataI_, dataI_, radius_, z_);

        scalarList dataIPlus1_ = mappingFileDataIPlus1_[mappingFieldsColumn_[fieldI]];
        scalar interpolatedDataIPlus1_ = bilinearInterpolation(radiusDataIPlus1_, zDataIPlus1_, dataIPlus1_, radius_, z_);

        scalarList fieldsRZ_;
        fieldsRZ_.append(interpolatedDataI_);
        fieldsRZ_.append(interpolatedDataIPlus1_);

        field_.boundaryFieldRef()[patchID][faceI].xx() = linearInterpolation(times_, fieldsRZ_, timeValue);
        field_.boundaryFieldRef()[patchID][faceI].yy() = field_.boundaryFieldRef()[patchID][faceI].xx();
        field_.boundaryFieldRef()[patchID][faceI].zz() = field_.boundaryFieldRef()[patchID][faceI].xx();
        field_.boundaryFieldRef()[patchID][faceI].xy() = 0;
        field_.boundaryFieldRef()[patchID][faceI].xz() = 0;
        field_.boundaryFieldRef()[patchID][faceI].yx() = 0;
        field_.boundaryFieldRef()[patchID][faceI].yz() = 0;
        field_.boundaryFieldRef()[patchID][faceI].zx() = 0;
        field_.boundaryFieldRef()[patchID][faceI].zy() = 0;

      }
    }
  }
}

void Foam::twoDaxiBoundaryMappingModel::write(Ostream& os) const
{
  os.writeKeyword("mappingType") << "\"2D-axi\"" << token::END_STATEMENT << nl;
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
  if (mappingAxis_!="none") {
    os.writeKeyword("mappingAxis") << mappingAxis_ << token::END_STATEMENT << nl;
  }
  if (mappingSymmetry_!="none") {
    os.writeKeyword("mappingSymmetry") << mappingSymmetry_ << token::END_STATEMENT << nl;
  }
}

// ************************************************************************* //

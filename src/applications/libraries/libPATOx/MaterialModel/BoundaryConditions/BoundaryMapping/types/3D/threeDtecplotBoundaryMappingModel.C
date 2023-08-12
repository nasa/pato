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

#include "threeDtecplotBoundaryMappingModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::threeDtecplotBoundaryMappingModel::threeDtecplotBoundaryMappingModel
(
    const fvMesh& mesh,
    const wordList& neededFields,
    const dictionary dict
)
  :
simpleBoundaryMappingModel(mesh, neededFields, dict),
mappingAxis_(dict_.lookupOrDefault<word>("mappingAxis","none")),
mappingSymmetry_(dict_.lookupOrDefault<word>("mappingSymmetry","none")),
updateTolTime_(dict_.lookupOrDefault<scalar>("updateTolTime",0.1)),
multiplier_(dict_.lookupOrDefault<scalar>("multiplier",1))
{
  Info << simpleModel::getTabLevel() << "3D-tecplot Boundary Mapping:" << endl;

  List<fileName> listFiles_ = filesInFolder(mappingFileName_.path());
  forAll(listFiles_, listI) {
    word name_= listFiles_[listI].name();
    int loc = name_.find(mappingFileName_.name()+"_");
    if (loc>=0) {
      word value_ = name_.replace(mappingFileName_.name()+"_","");
      if(isNumber(value_)) {
        timesMappingFileData_.append(stof(value_));
        listMappingFileData_.append(readTecplotFileData(listFiles_[listI]));
      }
    }
  }

  if (timesMappingFileData_.size()==0) {
    FatalErrorInFunction <<  "3D-tecplot method: input files not found using \"mappingFileName\".\nUsage example: mappingFileName \"" << (word) mappingFileName_ << "_10\"; //input file name for time = 10s"<< exit(FatalError);
  }
  if (timesMappingFileData_.size()<2) {
    FatalErrorInFunction <<  "3D-tecplot method needs at least 2 input files using \"mappingFileName\".\nUsage example: mappingFileName \"" << (word) mappingFileName_ << "_10\"; //input file name for time = 10s"<< exit(FatalError);
  }
  // Sort times and data
  labelList visitOrder;
  sortedOrder(timesMappingFileData_, visitOrder);
  timesMappingFileData_ = scalarList(timesMappingFileData_, visitOrder);
  listMappingFileData_ = List<List<scalarList> >(listMappingFileData_, visitOrder);

  wordList availableSymmetry_;
  availableSymmetry_.append("x");
  availableSymmetry_.append("y");
  availableSymmetry_.append("z");
  availableSymmetry_.append("none");
  if (!foundInList(mappingAxis_,availableSymmetry_)) {
    FatalErrorInFunction << "mappingSymmetry from BoundaryModel is not valid. Please use:\n("<< endl ;

    forAll(availableSymmetry_, wordI) {
      FatalErrorInFunction << "  " << availableSymmetry_[wordI] << endl;
    }
    FatalErrorInFunction  << ")"     << exit(FatalError);
  }

  forAll(listMappingFileData_, fileI) {
    List<scalarList>& mappingFileData_ = listMappingFileData_[fileI];
    if (mappingFileData_.size()< 4) {
      FatalErrorInFunction << mappingFileName_ << " column size is less than 4.\n The Tecplot file must have at least 4 variables: x,y,z,field1 \"" << exit(FatalError) ;
    }
  }
  simpleBoundaryMappingModel::init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::threeDtecplotBoundaryMappingModel::~threeDtecplotBoundaryMappingModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::threeDtecplotBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)
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
    if (mag(currentTimePatchesDataFields_[patchID][fieldI]-timeValue) < updateTolTime_) {
      return;
    }

    Info << "update " << mappingFieldsName_[fieldI] << " from \"" <<  (word) mappingFileName_ << "_*\""<<  endl;

    if(debug_) {
      Info << "--- currentTimePatchesDataFields_ ---  Foam::threeDtecplotBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)" << endl;
    }
    currentTimePatchesDataFields_[patchID][fieldI]=timeValue;
    int index_time=-1;

    if(debug_) {
      Info << "--- timesMappingFileData_ ---  Foam::threeDtecplotBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)" << endl;
    }
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

    if(debug_) {
      Info << "--- mappingFileDataI_ ---  Foam::threeDtecplotBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)" << endl;
    }
    List<scalarList>& mappingFileDataI_ = listMappingFileData_[index_time]; // listMappingFileData_[time][3(xyz)+fieldI][dataI(xyz)]
    List<scalarList>& mappingFileDataIPlus1_ = listMappingFileData_[index_time+1];

    if(debug_) {
      Info << "--- times_ ---  Foam::threeDtecplotBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)" << endl;
      Info << "--- timesMappingFileData_ = " << timesMappingFileData_ << " --- Foam::threeDtecplotBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)" << endl;
      Info << "--- timeValue = " << timeValue << " --- Foam::threeDtecplotBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)" << endl;
      Info << "--- index_time = " << index_time << " --- Foam::threeDtecplotBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)" << endl;
    }

    scalarList times_;
    times_.append(timesMappingFileData_[index_time]);
    times_.append(timesMappingFileData_[index_time+1]);

    if(debug_) {
      Info << "--- mesh_.objectRegistry::foundObject<volScalarField>(mappingFieldsName_[fieldI]) ---  Foam::threeDtecplotBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)" << endl;
    }
    if (mesh_.objectRegistry::foundObject<volScalarField>(mappingFieldsName_[fieldI])) {
      volScalarField& field_ = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        scalarList fieldsInterp_;
        if(debug_) {
          Info << "--- t = " << timesMappingFileData_[index_time] << ", face I = " << faceI << " ---  Foam::BoundaryMapping::interpolateNearest(label patchID, label faceI, label fieldI, List<scalarList>& mappingFileData)" << endl;
        }
        fieldsInterp_.append(interpolateNearest(patchID,faceI, fieldI, mappingFileDataI_));
        if(debug_) {
          Info << "--- t = " << timesMappingFileData_[index_time+1] << ", face I = " << faceI << " ---  Foam::BoundaryMapping::interpolateNearest(label patchID, label faceI, label fieldI, List<scalarList>& mappingFileData)" << endl;
        }
        fieldsInterp_.append(interpolateNearest(patchID,faceI, fieldI, mappingFileDataIPlus1_));
        field_.boundaryFieldRef()[patchID][faceI] = linearInterpolation(times_, fieldsInterp_, timeValue)*multiplier_;
      }
    }

    if(debug_) {
      Info << "--- mesh_.objectRegistry::foundObject<volVectorField>(mappingFieldsName_[fieldI]) ---  Foam::threeDtecplotBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)" << endl;
    }
    if (mesh_.objectRegistry::foundObject<volVectorField>(mappingFieldsName_[fieldI])) {
      volVectorField& field_ = const_cast<volVectorField&>(mesh_.objectRegistry::lookupObject<volVectorField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        scalarList fieldsInterp_;
        fieldsInterp_.append(interpolateNearest(patchID,faceI, fieldI, mappingFileDataI_));
        fieldsInterp_.append(interpolateNearest(patchID,faceI, fieldI, mappingFileDataIPlus1_));
        field_.boundaryFieldRef()[patchID][faceI].x() = linearInterpolation(times_, fieldsInterp_, timeValue)*multiplier_;
        field_.boundaryFieldRef()[patchID][faceI].y() = field_.boundaryFieldRef()[patchID][faceI].x();
        field_.boundaryFieldRef()[patchID][faceI].z() = field_.boundaryFieldRef()[patchID][faceI].x();
      }
    }

    if(debug_) {
      Info << "--- mesh_.objectRegistry::foundObject<volTensorField>(mappingFieldsName_[fieldI]) ---  Foam::threeDtecplotBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)" << endl;
    }
    if (mesh_.objectRegistry::foundObject<volTensorField>(mappingFieldsName_[fieldI])) {
      volTensorField& field_ = const_cast<volTensorField&>(mesh_.objectRegistry::lookupObject<volTensorField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        scalarList fieldsInterp_;
        fieldsInterp_.append(interpolateNearest(patchID,faceI, fieldI, mappingFileDataI_));
        fieldsInterp_.append(interpolateNearest(patchID,faceI, fieldI, mappingFileDataIPlus1_));
        field_.boundaryFieldRef()[patchID][faceI].xx() = linearInterpolation(times_, fieldsInterp_, timeValue)*multiplier_;
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


    if(debug_) {
      Info << "--- end ---  Foam::threeDtecplotBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)" << endl;
    }
  }
}

Foam::scalar Foam::threeDtecplotBoundaryMappingModel::interpolateNearest(label patchID, label faceI, label fieldI, List<scalarList>& mappingFileData)
{
  if(debug_) {
    Info << "--- start ---  Foam::BoundaryMapping::interpolateNearest(label patchID, label faceI, label fieldI, List<scalarList>& mappingFileData)" << endl;
  }
  if (mappingFileData.size()< 4) {
    FatalErrorInFunction << mappingFileName_ << " column size is less than 4.\n The Tecplot file must have at least 4 variables: x,y,z,field1 \"" << exit(FatalError) ;
  }
  if (mappingFieldsColumn_[fieldI] < 3) {
    FatalErrorInFunction <<  mappingFields_[fieldI] << " in mappingFields is not correct. ColumnI has to be more than 2." << exit(FatalError) ;
  }

  scalar value_ = 0;
  int num_i = 0;
  forAll(mesh_.boundaryMesh()[patchID][faceI], pI) {
    vector pointI = mesh_.points()[mesh_.boundaryMesh()[patchID][faceI][pI]];
    double xp = pointI[0];
    double yp = pointI[1];
    double zp = pointI[2];

    scalar distance_min = -1;
    label indexData = -1;

    scalarList xFileData_ = mappingFileData[0];
    scalarList yFileData_ = mappingFileData[1];
    scalarList zFileData_ = mappingFileData[2];

    forAll(mappingFileData[mappingFieldsColumn_[fieldI]], dataI) {
      double x = xFileData_[dataI];
      double y = yFileData_[dataI];
      double z = zFileData_[dataI];

      double dist = ::sqrt(::pow(x - xp, 2) + ::pow(y - yp, 2) + ::pow(z - zp, 2)) ;

      if (mappingSymmetry_!="none") {
        double symmetry_ = 0;
        double dist_symm = 0;
        if (mappingSymmetry_ == "x") {
          symmetry_ = -x;
          dist_symm = ::sqrt(::pow(symmetry_ - xp, 2) + ::pow(y - yp, 2) + ::pow(z - zp, 2)) ;
        }
        if (mappingSymmetry_ == "y") {
          symmetry_ = -y;
          dist_symm = ::sqrt(::pow(x - xp, 2) + ::pow(symmetry_ - yp, 2) + ::pow(z - zp, 2)) ;
        }
        if (mappingSymmetry_ == "z") {
          symmetry_ = -z;
          dist_symm = ::sqrt(::pow(x - xp, 2) + ::pow(y - yp, 2) + ::pow(symmetry_ - zp, 2)) ;
        }

        if (dist_symm < dist) {
          dist = dist_symm;
        }
      }

      if (distance_min < 0 || dist < distance_min) {
        distance_min = dist;
        indexData = dataI;
      }
    }

    if (indexData >= 0) {
      value_ += mappingFileData[mappingFieldsColumn_[fieldI]][indexData];
      if(debug_) {
        Info << "--- (xp=" << xp << ",yp=" << yp << ",zp=" << zp << ") => " << mappingFields_[fieldI].first() << "(" << xFileData_[indexData] << "," << yFileData_[indexData]  << "," << zFileData_[indexData] << ") = " <<  mappingFileData[mappingFieldsColumn_[fieldI]][indexData] <<
             " ---  Foam::BoundaryMapping::interpolateNearest(label patchID, label faceI, label fieldI, List<scalarList>& mappingFileData)" << endl;
      }
      num_i++;
    }
  }

  if (num_i != 0) {
    value_ /= num_i;
  }
  if(debug_) {
    Info << "--- end ---  Foam::BoundaryMapping::interpolateNearest(label patchID, label faceI, label fieldI, List<scalarList>& mappingFileData)" << endl;
  }

  return value_;
}


void Foam::threeDtecplotBoundaryMappingModel::write(Ostream& os) const
{
  os.writeKeyword("mappingType") << "\"3D-tecplot\"" << token::END_STATEMENT << nl;
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
  os.writeKeyword("updateTolTime") << updateTolTime_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //

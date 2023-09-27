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

#include "twoDaxiFluxMapTimeBoundaryMappingModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoDaxiFluxMapTimeBoundaryMappingModel::twoDaxiFluxMapTimeBoundaryMappingModel
(
    const fvMesh& mesh,
    const wordList& neededFields,
    const dictionary dict
)
  :
simpleBoundaryMappingModel(mesh, neededFields, dict),
energyModel(meshLookupOrConstructModel<simpleEnergyModel>(mesh,mesh.name(),simpleEnergyModel::modelName)),
fluxFactorMapFileName_(dict.lookup("fluxFactorMapFileName"))
{
  // Read the Boundary Conditions file with the maximum field values in function of time
  listMappingFileData_.append(readFileData(mappingFileName_));
  // Print the Boundary Mapping fields
  Info << simpleModel::getTabLevel() << "2D-axi_fluxMap Boundary Mapping: fields = ( ";
  forAll(mappingFieldsName_, fieldI) {
    Info <<   mappingFieldsName_[fieldI] << " ";
  }
  Info << ")" << endl;
  // Read the flux factor files (file name + "_" + time)
  scalarList indexes;
  List<fileName> listFiles_ = filesInFolder(fluxFactorMapFileName_.path());
  forAll(listFiles_, listI) {
    word name_= listFiles_[listI].name();
    int loc = name_.find(fluxFactorMapFileName_.name()+"_");
    if (loc>=0) {
      word value_ = name_.replace(fluxFactorMapFileName_.name()+"_","");
      if(isNumber(value_)) {
        fluxFactor_times.append(stof(value_));
        indexes.append(listI);
      }
    }
  }
  labelList visitOrder;
  sortedOrder(fluxFactor_times, visitOrder);   // Sort the times
  fluxFactor_times = scalarList(fluxFactor_times, visitOrder);
  indexes = scalarList(indexes, visitOrder);
  forAll(indexes, i) {
    dictionary new_dict=dict;
    new_dict.set("fluxFactorMapFileName",listFiles_[(int) indexes[i]]);
    word fluxMapName = "fluxMap_"+listFiles_[(int) indexes[i]].name();
    word pressureMapName = "empty_factor";
    fluxFactor_list.append(new fluxFactor(new_dict, energyModel, fluxMapName, pressureMapName));
  }
  if (fluxFactor_list.size()==0) {
    FatalErrorInFunction << "fluxFactor_list.size()==0" << exit(FatalError);
  }
  dictionary new_dict=dict;
  new_dict.set("fluxFactorMapFileName",listFiles_[(int) indexes[0]]);
  fileName fluxFactorMapFileName=dict.lookup("fluxFactorMapFileName");
  word fluxMapName = "fluxMap_"+fluxFactorMapFileName.name();
  word pressureMapName = "empty_factor";
  fluxFactor_ptr = new fluxFactor(new_dict, energyModel, fluxMapName, pressureMapName);

  // Verify the listMappingFileData_
  forAll(listMappingFileData_, fileI) {
    List<scalarList>& mappingFileData_ = listMappingFileData_[fileI];

    if (mappingFileData_.size()< 2) {
      FatalErrorInFunction << mappingFileName_ << " column size is less than 2.\n File format should be \"//time(s) field1 field2 ...\"" << exit(FatalError);
    }
  }
  simpleBoundaryMappingModel::init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoDaxiFluxMapTimeBoundaryMappingModel::~twoDaxiFluxMapTimeBoundaryMappingModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoDaxiFluxMapTimeBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)
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
    Info << "update " << mappingFieldsName_[fieldI]<< " from " <<  mappingFileName_ << " and \"" << (word) fluxFactorMapFileName_ <<  "_*\"" << endl;
    currentTimePatchesDataFields_[patchID][fieldI]=timeValue;

    // Set the time index
    label index_time=-1;
    forAll(fluxFactor_times, timeI) {
      if(fluxFactor_times[timeI] >= timeValue) {
        index_time=timeI - 1;
        break;
      }
      if (timeI == fluxFactor_times.size()-1) { // time value > last timesMappingFileData
        index_time = fluxFactor_times.size()- 2;
      }
    }
    if (index_time<0) { // first
      index_time=0;
    }
    if (index_time == fluxFactor_times.size()-1) { // last
      index_time = fluxFactor_times.size()- 2; // last - 1
    }

    // Flux factor data and times
    const volScalarField& fluxFactor_dataI = fluxFactor_list[index_time].fluxMap();
    const volScalarField& fluxFactor_dataIPlus1 = fluxFactor_list[index_time+1].fluxMap();
    scalarList f_times_;
    f_times_.append(fluxFactor_times[index_time]);
    f_times_.append(fluxFactor_times[index_time+1]);
    // Update the flux factor
    volScalarField& factor_ = fluxFactor_ptr->fluxMap();
    forAll(factor_.boundaryField()[patchID], faceI) {
      scalarList factor_data;
      factor_data.append(fluxFactor_dataI.boundaryField()[patchID][faceI]);
      factor_data.append(fluxFactor_dataIPlus1.boundaryField()[patchID][faceI]);
      factor_.boundaryFieldRef()[patchID][faceI] = linearInterpolation(f_times_, factor_data, timeValue);
    }

    List<scalarList>& mappingFileData_ = listMappingFileData_[0]; // only one file in constant type
    scalarList data_ = mappingFileData_[mappingFieldsColumn_[fieldI]];
    scalarList times_ = mappingFileData_[0];

    if (mesh_.objectRegistry::foundObject<volScalarField>(mappingFieldsName_[fieldI])) {
      volScalarField& field_ = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        field_.boundaryFieldRef()[patchID][faceI] = linearInterpolation(times_, data_, timeValue)*factor_.boundaryField()[patchID][faceI];
      }
    }

    if (mesh_.objectRegistry::foundObject<volVectorField>(mappingFieldsName_[fieldI])) {
      volVectorField& field_ = const_cast<volVectorField&>(mesh_.objectRegistry::lookupObject<volVectorField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        scalar value_ = linearInterpolation(times_, data_, timeValue)*factor_.boundaryField()[patchID][faceI];
        field_.boundaryFieldRef()[patchID][faceI].x() = value_;
        field_.boundaryFieldRef()[patchID][faceI].y() = value_;
        field_.boundaryFieldRef()[patchID][faceI].z() = value_;
      }
    }

    if (mesh_.objectRegistry::foundObject<volTensorField>(mappingFieldsName_[fieldI])) {
      volTensorField& field_ = const_cast<volTensorField&>(mesh_.objectRegistry::lookupObject<volTensorField>(mappingFieldsName_[fieldI]));
      forAll(field_.boundaryField()[patchID], faceI) {
        scalar value_ = linearInterpolation(times_, data_, timeValue)*factor_.boundaryField()[patchID][faceI];
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

void Foam::twoDaxiFluxMapTimeBoundaryMappingModel::write(Ostream& os) const
{
  os.writeKeyword("mappingType") << "\"2D-axi_fluxMapTime\"" << token::END_STATEMENT << nl;
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
  os.writeKeyword("fluxFactorNormal") << fluxFactor_ptr->fluxFactorNormal() << token::END_STATEMENT << nl;
  os.writeKeyword("fluxFactorCenter") << fluxFactor_ptr->fluxFactorCenter() << token::END_STATEMENT << nl;
  os.writeKeyword("fluxFactorProjection") << fluxFactor_ptr->fluxFactorProjection() << token::END_STATEMENT << nl;
  os.writeKeyword("fluxFactorMapFileName") << fluxFactorMapFileName_ << token::END_STATEMENT << nl;
  os.writeKeyword("pointMotionDirection") << fluxFactor_ptr->pointMotionDirection() << token::END_STATEMENT << nl;
  os.writeKeyword("fluxFactorThreshold") << fluxFactor_ptr->fluxFactorThreshold() << token::END_STATEMENT << nl;
  os.writeKeyword("moreThanThresholdFlag") << fluxFactor_ptr->moreThanThresholdFlag() << token::END_STATEMENT << nl;
}
// ************************************************************************* //

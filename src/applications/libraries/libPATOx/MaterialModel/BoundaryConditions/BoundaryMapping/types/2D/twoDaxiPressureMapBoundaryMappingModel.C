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

#include "twoDaxiPressureMapBoundaryMappingModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoDaxiPressureMapBoundaryMappingModel::twoDaxiPressureMapBoundaryMappingModel
(
    const fvMesh& mesh,
    const wordList& neededFields,
    const dictionary dict
)
  :
simpleBoundaryMappingModel(mesh, neededFields, dict),
massModel(meshLookupOrConstructModel<simpleMassModel>(mesh,mesh.name(),simpleMassModel::modelName)),
dynamicPressureFieldName_(dict_.lookup("dynamicPressureFieldName")),
p_dyn_(massModel.createVolField<scalar>(dynamicPressureFieldName_))
{
  listMappingFileData_.append(readFileData(mappingFileName_));
  Info << simpleModel::getTabLevel() << "2D-axi_pressureMap Boundary Mapping: fields = ( ";
  forAll(mappingFieldsName_, fieldI) {
    Info <<   mappingFieldsName_[fieldI] << " ";
  }
  Info << ")" << endl;
  fileName fluxFactorMapFileName=dict.lookup("fluxFactorMapFileName");
  word fluxMapName = "empty_flux_factor";
  word pressureMapName = "pressureMap_"+fluxFactorMapFileName.name();
  fluxFactor_ptr = new fluxFactor(dict, massModel, fluxMapName, pressureMapName);

  forAll(listMappingFileData_, fileI) {
    List<scalarList>& mappingFileData_ = listMappingFileData_[fileI];
    if (mappingFileData_.size()< 2) {
      FatalErrorInFunction << mappingFileName_ << " column size is less than 2.\n File format should be \"//time(s) field1 field2 ...\"" << exit(FatalError) ;
    }
  }

  simpleBoundaryMappingModel::init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoDaxiPressureMapBoundaryMappingModel::~twoDaxiPressureMapBoundaryMappingModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoDaxiPressureMapBoundaryMappingModel::update(scalar timeValue, label patchID, word fieldName)
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

    Info << "update " << mappingFieldsName_[fieldI]<< " from " <<  mappingFileName_ << " and " << fluxFactor_ptr->fluxFactorMapFileName() << endl;
    p_dyn_.correctBoundaryConditions();
    currentTimePatchesDataFields_[patchID][fieldI]=timeValue;
    fluxFactor& fluxFactor_ = *fluxFactor_ptr;
    volScalarField factor_ =  fluxFactor_.pressureMap();
    List<scalarList>& mappingFileData_ = listMappingFileData_[0]; // only one file in constant type
    scalarList data_ = mappingFileData_[mappingFieldsColumn_[fieldI]];
    scalarList times_ = mappingFileData_[0];

    if (mesh_.objectRegistry::foundObject<volScalarField>(mappingFieldsName_[fieldI])) {
      volScalarField& field_ = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(mappingFieldsName_[fieldI]));
      Info << mappingFieldsName_[fieldI] << endl;
      forAll(field_.boundaryField()[patchID], faceI) {
        scalar p_bf = linearInterpolation(times_, data_, timeValue);
        scalar p_dyn_bf = p_dyn_.boundaryField()[patchID][faceI];
        field_.boundaryFieldRef()[patchID][faceI] = (p_bf - p_dyn_bf) + p_dyn_bf*factor_.boundaryField()[patchID][faceI];
      }
    } else {

      FatalErrorInFunction << "2Daxi-PressureMap not implemented for vector or tensor" << exit(FatalError);
    }

  }
}

void Foam::twoDaxiPressureMapBoundaryMappingModel::write(Ostream& os) const
{
  os.writeKeyword("mappingType") << "\"2D-axi_pressureMap\"" << token::END_STATEMENT << nl;
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
  os.writeKeyword("fluxFactorMapFileName") << fluxFactor_ptr->fluxFactorMapFileName() << token::END_STATEMENT << nl;
  os.writeKeyword("pointMotionDirection") << fluxFactor_ptr->pointMotionDirection() << token::END_STATEMENT << nl;
  os.writeKeyword("fluxFactorThreshold") << fluxFactor_ptr->fluxFactorThreshold() << token::END_STATEMENT << nl;
  os.writeKeyword("moreThanThresholdFlag") << fluxFactor_ptr->moreThanThresholdFlag() << token::END_STATEMENT << nl;
  os.writeKeyword("dynamicPressureFieldName") << dynamicPressureFieldName_ << token::END_STATEMENT << nl;
}
// ************************************************************************* //

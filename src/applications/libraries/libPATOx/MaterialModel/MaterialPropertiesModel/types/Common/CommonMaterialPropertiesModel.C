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
    along with OpenFOAM.  If Porous_t, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "CommonMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CommonMaterialPropertiesModel::CommonMaterialPropertiesModel(simpleModel& simpleModel, const word& modelName, const word& typeName)
  :
simpleModel_(simpleModel),
modelName_(modelName),
typeName_(typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CommonMaterialPropertiesModel::~CommonMaterialPropertiesModel()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CommonMaterialPropertiesModel::initialize()
{
  word bold_on="\e[1m";
  word bold_off="\e[0m";
  const fvMesh& mesh_ = simpleModel_.mesh();
  Info << simpleModel_.getTabLevel("| ") << bold_on << "Initialize " << (fileName) modelName_
       << " type " << (fileName) typeName_ <<  bold_off << endl;
  simpleModel_.tabLevel_++;
  forAll(listMapFields,i) {
    word fieldName = listMapFields[i].fieldName_;
    word typeName = listMapFields[i].fieldTypeName_;
    mapFields[typeName+"("+fieldName+")"]=&listMapFields[i];
    int count=0;
    if (typeName=="volScalarField") {
      if(mesh_.objectRegistry::foundObject<volScalarField>(fieldName)) {
        activeFieldIndexes_.append(i);
        scalarFields_.insert(fieldName,meshLookup<volScalarField>(mesh_,fieldName,""));
        listMapFields[i].f_init_();
        count++;
      }
    }
    if (typeName=="volVectorField") {
      if(mesh_.objectRegistry::foundObject<volVectorField>(fieldName)) {
        activeFieldIndexes_.append(i);
        vectorFields_.insert(fieldName,meshLookup<volVectorField>(mesh_,fieldName,""));
        listMapFields[i].f_init_();
        count++;
      }
    }
    if (typeName=="volTensorField") {
      if(mesh_.objectRegistry::foundObject<volTensorField>(fieldName)) {
        activeFieldIndexes_.append(i);
        tensorFields_.insert(fieldName,meshLookup<volTensorField>(mesh_,fieldName,""));
        listMapFields[i].f_init_();
        count++;
      }
    }
    if (count==0) {
      Info << simpleModel_.getTabLevel() << (fileName) fieldName << " not found in mesh. "
           << "It will not be updated." << endl;
    }
    count=0;
  }
  simpleModel_.tabLevel_--;
}

void Foam::CommonMaterialPropertiesModel::update()
{
  if (simpleModel_.debug()) {
    word bold_on="\e[1m";
    word bold_off="\e[0m";
    Info << simpleModel_.getTabLevel("| ") << bold_on << "Update " << (fileName) modelName_
         << " type " << (fileName) typeName_ <<  bold_off << endl;
    simpleModel_.tabLevel_++;
  }
  forAll(activeFieldIndexes_,i) {
    if (simpleModel_.debug()) {
      Info << simpleModel_.getTabLevel() << "Update " << (fileName) listMapFields[activeFieldIndexes_[i]].fieldName_ << endl;
    }
    listMapFields[activeFieldIndexes_[i]].f_update_();
  }
  if (simpleModel_.debug()) simpleModel_.tabLevel_--;
}



// ************************************************************************* //

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

#include "Porous_factorMaterialPropertiesModel.H"

#include "LinearInterpolation_factor_MaterialPropertiesObject.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Porous_factorMaterialPropertiesModel::Porous_factorMaterialPropertiesModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
PorousMaterialPropertiesModel(mesh, dictName),
startModelInit_(startModelInit()),
fieldFactors_
(
    simpleMaterialPropertiesModel
    ::materialDict_.subDict
    (
        "MaterialProperties"
    ).lookupOrDefault<List<Tuple2<word,scalar>>>
    (
        "fieldFactors",
        List<Tuple2<word,scalar>>(0)
    )

)
{
  Info << getTabLevel() << "fieldFactors = (" << nl;
  forAll(fieldFactors_, i) {
    Info << getTabLevel() << fieldFactors_[i] << endl;
  }
  Info << getTabLevel() << ")" << endl;
  wordList matPropFieldNames = {
    "cp_v", "h_v", "ki_v", "kj_v", "kk_v", "emissivity_v",
    "absorptivity_v", "cp_c", "h_c", "ki_c", "kj_c", "kk_c",
    "emissivity_c", "absorptivity_c"
  };

  // Initialize virgin and char factor lists
  int factorList_size = matPropFieldNames.size()/2;
  for (int i = 0; i < factorList_size; i++ ) {
    factorList_char_.append(1);
    factorList_virgin_.append(1);
  }

  // Check that correct fields are entered
  forAll(fieldFactors_, i) {
    if (!foundInList(fieldFactors_[i].first(), matPropFieldNames)) {
      FatalError
          << fieldFactors_[i].first()
          << " field not found in the material properties fields: "
          << matPropFieldNames << exit(FatalError);
    }
  }

  // Assign factors to virgin list in correct order
  forAll(factorList_virgin_, i) {
    forAll(fieldFactors_, j) {
      if (fieldFactors_[j].first() == matPropFieldNames[i]) {
        factorList_virgin_[i] = fieldFactors_[j].second();
      }
    }
  }

  // Assign factors to char list in correct order
  forAll(factorList_char_, i) {
    forAll(fieldFactors_, j) {
      if
      (
          fieldFactors_[j].first()
          == matPropFieldNames[i+factorList_virgin_.size()]
      ) {
        factorList_char_[i] = fieldFactors_[j].second();
      }
    }
  }
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Porous_factorMaterialPropertiesModel::~Porous_factorMaterialPropertiesModel()
{}

// * * * * * * * * * * * * * * * * Memeber Functions  * * * * * * * * * * * * * * * //

void Foam::Porous_factorMaterialPropertiesModel::initialize()
{
  PorousMaterialPropertiesModel::initialize();

  // Create new virgin and char objects with factors applied
  virginObject.reset
  (
      new LinearInterpolation_Factor_MaterialPropertiesObject
      (
          readCharredVirginDict("virgin"),
          factorList_virgin_
      )
  );
  charObject.reset
  (
      new LinearInterpolation_Factor_MaterialPropertiesObject
      (
          readCharredVirginDict("char"),
          factorList_char_
      )
  );
}

void Foam::Porous_factorMaterialPropertiesModel::update()
{
  PorousMaterialPropertiesModel::update();
}


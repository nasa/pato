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

#include "simpleModel.H"


// * * * * * * * * * * * * * * * * Protected Static Members  * * * * * * * * * * * * * * //

int Foam::simpleModel::tabLevel_=0;
Foam::PMap<Foam::wordList> Foam::simpleModel::wordPropNames_;
Foam::PMap<Foam::wordList> Foam::simpleModel::fileNamePropNames_;
Foam::PMap<Foam::wordList> Foam::simpleModel::switchPropNames_;
Foam::PMap<Foam::wordList> Foam::simpleModel::scalarPropNames_;
Foam::PMap<Foam::wordList> Foam::simpleModel::vectorPropNames_;
Foam::PMap<Foam::wordList> Foam::simpleModel::tensorPropNames_;
Foam::PMap<Foam::wordList> Foam::simpleModel::dimScalarPropNames_;
Foam::PMap<Foam::wordList> Foam::simpleModel::dimVectorPropNames_;
Foam::PMap<Foam::wordList> Foam::simpleModel::dimTensorPropNames_;
Foam::PMap<Foam::PList<Foam::word>> Foam::simpleModel::wordProps_;
Foam::PMap<Foam::PList<Foam::fileName>> Foam::simpleModel::fileNameProps_;
Foam::PMap<Foam::PList<Foam::Switch>> Foam::simpleModel::switchProps_;
Foam::PMap<Foam::PList<Foam::scalar>> Foam::simpleModel::scalarProps_;
Foam::PMap<Foam::PList<Foam::vector>> Foam::simpleModel::vectorProps_;
Foam::PMap<Foam::PList<Foam::tensor>> Foam::simpleModel::tensorProps_;
Foam::PMap<Foam::PList<Foam::dimensionedScalar>> Foam::simpleModel::dimScalarProps_;
Foam::PMap<Foam::PList<Foam::dimensionedVector>> Foam::simpleModel::dimVectorProps_;
Foam::PMap<Foam::PList<Foam::dimensionedTensor>> Foam::simpleModel::dimTensorProps_;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleModel::simpleModel
(
    const fvMesh& mesh,
    const word& regionName,
    const word& simpleModelName
):
IOdictionary
(
    IOobject
    (
        regionName+simpleModelName+"Model",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
mesh_(mesh),
regionName_(regionName),
simpleModelName_(simpleModelName),
materialDict_(
    IOobject
    (
        regionName_+"Properties",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    )
),
subDictPath_(mesh.time().constant()+"/"+regionName_+"/"+regionName_+"Properties."+simpleModelName_),
debug_(materialDict_.lookupOrDefault<Switch>("debug", false)),
dynamicMesh_(isA<dynamicFvMesh>(mesh)),
materialDirectory_(materialDict_.isDict("MaterialProperties")?fileName(materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand():""),
constantPropertiesFile_(materialDirectory_+"/constantProperties"),
constantPropertiesDictionary_(
    IOobject
    (
        constantPropertiesFile_,
        mesh.time().db().parent(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
nPatches_(mesh.boundaryMesh().size()),
createFields_(materialDict_.isDict(simpleModelName)?materialDict_.subDict(simpleModelName).lookupOrDefault<List<Tuple2<word,word>>>("createFields",List<Tuple2<word,word>>()):List<Tuple2<word,word>>()),
modelInitialized_("no")
{
  // Read Fields
  tabLevel_++;
  wordList authorizedTypes= {"volScalarField","volVectorField","volTensorField",
                             "surfaceScalarField","surfaceVectorField","surfaceTensorField",
                             "uniformDimensionedScalarField",
                             "uniformDimensionedVectorField",
                             "uniformDimensionedTensorField"
                            };
  forAll(createFields_, fieldI) {
    const word fieldName = createFields_[fieldI].first();
    const word typeName = createFields_[fieldI].second();
    if (!foundInList(typeName, authorizedTypes)) {
      FatalError << "Error in readFields from " << simpleModelName << " Model:"
                 << typeName << " not found in " << authorizedTypes << exit(FatalError);
    }
    if (typeName=="volScalarField") {
      if (mesh_.objectRegistry::foundObject<volScalarField>(fieldName)) {
        createdVolFields_.append(fieldName);
        refVolField<scalar>(fieldName);
      } else {
        createVolField<scalar>(fieldName);
      }
    }
    if (typeName=="volVectorField") {
      if (mesh_.objectRegistry::foundObject<volVectorField>(fieldName)) {
        createdVolFields_.append(fieldName);
        refVolField<vector>(fieldName);
      } else {
        createVolField<vector>(fieldName);
      }
    }
    if (typeName=="volTensorField") {
      if (mesh_.objectRegistry::foundObject<volTensorField>(fieldName)) {
        createdVolFields_.append(fieldName);
        refVolField<tensor>(fieldName);
      } else {
        createVolField<tensor>(fieldName);
      }
    }
    if (typeName=="surfaceScalarField") {
      if (mesh_.objectRegistry::foundObject<surfaceScalarField>(fieldName)) {
        createdSurfaceFields_.append(fieldName);
        refSurfaceField<scalar>(fieldName);
      } else {
        createSurfaceField<scalar>(fieldName);
      }
    }
    if (typeName=="surfaceVectorField") {
      if (mesh_.objectRegistry::foundObject<surfaceVectorField>(fieldName)) {
        createdSurfaceFields_.append(fieldName);
        refSurfaceField<vector>(fieldName);
      } else {
        createSurfaceField<vector>(fieldName);
      }
    }
    if (typeName=="surfaceTensorField") {
      if (mesh_.objectRegistry::foundObject<surfaceTensorField>(fieldName)) {
        createdSurfaceFields_.append(fieldName);
        refSurfaceField<tensor>(fieldName);
      } else {
        createSurfaceField<tensor>(fieldName);
      }
    }
    if (typeName=="uniformDimensionedScalarField") {
      if (mesh_.objectRegistry::foundObject<uniformDimensionedScalarField>(fieldName)) {
        createdUniformFields_.append(fieldName);
        refUniformField<scalar>(fieldName);
      } else {
        createUniformField<scalar>(fieldName);
      }
    }
    if (typeName=="uniformDimensionedVectorField") {
      if (mesh_.objectRegistry::foundObject<uniformDimensionedVectorField>(fieldName)) {
        createdUniformFields_.append(fieldName);
        refUniformField<vector>(fieldName);
      } else {
        createUniformField<vector>(fieldName);
      }
    }
    if (typeName=="uniformDimensionedTensorField") {
      if (mesh_.objectRegistry::foundObject<uniformDimensionedTensorField>(fieldName)) {
        createdUniformFields_.append(fieldName);
        refUniformField<tensor>(fieldName);
      } else {
        createUniformField<tensor>(fieldName);
      }
    }
  }
  tabLevel_--;

  // Add wordList for new regions
  if (!wordPropNames_.found(regionName)) wordPropNames_.insert(regionName, new List<word>());
  if (!fileNamePropNames_.found(regionName)) fileNamePropNames_.insert(regionName, new List<word>());
  if (!switchPropNames_.found(regionName)) switchPropNames_.insert(regionName, new List<word>());
  if (!scalarPropNames_.found(regionName)) scalarPropNames_.insert(regionName, new List<word>());
  if (!vectorPropNames_.found(regionName)) vectorPropNames_.insert(regionName, new List<word>());
  if (!tensorPropNames_.found(regionName)) tensorPropNames_.insert(regionName, new List<word>());
  if (!dimScalarPropNames_.found(regionName)) dimScalarPropNames_.insert(regionName, new List<word>());
  if (!dimVectorPropNames_.found(regionName)) dimVectorPropNames_.insert(regionName, new List<word>());
  if (!dimTensorPropNames_.found(regionName)) dimTensorPropNames_.insert(regionName, new List<word>());

  // Add PropList for new regions
  if (!wordProps_.found(regionName)) wordProps_.insert(regionName, new PList<word>());
  if (!fileNameProps_.found(regionName)) fileNameProps_.insert(regionName, new PList<fileName>());
  if (!switchProps_.found(regionName)) switchProps_.insert(regionName, new PList<Switch>());
  if (!scalarProps_.found(regionName)) scalarProps_.insert(regionName, new PList<scalar>());
  if (!vectorProps_.found(regionName)) vectorProps_.insert(regionName, new PList<vector>());
  if (!tensorProps_.found(regionName)) tensorProps_.insert(regionName, new PList<tensor>());
  if (!dimScalarProps_.found(regionName)) dimScalarProps_.insert(regionName, new PList<dimensionedScalar>());
  if (!dimVectorProps_.found(regionName)) dimVectorProps_.insert(regionName, new PList<dimensionedVector>());
  if (!dimTensorProps_.found(regionName)) dimTensorProps_.insert(regionName, new PList<dimensionedTensor>());

  // Increase the tab level for Info
  tabLevel_++;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleModel::~simpleModel()
{}


// ************************************************************************* //


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

#include "WriteControlIOModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::WriteControlIOModel::WriteControlIOModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleIOModel(mesh, regionName),
outputList_(simpleIOModel::outputList_),
materialPropertiesModel_(refModel<simpleMaterialPropertiesModel>()),
initMatProp_(initMatProp()),
initOutput_(simpleIOModel::initOutput()),
userWriteTimes(materialDict_.subDict("IO").lookupOrDefault<scalarList>("additionalWriteTimes",scalarList())),
userWriteFlags(userWriteTimes.size(), true),
userProbeTimes(materialDict_.subDict("IO").lookupOrDefault<scalarList>("additionalProbingTimes",scalarList())),
userProbeFlags(userProbeTimes.size(), true)
{
  wordList outputList_tmp(materialDict_.subDict("IO").lookupOrDefault<wordList>("writeFields",nullList));
  wordList& outputList = const_cast<wordList&>(outputList_);
  outputList = outputList_tmp;
  foundFieldsInMeshAll(mesh_,outputList);
  modelInitialized();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::WriteControlIOModel::~WriteControlIOModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Switch Foam::WriteControlIOModel::initMatProp()
{
  materialPropertiesModel_.initialize();
  materialPropertiesModel_.update();
  return true;
}

void Foam::WriteControlIOModel::update()
{
//  if (userWriteTimes.size()!=0) {
  forAll(userWriteTimes, timeI) {
    if (mesh_.time().value() >= userWriteTimes[timeI] && userWriteFlags[timeI] ) {
      if (dynamicMesh_) {
        dynamicFvMesh& mesh = (dynamicFvMesh&)(this->mesh_);
        mesh.write();
      } else {
        mesh_.write();
      }
      forAll(outputList_, fieldI) {
        if (mesh_.objectRegistry::foundObject<volScalarField>(outputList_[fieldI])) {
          mesh_.objectRegistry::lookupObject<volScalarField>(outputList_[fieldI]).write();
        }
        if (mesh_.objectRegistry::foundObject<volVectorField>(outputList_[fieldI])) {
          mesh_.objectRegistry::lookupObject<volVectorField>(outputList_[fieldI]).write();
        }
        if (mesh_.objectRegistry::foundObject<volTensorField>(outputList_[fieldI])) {
          mesh_.objectRegistry::lookupObject<volTensorField>(outputList_[fieldI]).write();
        }
        if (mesh_.objectRegistry::foundObject<surfaceScalarField>(outputList_[fieldI])) {
          mesh_.objectRegistry::lookupObject<surfaceScalarField>(outputList_[fieldI]).write();
        }
        if (mesh_.objectRegistry::foundObject<surfaceVectorField>(outputList_[fieldI])) {
          mesh_.objectRegistry::lookupObject<surfaceVectorField>(outputList_[fieldI]).write();
        }
        if (mesh_.objectRegistry::foundObject<surfaceTensorField>(outputList_[fieldI])) {
          mesh_.objectRegistry::lookupObject<surfaceTensorField>(outputList_[fieldI]).write();
        }
      }
      userWriteFlags[timeI]=false;
      break;
    }
  }
//    }
//  if (userProbeTimes.size()!=0) {
  forAll(userProbeTimes, timeI) {
    if (mesh_.time().value() >= userProbeTimes[timeI] && userProbeFlags[timeI]) {
      forAll(probingDictNames_, nameI) {
        sampleFunctions_[nameI].writeOutput();
      }
      userProbeFlags[timeI]=false;
      break;
    }
  }
//    }
  if (mesh_.time().outputTime()) {
    forAll(probingDictNames_, nameI) {
      sampleFunctions_[nameI].writeOutput();
    }
    if(writeFields_) {
      forAll(outputList_,fieldI) {
        if (mesh_.objectRegistry::foundObject<volScalarField>(outputList_[fieldI])) {
          mesh_.objectRegistry::lookupObject<volScalarField>(outputList_[fieldI]).write();
        }
        if (mesh_.objectRegistry::foundObject<volVectorField>(outputList_[fieldI])) {
          mesh_.objectRegistry::lookupObject<volVectorField>(outputList_[fieldI]).write();
        }
        if (mesh_.objectRegistry::foundObject<volTensorField>(outputList_[fieldI])) {
          mesh_.objectRegistry::lookupObject<volTensorField>(outputList_[fieldI]).write();
        }
        if (mesh_.objectRegistry::foundObject<surfaceScalarField>(outputList_[fieldI])) {
          mesh_.objectRegistry::lookupObject<surfaceScalarField>(outputList_[fieldI]).write();
        }
        if (mesh_.objectRegistry::foundObject<surfaceVectorField>(outputList_[fieldI])) {
          mesh_.objectRegistry::lookupObject<surfaceVectorField>(outputList_[fieldI]).write();
        }
        if (mesh_.objectRegistry::foundObject<surfaceTensorField>(outputList_[fieldI])) {
          mesh_.objectRegistry::lookupObject<surfaceTensorField>(outputList_[fieldI]).write();
        }
      }
    }
  }
}

void Foam::WriteControlIOModel::writeOutput()
{
  forAll(probingDictNames_, nameI) {
    sampleFunctions_[nameI].writeOutput();
  }
  forAll(outputList_, fieldI) {
    if (mesh_.objectRegistry::foundObject<volScalarField>(outputList_[fieldI])) {
      mesh_.objectRegistry::lookupObject<volScalarField>(outputList_[fieldI]).write();
    }
    if (mesh_.objectRegistry::foundObject<volVectorField>(outputList_[fieldI])) {
      mesh_.objectRegistry::lookupObject<volVectorField>(outputList_[fieldI]).write();
    }
    if (mesh_.objectRegistry::foundObject<volTensorField>(outputList_[fieldI])) {
      mesh_.objectRegistry::lookupObject<volTensorField>(outputList_[fieldI]).write();
    }
    if (mesh_.objectRegistry::foundObject<volSymmTensorField>(outputList_[fieldI])) {
      mesh_.objectRegistry::lookupObject<volSymmTensorField>(outputList_[fieldI]).write();
    }
    if (mesh_.objectRegistry::foundObject<surfaceScalarField>(outputList_[fieldI])) {
      mesh_.objectRegistry::lookupObject<surfaceScalarField>(outputList_[fieldI]).write();
    }
    if (mesh_.objectRegistry::foundObject<surfaceVectorField>(outputList_[fieldI])) {
      mesh_.objectRegistry::lookupObject<surfaceVectorField>(outputList_[fieldI]).write();
    }
    if (mesh_.objectRegistry::foundObject<surfaceTensorField>(outputList_[fieldI])) {
      mesh_.objectRegistry::lookupObject<surfaceTensorField>(outputList_[fieldI]).write();
    }
  }
}

// ************************************************************************* //

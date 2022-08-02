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
    const word& dictName
)
  :
simpleIOModel(mesh, dictName),
mesh_(mesh),
phaseName_(""),
dictName_(dictName),
outputList_(simpleIOModel::outputList_),
initOutput_(simpleIOModel::initOutput()),
writeControlTimes(materialDict_.subDict("IO").lookup("writeControlTimes")),
writeControlFlags(writeControlTimes.size(), true)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::WriteControlIOModel::~WriteControlIOModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::WriteControlIOModel::update()
{
  forAll(writeControlTimes, timeI) {
    if (mesh_.time().value() > writeControlTimes[timeI] && writeControlFlags[timeI] ) {
      mesh_.write();
      writeOutput();
      writeControlFlags[timeI]=false;
      break;
    }
  }
  if (mesh_.time().outputTime()) {
    writeOutput();
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

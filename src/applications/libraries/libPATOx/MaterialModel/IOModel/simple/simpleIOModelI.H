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

#include "simpleIOModel.H"

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::label& Foam::simpleIOModel::infoDebug() const
{
  return const_cast<label&>(infoDebug_);
}

inline const Foam::List<Foam::fileName>& Foam::simpleIOModel::readFilesList() const
{
  return const_cast<List<fileName>&>(readFilesList_);
}

inline const Foam::List<Foam::List<Foam::scalarList> >& Foam::simpleIOModel::filesData() const
{
  return const_cast<List<List<scalarList> >&>(filesData_);
}

inline bool Foam::simpleIOModel::initOutput()
{

  if (materialDict_.isDict("IO")) {

    forAll(readFilesList_, fileI) {
      filesData_.append(readFileData(readFilesList_[fileI]));
    }

    wordList outputList_tmp(materialDict_.subDict("IO").lookupOrDefault<wordList>("writeFields",nullList));
    wordList probingDictNames_tmp(materialDict_.subDict("IO").lookupOrDefault<wordList>("probingFunctions",nullList));
    outputList_=outputList_tmp;
    probingDictNames_=probingDictNames_tmp;

    forAll(probingDictNames_, nameI) {
      sampleFunctions_.append(new SampleFunction(mesh_, phaseName_, dictName_, probingDictNames_[nameI]));
    }

    foundFieldsInMeshAll(mesh_,outputList_);
  }
  return true;
}


// ************************************************************************* //


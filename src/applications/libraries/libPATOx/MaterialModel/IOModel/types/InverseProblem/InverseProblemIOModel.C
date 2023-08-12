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

#include "InverseProblemIOModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::InverseProblemIOModel::InverseProblemIOModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
noIOModel(mesh, regionName),
startModelInit_(startModelInit()),
runTimeModifiable_(mesh.time().controlDict().lookup("runTimeModifiable")),
writeInterval_(mesh.time().controlDict().lookupOrDefault<scalar>("writeInterval",1e100)),
writeControl_(mesh.time().controlDict().lookupOrDefault<word>("writeControl","timeStep")),
outputList_(noIOModel::outputList_),
readFilesList_(noIOModel::readFilesList_),
filesData_(noIOModel::filesData_),
compareList_(materialDict_.subDict("IO").template lookupOrDefault<List<Tuple2<Tuple2<word,scalar>,Tuple2<word,vector>> >>("compareList",List<Tuple2<Tuple2<word,scalar>,Tuple2<word,vector>> >(0))),
dakotaMethod_(materialDict_.subDict("IO").template lookupOrDefault<word>("dakotaMethod",word::null)),
dakotaResultsFileName_(materialDict_.subDict("IO").template lookupOrDefault<fileName>("dakotaResultsFileName",fileName("results/results").expand())),
indexSamplingCoord_(compareList_.size(),-1),
indexSamplingFields_(compareList_.size(),-1),
outputDeltaTime_(materialDict_.subDict("IO").template lookupOrDefault<scalar>("outputDeltaTime",-1)),
outputTime_(0),
energyModel_(refModel<simpleEnergyModel>()),
Ta_(energyModel_.refVolField<scalar>("Ta")),
initOutput_(noIOModel::initOutput())
{
  dakotaResultsFileName_=changeEnviVar(dakotaResultsFileName_); // change environment variable to full path in the file name
  readInput(); // check the readFiles and fill the data accordingly

  if (compareList_.size() > 0 && dakotaMethod_ == word::null) {
    FatalErrorInFunction << "dakotaMethod not found in \"IO\"" << exit(FatalError);
  }
  // Dakota method check
  if (dakotaMethod_ != word::null) {
    dakotaInit();
  }
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::InverseProblemIOModel::~InverseProblemIOModel()
{
  OFstream os_out(dakotaResultsFileName_);
  forAll(dakotaResults_, resI) {
    os_out << dakotaResults_[resI] << endl;
  }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::InverseProblemIOModel::update()
{
  updateDakotaResults();
  noIOModel::writeOutput();
}

inline void Foam::InverseProblemIOModel::updateDakotaResults()
{
  // Compute dakota results
  if (dakotaMethod_ != word::null) {
    const scalar currentTime_ = mesh_.time().value(); // current time

    // update results if outputFlag = true
    // outputFlag = true if the currentTime is more tham the outputTime
    bool outputFlag = false;
    if (outputDeltaTime_>0) {
      if (currentTime_>=outputTime_) {
        outputFlag = true;
        outputTime_+=outputDeltaTime_;
      }
    }

    forAll(compareList_, compI) {
      const scalarList& times_ = compareReadFileData_[compI][0];
//      if (currentTime_>=min(times_) && currentTime_ <= max(times_) ) {
      if (mesh_.time().outputTime()||outputFlag) {
        compareSampleFunctions_[compI].updateProbing();

        // probingSimulation_[fieldI][dictI][coordI]
        List<List<Field<scalar> > >& probingSimulation_ = compareSampleFunctions_[compI].scalarMasterFields()();

        const label fieldI = indexSamplingFields_[compI];
        const label dictI = 0; // only 1 dict allowed in dakotaInit()
        const label coordI = indexSamplingCoord_[compI];

        // sampling value
        const scalar simulation_ = probingSimulation_[fieldI][dictI][coordI];

        const label columnI = compareReadFileColumn_[compI];

        const scalarList& data_ = compareReadFileData_[compI][columnI];
        // data value from readFiles
        const scalar dataValue_ = linearInterpolation(times_, data_,currentTime_);

        const scalar result_ = sqrt(pow((simulation_-dataValue_),2));

        if (dakotaMethod_ == "MOGA") {
          dakotaResults_[compI]+=result_;
        }

        if (dakotaMethod_ == "nl2sol") {
          dakotaResults_.append(result_);
        }
      }
    }
  }
}

inline void Foam::InverseProblemIOModel::dakotaInit()
{
  wordList availableMethod_;
  availableMethod_.append("nl2sol");
  availableMethod_.append("MOGA");
  Switch correctMethod_ = "no";
  forAll(availableMethod_, methodI) {
    if (dakotaMethod_==availableMethod_[methodI]) {
      correctMethod_ = "yes";
    }
  }
  // Verify dakota method
  if (!correctMethod_) {
    FatalErrorInFunction << dakotaMethod_ << " method not found. The available methods are " << availableMethod_ << exit(FatalError);
  }

  // Verify compare list
  forAll(compareList_, listI) {
    const Tuple2<word,scalar>& compare_ = compareList_[listI].first();
    const Tuple2<word,vector>& read_ = compareList_[listI].second();
    compareReadFile_.append(compare_.first());
    compareFields_.append(read_.first());
    compareReadFileColumn_.append(compare_.second());
    compareFieldsCoord_.append(read_.second());

    // Verify compareReadFile_
    Switch inReadFiles_ = "no";
    forAll(readFilesList_, readI) {
      if (readFilesList_[readI].name() == compareReadFile_[listI] ) {
        inReadFiles_ = "yes";
        compareReadFileData_.append(filesData_[readI]);
      }
    }

    if(!inReadFiles_) {
      FatalErrorInFunction << compareReadFile_[listI] << " from \"compareList\" not found in \"readFiles\"" << exit(FatalError);
    }

    // Verify compareReadFileColumn_
    if (compareReadFileColumn_[listI] >= compareReadFileData_[listI].size()) {
      FatalErrorInFunction << compareReadFile_[listI]  << " column number from \"compareList\" is not correct " << exit(FatalError);
    }
  }

  // Verify field names from compare list are mesh scalarFields
  forAll(compareFields_, wordI) {
    if (!mesh_.objectRegistry::foundObject<volScalarField>(compareFields_[wordI])) {

      FatalErrorInFunction
          << compareFields_[wordI] << " from compareList not found in the mesh scalarFields" << nl
          << "\e[1mvolScalarField available: \e[0m "
          << mesh_.objectRegistry::sortedNames("volScalarField") << nl
          << exit(FatalError);
    }
  }

  forAll(compareList_,compI) {
    // Find in probingFunctions the fields and probes from compareList
    label funcLabel_ = -1;
    forAll(noIOModel::probingDictNames_, dictI) {
      Switch fieldFound_ = "no";
      Switch coordFound_ = "no";

      SampleFunction func_(mesh_, phaseName_, regionName_, noIOModel::probingDictNames_[dictI]);
      PtrList<List<vector> >& probingPoints_(func_.probingPoints());
      PtrList<word>& probingFields_scalar_(func_.probingFields_scalar());

      // FieldFound
      forAll(probingFields_scalar_,fieldI) {

        if (probingFields_scalar_[fieldI] == compareFields_[compI]) {
          fieldFound_="yes";
        }
      }
      // coordFound
      forAll(probingPoints_[0], vectI) {
        if (compareFieldsCoord_[compI] == probingPoints_[0][vectI]) {
          coordFound_ = "yes";
        }
      }

      if (coordFound_ && fieldFound_) {
        funcLabel_=dictI;
        break;
      }
    }

    if(funcLabel_<0) {
      FatalErrorInFunction <<  compareFields_[compI] << " and/or " << compareFieldsCoord_[compI] << " not found in probingFunctions: " << noIOModel::probingDictNames_ << exit(FatalError);
    }

    compareSampleFunctions_.append(new SampleFunction(mesh_, phaseName_, regionName_, noIOModel::probingDictNames_[funcLabel_]));
    PtrList<List<vector> >& probingPoints_(compareSampleFunctions_[compI].probingPoints());
    PtrList<word>& probingFields_scalar_(compareSampleFunctions_[compI].probingFields_scalar());

    if (probingPoints_.size()!=1) {
      FatalErrorInFunction << "Size of \"sets\" != 1 in " << compareSampleFunctions_[compI].dictPath() << exit(FatalError);
    }

    // index of the sampling fields from compareList_ to probingFunctions_
    forAll(probingFields_scalar_,fieldI) {
      if (probingFields_scalar_[fieldI] == compareFields_[compI]) {
        indexSamplingFields_[compI] = fieldI;
      }
    }

    // index of the sampling coord from compareList_ to probingFunctions_
    forAll(probingPoints_[0], vectI) {
      if (compareFieldsCoord_[compI] == probingPoints_[0][vectI]) {
        indexSamplingCoord_[compI] = vectI;
      }
    }
  }

  // Verify method
  if (dakotaMethod_ == "nl2sol") {
    if (compareList_.size() > 1) {
      FatalErrorInFunction << "compareList size has to be 1 with dakota method = nl2sol." << exit(FatalError);
    }
  }
  if (dakotaMethod_ == "MOGA") {
    dakotaResults_.resize(compareList_.size(),0);
  }
}

inline void Foam::InverseProblemIOModel::readInput()
{
  // verify readFiles do not have the same name
  forAll(readFilesList_, fileI) {
    forAll(readFilesList_, fileJ) {
      if (fileI != fileJ && readFilesList_[fileI].name() == readFilesList_[fileJ].name()) {
        FatalErrorInFunction << readFilesList_[fileI] << " has the same file name than " << readFilesList_[fileJ] <<  exit(FatalError);
      }
    }
  }

  forAll(readFilesList_, fileI) {
    word fileName_ = readFilesList_[fileI].name();
    filesData_[fileI];
  }
}

// ************************************************************************* //

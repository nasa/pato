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
    along with OpenFOAM.  If PurePyrolysist, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "PurePyrolysisIOModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PurePyrolysisIOModel::PurePyrolysisIOModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleIOModel(mesh, dictName),
runTimeModifiable_(mesh.time().controlDict().lookup("runTimeModifiable")),
writeInterval_(mesh.time().controlDict().lookupOrDefault<scalar>("writeInterval",1e100)),
writeControl_(mesh.time().controlDict().lookupOrDefault<word>("writeControl","timeStep")),
Ta_(meshLookupOrConstructScalar(mesh, "Ta")),
MaterialChemistryModel_(meshLookupOrConstructModel<simpleMaterialChemistryModel>(mesh_,dictName_,"MaterialChemistry")),
elementNames_(MaterialChemistryModel_.elementNames()),
speciesNames_(MaterialChemistryModel_.speciesNames()),
pyrolysisModel_(meshLookupOrConstructModel<simplePyrolysisModel>(mesh_,dictName_,"Pyrolysis")),
pi_(pyrolysisModel_.pi()),
piTotal_(pyrolysisModel_.piTotal()),
materialPropertiesModel_(meshLookupOrConstructModel<simpleMaterialPropertiesModel>(mesh_,dictName_,"MaterialProperties")),
rho_s_(materialPropertiesModel_.rho_s()),
initial_rho_s_(dimensionedScalar("0", dimMass/pow3(dimLength), scalar(rho_s_[0]))),
outputList_(simpleIOModel::outputList_),
readFilesList_(simpleIOModel::readFilesList_),
filesData_(simpleIOModel::filesData_),
rhoRatioFileName_(simpleIOModel::materialDict_.subDict("IO").template lookupOrDefault<word>("rhoRatioFileName","rhoRatio.txt")),
temperatureFileName_(simpleIOModel::materialDict_.subDict("IO").template lookupOrDefault<word>("temperatureFileName","Temperature.txt")),
normalizedPiTotalFileName_(simpleIOModel::materialDict_.subDict("IO").template lookupOrDefault<word>("normalizedPiTotal","normalizedPiTotal.txt")),
rhoRatio_
(
    IOobject
    (
        "rhoRatio",
        mesh_.time().timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh_,
    dimensionedScalar("0", dimless, scalar(1))
),
normalizedPiTotal_
(
    IOobject
    (
        "normalizedPiTotal",
        mesh_.time().timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", pow(dimTime,-1) ,scalar(0.0))
),
massFractions_(MaterialChemistryModel_.massFractions()),
initNormalizedPi_(initNormalizedPi()),
compareList_(simpleIOModel::materialDict_.subDict("IO").template lookupOrDefault<List<wordList> >("compareList",nullListList_)),
compareFirstTime_(0),
dakotaResultsFileName_(simpleIOModel::materialDict_.subDict("IO").template lookupOrDefault<fileName>("dakotaResultsFileName",fileName("results/results").expand())),
dakotaResponses_(simpleIOModel::materialDict_.subDict("IO").template lookupOrDefault<word>("dakotaResponses",word::null)),
preRunNumber_(simpleIOModel::materialDict_.subDict("IO").template lookupOrDefault<scalar>("preRunNumber", 20)),
initOutput_(simpleIOModel::initOutput())
{
  readInput(); // check the readFiles and fill the data accordingly

  // Dakota method check
  if (dakotaResponses_ != word::null) {
    wordList availableMethod_;
    availableMethod_.append("calibration");
    availableMethod_.append("objective");
    Switch correctMethod_ = "no";
    forAll(availableMethod_, methodI) {
      if (dakotaResponses_==availableMethod_[methodI]) {
        correctMethod_ = "yes";
      }
    }
    // Verify dakota method
    if (!correctMethod_) {
      FatalErrorInFunction << dakotaResponses_ << " method not found. The available methods are " << availableMethod_ << exit(FatalError);
    }
    // Verify system/controlDict
    if (!runTimeModifiable_) {
      FatalErrorInFunction << "runTimeModifiable (system/controlDict) must be \"yes\"" << exit(FatalError);
    }
    if (writeControl_!="timeStep") {
      FatalErrorInFunction << "writeControl (system/controlDict) must be \"timeStep\"" << exit(FatalError);
    }

    // Verify compare list
    if (compareList_.size() == 0) {
      FatalErrorInFunction << "compareList not found." << exit(FatalError);
    }
    forAll(compareList_, coupleI) {
      if (compareList_[coupleI].size() != 2) {
        FatalErrorInFunction << "compareList must be a list of 2 names (data field). \nExample:" << nl
                             << "IO" << nl << "{" << nl << "\tcompareMethod ( (dataFile1 fieldName1) );"
                             << nl << "}"<< exit(FatalError);
      }
      const word& dataFileName_ = compareList_[coupleI][0];
      const word& fieldName_ = compareList_[coupleI][1];
      compareData_.append(dataFileName_);
      compareFields_.append(fieldName_);
    }
    // Verify field names from compare list
    foundFieldsInMesh(mesh_,compareFields_);

    // Verify dataFileName from compare list
    List<int> indexCompareToReadData_(compareData_.size(),-1);
    forAll(compareData_, compareFileI) {
      forAll(readFilesList_, readFileI) {
        if (compareData_[compareFileI] == readFilesList_[readFileI].name()) {
          indexCompareToReadData_[compareFileI] = readFileI;
        }
      }
    }
    forAll(indexCompareToReadData_, compareFileI) {
      if (indexCompareToReadData_[compareFileI] < 0) {
        FatalErrorInFunction << compareData_[compareFileI] << " not found in \"readFiles\"" << exit(FatalError);
      }
    }

    // Verify the times in compare list are the same
    scalarList firstTime_ = filesData_[indexCompareToReadData_[0]][0];
    compareFirstTime_ = min(firstTime_);
    forAll(indexCompareToReadData_, indexI) {
      int index_ = indexCompareToReadData_[indexI];
      scalarList timesData_ = filesData_[index_][0];
      if (firstTime_ != timesData_) {
        FatalErrorInFunction << "Times are different in the compare list" << exit(FatalError);
      }
    }

    // Take the index of the species data (indexSpeciesCompareData_)
    forAll(compareData_, nameI) {
      int indexSpecies_ = -1;
      forAll(speciesNames_, specieI) {
        word fileName_ = compareData_[nameI];
        fileName_.replaceAll(".txt","");
        if (fileName_ == speciesNames_[specieI]) {
          indexSpecies_ = specieI;
        }
      }
      if(indexSpecies_>=0) {
        indexSpeciesCompareData_.append(indexSpecies_);
      }
    }

    // Verify method
    if (dakotaResponses_ == "calibration") {
      if (compareList_.size() > 1) {
        FatalErrorInFunction << "compareList size has to be 1 with dakota method = calibration." << exit(FatalError);
      }
    }
    if (dakotaResponses_ == "objective") {
      dakotaResults_.resize(compareList_.size(),0);
    }

    // Verify first time
    if (compareFirstTime_ <= 0) {
      FatalErrorInFunction << "compareList first time has to be more than 0" << exit(FatalError);
    }
    if ( mesh.time().startTime().value() != 0) {
      FatalErrorInFunction << "start time has to be 0" << exit(FatalError);
    }

    // Update delta time
    const_cast<Time&>(mesh_.time()).setDeltaT(compareFirstTime_/preRunNumber_);

    // Update EndTime
    const_cast<Time&>(mesh_.time()).setEndTime(firstTime_[firstTime_.size()-1]);
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PurePyrolysisIOModel::~PurePyrolysisIOModel()
{
  OFstream os_out(dakotaResultsFileName_);
  forAll(dakotaResults_, resI) {
    os_out << dakotaResults_[resI] << endl;
  }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PurePyrolysisIOModel::update()
{
  updateSampling();
  writeOutput();
}

void Foam::PurePyrolysisIOModel::writeOutput()
{
  if (mesh_.time().outputTime()) {
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
}

void Foam::PurePyrolysisIOModel::updateDeltaTime
(
    const scalarList& times
)
{
  const scalar timeValue = mesh_.time().value();
  scalar deltaTime_=0;
  forAll(times, timeI) {
    if (times[timeI] > timeValue) {
      deltaTime_ = times[timeI] - timeValue;
      if (deltaTime_ > writeInterval_) {
        FatalErrorInFunction
            << "writeInterval (system/controlDict) is too small compared to the times in the data files from ."
            << exit(FatalError);
      }
      const_cast<Time&>(mesh_.time()).setDeltaT(deltaTime_);
      break;
    }
  }
}

void Foam::PurePyrolysisIOModel::updateSampling()
{
  const scalarList& times = temperatureData_[0];
  const scalarList& temps = temperatureData_[1];
  const scalar timeValue = mesh_.time().value();

  if ( timeValue > max(times) || timeValue < min(times)) {
    FatalErrorInFunction << "current time (" << timeValue << " s) is outside of the time table from " << temperatureFileName_  << error(FatalError);
  }
  // update temperature from temperatureFileName_
  forAll(Ta_, cellI) {
    Ta_[cellI]= linearInterpolation(times, temps, timeValue);
  }
  Ta_.correctBoundaryConditions();

  // update rhoRatio, normalizedPiTotal, massFractions and normalizedPi
  rhoRatio_ == rho_s_/initial_rho_s_;
  normalizedPiTotal_ == piTotal_/initial_rho_s_;

  forAll(massFractions_, specieI) {
    forAll(massFractions_[specieI], cellI) {
      scalar piTot = piTotal_[cellI];

      if (piTot != 0) {
        massFractions_[specieI][cellI] = pi_[specieI][cellI]/piTot;
      } else {
        massFractions_[specieI][cellI] = 0;
      }
      normalizedPi_[specieI][cellI] = piTotal_[cellI]*massFractions_[specieI][cellI]/initial_rho_s_.value();
    }
  }

  // Compute dakota results
  if (dakotaResponses_ != word::null) {
    if (timeValue < compareFirstTime_) {
      const_cast<Time&>(mesh_.time()).setDeltaT(compareFirstTime_/preRunNumber_);
    } else {
      Info << "Writing in dakota results at time = " << timeValue << " s" << endl;
      forAll(compareList_, nameI) {
        word fileName_ = compareList_[nameI][0];
        if (normalizedPiTotalFileName_==fileName_) {
          scalar simulation_ = normalizedPiTotal_[0];
          scalar data_ = linearInterpolation(normalizedPiTotalData_[0], normalizedPiTotalData_[1], timeValue);
          if (dakotaResponses_ == "objective") {
            dakotaResults_[nameI]+=sqrt(pow((simulation_-data_)/max(normalizedPiTotalData_[1]),2));
          }
          if (dakotaResponses_ == "calibration") {
            dakotaResults_.append(sqrt(pow((simulation_-data_),2))); // [1/s]
          }
        }

        if (rhoRatioFileName_==fileName_) {
          Info << "compare rhoRatio" << endl;
          scalar simulation_ = rhoRatio_[0];
          scalar data_ = linearInterpolation(rhoRatioData_[0], rhoRatioData_[1], timeValue);
          if (dakotaResponses_ == "objective") {
            dakotaResults_[nameI]+=sqrt(pow((simulation_-data_),2));
          }
          if (dakotaResponses_ == "calibration") {
            dakotaResults_.append(sqrt(pow((simulation_-data_)/max(rhoRatioData_[1]),2))); // [-]
          }
        }
      }

      forAll(indexSpeciesCompareData_, nameI) {
        label specieI = indexSpeciesCompareData_[nameI];
        scalar simulation_ = normalizedPi_[specieI][0];
        scalar data_ = linearInterpolation(speciesData_[specieI][0], speciesData_[specieI][1], timeValue);
        if (dakotaResponses_ == "objective") {
          dakotaResults_[nameI]+=sqrt(pow((simulation_-data_),2));
        }
        if (dakotaResponses_ == "calibration") {
          dakotaResults_.append(sqrt(pow((simulation_-data_),2))); // [1/s]
        }
      }
      // upadte delta time
      updateDeltaTime(times);
    } // timeValue != 0
  }

}

void Foam::PurePyrolysisIOModel::readInput()
{
  // verify readFiles do not have the same name
  forAll(readFilesList_, fileI) {
    forAll(readFilesList_, fileJ) {
      if (fileI != fileJ && readFilesList_[fileI].name() == readFilesList_[fileJ].name()) {
        FatalErrorInFunction << readFilesList_[fileI] << " has the same file name than " << readFilesList_[fileJ] <<  exit(FatalError);
      }
    }
  }

  int indexTemperature=-1;

  forAll(speciesNames_, specieI) {
    speciesData_.append(new List<scalarList>);
  }
  forAll(readFilesList_, fileI) {
    word fileName_ = readFilesList_[fileI].name();
    if (fileName_ == normalizedPiTotalFileName_) {
      normalizedPiTotalData_ = filesData_[fileI];
    }
    if (fileName_ == rhoRatioFileName_) {
      rhoRatioData_ = filesData_[fileI];
    }
    if (fileName_ == temperatureFileName_) {
      temperatureData_ = filesData_[fileI];
      indexTemperature = fileI;
    }
    fileName_.replaceAll(".txt","");
    forAll(speciesNames_, specieI) {
      if (fileName_ == speciesNames_[specieI]) {
        speciesData_[specieI] = filesData_[fileI];
      }
    }
  }

  // Need always temperature data. It will determine the time steps and the temperature of the cell.
  if (temperatureData_.size() == 0) {
    FatalErrorInFunction << temperatureFileName_ << " not found in " << simpleIOModel::materialDict_.path().name() << "/" <<  simpleIOModel::materialDict_.name()
                         << ":" << nl << "IO\n{\n  readFiles\n  {\n   ...\n  }\n}"
                         << exit(FatalError);
  }
  if (temperatureData_.size() < 2) {
    FatalErrorInFunction << "Problem reading " << readFilesList_[indexTemperature]
                         << nl
                         << "The file must have at least two columns."
                         << exit(FatalError);
  }
}

Switch Foam::PurePyrolysisIOModel::initNormalizedPi()
{
  normalizedPi_.resize(massFractions_.size());
  forAll(normalizedPi_, specieI) {
    normalizedPi_.set
    (specieI,
     new volScalarField
     (
         IOobject
         (
             "normalizedPi["+speciesNames_[specieI]+"]",
             mesh_.time().timeName(),
             mesh_,
             IOobject::READ_IF_PRESENT,
             IOobject::NO_WRITE
         ),
         mesh_,
         dimensionedScalar("0", pow(dimTime,-1), scalar(0))
     )
    );
  }
  return true;
}

// ************************************************************************* //

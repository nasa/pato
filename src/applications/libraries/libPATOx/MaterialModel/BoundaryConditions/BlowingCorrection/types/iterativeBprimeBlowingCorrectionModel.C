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

#include "iterativeBprimeBlowingCorrectionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::iterativeBprimeBlowingCorrectionModel::iterativeBprimeBlowingCorrectionModel
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& regionName,
    const label& currentPatchID,
    const dictionary dict
)
  :
simpleBlowingCorrectionModel(mesh, phaseName, regionName, currentPatchID, dict),
energyModel(meshLookupOrConstructModel<simpleEnergyModel>(mesh_,regionName_,simpleEnergyModel::modelName)),
folderPreviousBprime_(changeEnviVar(dict.lookup("folderPreviousBprime")))
{
  // Verify the folder exists
  if (!isDir(folderPreviousBprime_)) {
    FatalErrorInFunction << folderPreviousBprime_ << " not found." << exit(FatalError);
  }

  // Read the previous Bprime values
  Foam::Time runTimePreviousBprime_(folderPreviousBprime_.path(),folderPreviousBprime_.name()); // time in folderPreviousBprime
  Foam::fvMesh meshPreviousBprime_ // mesh in folderPreviousBprime
  (
      Foam::IOobject
      (
          regionName,
          runTimePreviousBprime_.timeName(),
          runTimePreviousBprime_,
          Foam::IOobject::MUST_READ
      )
  );

  // Get the time directories
  fileNameList dirs=readDir(folderPreviousBprime_, fileType::directory);
  forAll(dirs, dirI) {
    if (isNumber(dirs[dirI])) {
      timesPreviousBprime_.append(std::stof(dirs[dirI]));
    }
  }

  // Sorted order of the times
  labelList visitOrder;
  sortedOrder(timesPreviousBprime_, visitOrder);
  timesPreviousBprime_ = scalarList(timesPreviousBprime_, visitOrder);
  valuesPreviousBprime_.resize(meshPreviousBprime_.boundaryMesh()[currentPatchID].size());
  forAll(valuesPreviousBprime_, faceI) {
    valuesPreviousBprime_[faceI].resize(timesPreviousBprime_.size()); // valuesPreviousBprime_[faceI][timeI]
  }

  // Compute the previous Bprime values
  forAll(timesPreviousBprime_, timeI) {
    readOption ro = IOobject::MUST_READ;
    if (timeI==0) {
      ro = IOobject::READ_IF_PRESENT;
    }
    runTimePreviousBprime_.setTime(timesPreviousBprime_[timeI], timeI);
    volScalarField previousBprimeCw
    (
        IOobject
        (
            "BprimeCw",
            runTimePreviousBprime_.timeName(),
            meshPreviousBprime_,
            ro,
            IOobject::NO_WRITE
        ),
        meshPreviousBprime_,
        dimensionedScalar("0", dimless, scalar(0.0))
    );
    volScalarField previousBprimeGw
    (
        IOobject
        (
            "BprimeGw",
            runTimePreviousBprime_.timeName(),
            meshPreviousBprime_,
            ro,
            IOobject::NO_WRITE
        ),
        meshPreviousBprime_,
        dimensionedScalar("0", dimless, scalar(0.0))
    );
    scalarList previousBprime = previousBprimeCw.boundaryField()[currentPatchID_]+previousBprimeGw.boundaryField()[currentPatchID_];
    forAll(previousBprime, faceI) {
      valuesPreviousBprime_[faceI][timeI]=previousBprime[faceI];
    }
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::iterativeBprimeBlowingCorrectionModel::~iterativeBprimeBlowingCorrectionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::iterativeBprimeBlowingCorrectionModel::initialize()
{
  word bold_on="\e[1m";
  word bold_off="\e[0m";
  Info << energyModel.getTabLevel("| ") << bold_on << "Initialize iterativeBprimeBlowingCorrectionModel" <<  bold_off << endl;
  energyModel.tabLevel_++;
  // References to models
  const simpleMassModel& massModel=meshLookupOrConstructModel<simpleMassModel>(mesh_,regionName_,simpleMassModel::modelName);
  // References to fields from Energy Model
  scalarFields_.insert("blowingCorrection",energyModel.refVolField<scalar>("blowingCorrection"));
  scalarFields_.insert("rhoeUeCH",energyModel.refVolField<scalar>("rhoeUeCH"));
  scalarFields_.insert("lambda",energyModel.refVolField<scalar>("lambda"));
  scalarFields_.insert("chemistryOn",energyModel.refVolField<scalar>("chemistryOn"));
  scalarFields_.insert("BprimeCw",energyModel.refVolField<scalar>("BprimeCw"));
  // References to fields from Mass Model
  scalarFields_.insert("mDotGw",massModel.refVolField<scalar>("mDotGw"));
  // Initialized flag
  initialized_=true;
  energyModel.tabLevel_--;
}

void Foam::iterativeBprimeBlowingCorrectionModel::update()
{
  if (!initialized_) {
    initialize();
  }

  scalar min_rhoeUeCH = 1e-6;
  forAll(mesh_.boundaryMesh()[currentPatchID_], faceI) {

    // updated by this function
    scalar& blowingCorrection_BF = scalarFields_["blowingCorrection"].boundaryFieldRef()[currentPatchID_][faceI];
    scalar& rhoeUeCH_BF = scalarFields_["rhoeUeCH"].boundaryFieldRef()[currentPatchID_][faceI];

    // update the previous Bprime value
    const scalar previousBprime_BF = linearInterpolation(timesPreviousBprime_, valuesPreviousBprime_[faceI], mesh_.time().value());

    // constant fields already updated
    const scalar mDotGw_BF = scalarFields_["mDotGw"].boundaryField()[currentPatchID_][faceI];
    const scalar lambda_BF = scalarFields_["lambda"].boundaryField()[currentPatchID_][faceI];
    const scalar chemistryOn_BF = scalarFields_["chemistryOn"].boundaryField()[currentPatchID_][faceI];
    const scalar BprimeCw_BF = scalarFields_["BprimeCw"].boundaryField()[currentPatchID_][faceI];
    const scalar Bprime_BF = BprimeCw_BF + mDotGw_BF / (rhoeUeCH_BF * blowingCorrection_BF);

    if ((int) chemistryOn_BF) {

      // minimum value of rhoeUeCH
      if (rhoeUeCH_BF < min_rhoeUeCH) {
        rhoeUeCH_BF=min_rhoeUeCH;
      }

      // blowing correction condition
      bool condition_ = (rhoeUeCH_BF <= 1e-4);
      if (energyModel.materialDict().isDict("Pyrolysis")) {
        word typePyro = energyModel.materialDict().subDict("Pyrolysis").lookupOrDefault<word>("PyrolysisType", "noPyrolysisSolver<specifiedPyrolysisModel>");
        typePyro.replaceAll("PyrolysisSolver<specifiedPyrolysisModel>", "");
        if (typePyro != "no") {
          condition_ = (mDotGw_BF <= 1e-6) || (rhoeUeCH_BF <= 1e-4);
        }
      }

      if (debug_) {
        Info << energyModel.getTabLevel("| ")  << "--- update blowing correction --- Foam::iterativeBprimeBlowingCorrectionModel::update()" << endl;
        Info << energyModel.getTabLevel("| ")  << "--- rhoeUeCH=" << rhoeUeCH_BF << " --- Foam::iterativeBprimeBlowingCorrectionModel::update()" << endl;
        Info << energyModel.getTabLevel("| ")  << "--- BprimeCw_BF=" << BprimeCw_BF << " --- Foam::iterativeBprimeBlowingCorrectionModel::update()" << endl;
        Info << energyModel.getTabLevel("| ")  << "--- mDotGw_BF=" << mDotGw_BF << " --- Foam::iterativeBprimeBlowingCorrectionModel::update()" << endl;
        Info << energyModel.getTabLevel("| ")  << "--- condition_=" << condition_ << " --- Foam::iterativeBprimeBlowingCorrectionModel::update()" << endl;
      }

      if (((mDotGw_BF / (rhoeUeCH_BF * blowingCorrection_BF) + BprimeCw_BF )==0)||lambda_BF==0) {
        condition_=1;
      }

      if (condition_) {
        blowingCorrection_BF= 1.0;
      } else {
        blowingCorrection_BF = (previousBprime_BF/Bprime_BF)*::log(1+2*lambda_BF*Bprime_BF)/::log(1+2*lambda_BF*previousBprime_BF);
      }
    }
  }
}

void Foam::iterativeBprimeBlowingCorrectionModel::write(Ostream& os) const
{
  os.writeKeyword("blowingCorrectionType") << simpleBlowingCorrectionModel::blowingCorrectionType_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //

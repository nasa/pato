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

#include "constantLambdaBlowingCorrectionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantLambdaBlowingCorrectionModel::constantLambdaBlowingCorrectionModel
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& regionName,
    const label& currentPatchID,
    const dictionary dict
)
  :
simpleBlowingCorrectionModel(mesh, phaseName, regionName, currentPatchID, dict),
energyModel(meshLookupOrConstructModel<simpleEnergyModel>(mesh_,regionName_,simpleEnergyModel::modelName))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantLambdaBlowingCorrectionModel::~constantLambdaBlowingCorrectionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constantLambdaBlowingCorrectionModel::initialize()
{
  word bold_on="\e[1m";
  word bold_off="\e[0m";
  Info << energyModel.getTabLevel("| ") << bold_on << "Initialize constantLambdaBlowingCorrectionModel" <<  bold_off << endl;
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

void Foam::constantLambdaBlowingCorrectionModel::update()
{
  if (!initialized_) {
    initialize();
  }

  scalar min_rhoeUeCH = 1e-6;
  forAll(mesh_.boundaryMesh()[currentPatchID_], faceI) {

    // updated by this function
    scalar& blowingCorrection_BF = scalarFields_["blowingCorrection"].boundaryFieldRef()[currentPatchID_][faceI];
    scalar& rhoeUeCH_BF = scalarFields_["rhoeUeCH"].boundaryFieldRef()[currentPatchID_][faceI];

    // constant fields already updated
    const scalar mDotGw_BF = scalarFields_["mDotGw"].boundaryField()[currentPatchID_][faceI];
    const scalar lambda_BF = scalarFields_["lambda"].boundaryField()[currentPatchID_][faceI];
    const scalar chemistryOn_BF = scalarFields_["chemistryOn"].boundaryField()[currentPatchID_][faceI];
    const scalar BprimeCw_BF = scalarFields_["BprimeCw"].boundaryField()[currentPatchID_][faceI];

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
        Info << energyModel.getTabLevel("| ")  << "--- update blowing correction --- Foam::constantLambdaBlowingCorrectionModel::update()" << endl;
        Info << energyModel.getTabLevel("| ")  << "--- rhoeUeCH=" << rhoeUeCH_BF << " --- Foam::constantLambdaBlowingCorrectionModel::update()" << endl;
        Info << energyModel.getTabLevel("| ")  << "--- BprimeCw_BF=" << BprimeCw_BF << " --- Foam::constantLambdaBlowingCorrectionModel::update()" << endl;
        Info << energyModel.getTabLevel("| ")  << "--- mDotGw_BF=" << mDotGw_BF << " --- Foam::constantLambdaBlowingCorrectionModel::update()" << endl;
        Info << energyModel.getTabLevel("| ")  << "--- condition_=" << condition_ << " --- Foam::constantLambdaBlowingCorrectionModel::update()" << endl;
      }

      if (((mDotGw_BF / (rhoeUeCH_BF * blowingCorrection_BF) + BprimeCw_BF )==0)||lambda_BF==0) {
        condition_=1;
      }

      if (condition_) {
        blowingCorrection_BF= 1.0;
      } else {
        blowingCorrection_BF =
            ::log
            (
                1. + 2. * lambda_BF *
                (
                    mDotGw_BF / (rhoeUeCH_BF * blowingCorrection_BF)
                    + BprimeCw_BF
                )
            ) /
            (
                2. * lambda_BF *
                (
                    mDotGw_BF / (rhoeUeCH_BF * blowingCorrection_BF)
                    + BprimeCw_BF
                )
            );
      }
    }
  }
}

void Foam::constantLambdaBlowingCorrectionModel::write(Ostream& os) const
{
  os.writeKeyword("blowingCorrectionType") << simpleBlowingCorrectionModel::blowingCorrectionType_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //

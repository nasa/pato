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

#include "GradedPorousMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GradedPorousMaterialPropertiesModel::GradedPorousMaterialPropertiesModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
PorousMaterialPropertiesModel(mesh, dictName),
startModelInit_(startModelInit()),
nSolidPhases_(int(createScalarProp("nSolidPhases"))),
pyrolysisModel_(const_cast<simplePyrolysisModel&>(PorousMaterialPropertiesModel::pyrolysisModel_)),
solidEpsI_(pyrolysisModel_.refVolFieldPList<scalar>("epsI_s",1)),
solidEps_(pyrolysisModel_.refVolFieldPList<scalar>("eps_s",1)),
solidRhoI_(pyrolysisModel_.refVolFieldPList<scalar>("rhoI_s",1)),
solidRho_(pyrolysisModel_.refVolFieldPList<scalar>("rho_s",1))
{
  initSolidEps();
  pyrolysisModel_.initialize();
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GradedPorousMaterialPropertiesModel::~GradedPorousMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void Foam::GradedPorousMaterialPropertiesModel::initSolidEps()
{
  if (solidEps_.size() != nSolidPhases_ && solidRho_.size() != nSolidPhases_) {
    ExitError("solidEps_.size() != nSolidPhases_ && solidRho_.size() .size() != nSolidPhases");
  }

  gradingFiles_.resize(nSolidPhases_);
  gradingPatchNames_.resize(nSolidPhases_);
  gradingPatchID_.resize(nSolidPhases_);

  for (int i = 0; i < nSolidPhases_; i++) { // the solid phases are defined from value 1 (0 is attributed to the gas).
    gradingFiles_.set(i,createFileNameProp("gradingFileEpsI["+std::to_string(i+1)+"]","yes",""));
    if (gradingFiles_[i]=="") {
      continue;
    }
    // Name of the patches for the grading distance calculation
    gradingPatchNames_.set(i,createWordProp("gradingPatchNameEpsI["+std::to_string(i+1)+"]"));
    label patchID=mesh_.boundaryMesh().findPatchID(gradingPatchNames_[i]);
    if (patchID<0) {
      ExitError("The grading patch "+gradingPatchNames_[i]+" is not found in mesh.");
    }
    gradingPatchID_.set(i, new label(patchID)); // ID of the grading patches
    List<scalarList> data_(readFileData(gradingFiles_[i])); // distance and phase volume fraction
    if(data_.size()!=2) {
      FatalErrorInFunction << gradingFiles_[i] << " must have two columns (distance and phase volume fraction)" << exit(FatalError);
    }
    if(data_[0].size()!=data_[1].size()) {
      FatalErrorInFunction << gradingFiles_[i] << " must have two columns (distance and phase volume fraction) of the same length" << exit(FatalError);
    }
    List<vector> global_Cf = globalFaceCenters(mesh_, gradingPatchNames_[i], this->regionName_); // global face centers

    forAll(solidEps_[i], cellI) {
      scalar distance_=-1;
      forAll(global_Cf, faceI) {
        scalar d = dist(mesh_.C()[cellI], global_Cf[faceI]);
        if (d<distance_ || distance_<0) {
          distance_=d;
        }
      }
      if ( (distance_>max(data_[0])) || (distance_ < (min(data_[0])))) {
        solidEpsI_[i][cellI] = 0;
      } else {
        solidEpsI_[i][cellI] = linearInterpolation(data_[0],data_[1],distance_);
      }
    }

    forAll(solidEpsI_[i].boundaryField(), patchI) {
      forAll(solidEpsI_[i].boundaryField()[patchI], faceI) {
        scalar distance_=-1;
        if (patchI==gradingPatchID_[i]) {
          distance_=0;
        } else {
          forAll(global_Cf, faceJ) {
            scalar d = dist(mesh_.Cf().boundaryField()[patchI][faceI], global_Cf[faceJ]);
            if (d<distance_ || distance_<0) {
              distance_=d;
            }
          }
        }
        if ( (distance_>max(data_[0])) || (distance_ < (min(data_[0])))) {
          solidEpsI_[i].boundaryFieldRef()[patchI][faceI] = 0;
        } else {
          solidEpsI_[i].boundaryFieldRef()[patchI][faceI] = linearInterpolation(data_[0],data_[1],distance_);
        }
      }
    }
    forAll(solidEps_[i], cellI) {
      solidEps_[i][cellI]=solidEpsI_[i][cellI];
    }
    forAll(solidEps_[i].boundaryField(), patchI) {
      forAll(solidEps_[i].boundaryField()[patchI], faceI) {
        solidEps_[i].boundaryFieldRef()[patchI][faceI]=solidEpsI_[i].boundaryField()[patchI][faceI];
      }
    }
    Info << getTabLevel() << "epsI[" << i+1 << "] graded using " << gradingFiles_[i] <<  nl;
  }
}


// ************************************************************************* //

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

#include "volumeAverageMaterialFailureMassRemovalModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volumeAverageMaterialFailureMassRemovalModel::volumeAverageMaterialFailureMassRemovalModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleMaterialFailureMassRemovalModel(mesh, dictName),
failureCriteriaBC_(failureCriteriaBC()),
solidMechanicsModel_(meshLookupOrConstructModel<simpleSolidMechanicsModel>(mesh_,dictName,simpleSolidMechanicsModel::modelName)),
failureCriteria_(solidMechanicsModel_.createVolFieldIfNotFound<scalar>(solidMechanicsModel_,"failureCriteria", dimensionedScalar("0",dimless, 0),failureCriteriaBC_)),
failureMassLoss_(solidMechanicsModel_.createVolField<scalar>("failureMassLoss", dimensionedScalar("0",dimMass, 0))),
energyModel_(solidMechanicsModel_.refModel<simpleEnergyModel>()),
rho_s_(energyModel_.refVolField<scalar>("rho_s")),
cellMotionU_ptr(dynamicMesh_?&const_cast<volVectorField&>(mesh_.objectRegistry::lookupObject<volVectorField>("cellMotionU")):nullptr),
totalFailureMassLoss_(0),
failureRecessionRateOutputName_(materialDict_.subDict("SolidMechanics").subDict("MaterialFailure").lookupOrDefault<fileName>("failureRecessionRateOutputName","none")),
width_(materialDict_.subDict("SolidMechanics").subDict("MaterialFailure").lookupOrDefault<int>("width",15)),
precision_(materialDict_.subDict("SolidMechanics").subDict("MaterialFailure").lookupOrDefault<int>("precision",3))
{
  if (isDir(failureRecessionRateOutputName_)) {
    FatalError << failureRecessionRateOutputName_ << " exists already." << exit(FatalError);
  } else {
    system("mkdir -p " + failureRecessionRateOutputName_);
  }
  initialTotalMass_=0;
  forAll(rho_s_, cellI) {
    initialTotalMass_ += rho_s_[cellI] * mesh_.V()[cellI];
  }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::volumeAverageMaterialFailureMassRemovalModel::~volumeAverageMaterialFailureMassRemovalModel()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::volumeAverageMaterialFailureMassRemovalModel::updateMassRemoval()
{
  if(dynamicMesh_) {
    cleanFields();
    volVectorField& cellMotionU_ = *cellMotionU_ptr;
    labelList failPatches;
    forAll(failureCriteria_.boundaryField(), patchI) {
      if (failureCriteriaBC_[patchI]=="wedge") continue;
      bool failPatchFound = false;
      forAll(failureCriteria_.boundaryField()[patchI], faceI) {
        if (failureCriteria_.boundaryField()[patchI][faceI]==FAILS) {
          failPatchFound = true;
          scalar distance = distanceFromSurface(patchI, faceI);
          const vector nf =
              - mesh_.Sf().boundaryField()[patchI][faceI]
              / mesh_.magSf().boundaryField()[patchI][faceI]; // normal to the surface
          cellMotionU_.boundaryFieldRef()[patchI][faceI] += nf * distance / mesh_.time().deltaT().value();
          failureMassLoss_.boundaryFieldRef()[patchI][faceI] = rho_s_.boundaryField()[patchI][faceI] * \
              mag(cellMotionU_.boundaryField()[patchI][faceI]) * mesh_.magSf().boundaryField()[patchI][faceI] * \
              mesh_.time().deltaT().value();
          totalFailureMassLoss_+= failureMassLoss_.boundaryFieldRef()[patchI][faceI];
        }
      }
      if (failPatchFound) failPatches.append(patchI);
    }
    Info << "\ttotalFailureMassLoss = "<<totalFailureMassLoss_<<" kg = " << 100*totalFailureMassLoss_/initialTotalMass_ << " % initial mass" << endl;

    if (failureRecessionRateOutputName_ != "none") {
      fileName failureRecessionRateFileName = failureRecessionRateOutputName_ +"/t_"+ mesh_.time().timeName();
      Info << "Writing " << failureRecessionRateFileName << endl;
      os_failureRecessionRate_ptr = new OFstream(failureRecessionRateFileName);
      os_failureRecessionRate_ptr->setf(ios_base::scientific, ios_base::floatfield);
      os_failureRecessionRate_ptr->setf(ios_base::left);
      wordList titles = {"// x[m]","y[m]","z[m]","recessionRate[m/s]"};
      forAll(titles, i) {
        *os_failureRecessionRate_ptr << setw(width_) << titles[i];
      }
      *os_failureRecessionRate_ptr << endl;
      forAll(failPatches, patchI) {
        forAll(mesh_.boundaryMesh()[failPatches[patchI]], faceI) {
          vector xyz = mesh_.Cf().boundaryField()[failPatches[patchI]][faceI];
          scalar recessionRate = mag(cellMotionU_.boundaryField()[failPatches[patchI]][faceI]);
          *os_failureRecessionRate_ptr << setw(width_) << setprecision(precision_) << xyz[0];
          *os_failureRecessionRate_ptr << setw(width_) << setprecision(precision_) << xyz[1];
          *os_failureRecessionRate_ptr << setw(width_) << setprecision(precision_) << xyz[2];
          *os_failureRecessionRate_ptr << setw(width_) << setprecision(precision_) << recessionRate;
          *os_failureRecessionRate_ptr << endl;
        }
      }
    }
  }
}

scalar Foam::volumeAverageMaterialFailureMassRemovalModel::distanceFromSurface(label patchI, label faceI)
{
  if (debug_) {
    Info << "Foam::volumeAverageMaterialFailureMassRemovalModel::distanceFromSurface(" << patchI << ", " << faceI << ")" << endl;
  }
  const scalar delta_x = 0.5/mesh_.deltaCoeffs().boundaryField()[patchI][faceI]; // Distance from face to first cell center
  scalar distance = 0; // m
  int iter_max = 1e4;
  int iter = 0;
  const vector face_center = mesh_.Cf().boundaryField()[patchI][faceI];
  const vector nf =
      - mesh_.Sf().boundaryField()[patchI][faceI]
      / mesh_.magSf().boundaryField()[patchI][faceI]; // normal to the surface
  while (true) {
    iter++;
    if (iter > iter_max) {
      FatalError << "Foam::volumeAverageMaterialFailureMassRemovalModel::distanceFromSurface(" << patchI <<
                 ", " << faceI << "): iter > iter_max" << exit(FatalError);
    }
    distance += delta_x;
    vector probing_point=face_center+nf*distance;
    meshSearchPtr_.reset(new meshSearch(mesh_));
    int cell_id = meshSearchPtr_().findCell(probing_point);
    if (cell_id<0) {
      distance -= delta_x;
      break;
    } else {
      if (failureCriteria_[cell_id]==FAILS) {
        continue;
      } else {
        distance -= delta_x;
        break;
      }
    }
  }
  if (debug_) {
    Info << "Foam::volumeAverageMaterialFailureMassRemovalModel::distanceFromSurface(" << patchI << ", " << faceI \
         << "): Distance =" << distance << " m, No Iterations " << iter << endl;
  }
  return distance;
}

wordList Foam::volumeAverageMaterialFailureMassRemovalModel::failureCriteriaBC()
{
  wordList failureCriteriaBC(mesh_.boundaryMesh().size());
  forAll(mesh_.boundaryMesh(), patchI) {
    if (isA<wedgePolyPatch>(mesh_.boundaryMesh()[patchI])) {
      failureCriteriaBC[patchI]="wedge";
    } else {
      failureCriteriaBC[patchI]="zeroGradient";
    }
  }
  return failureCriteriaBC;
}

void Foam::volumeAverageMaterialFailureMassRemovalModel::cleanFields()
{
  totalFailureMassLoss_ = 0;
  forAll(failureMassLoss_.boundaryField(), patchI) {
    forAll(failureMassLoss_.boundaryField()[patchI], faceI) {
      failureMassLoss_.boundaryFieldRef()[patchI][faceI] = 0;
    }
  }
}

// ************************************************************************* //

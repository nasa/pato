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
    along with OpenFOAM.  If Bprimet, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "erosionModelBoundaryConditions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::erosionModelBoundaryConditions::erosionModelBoundaryConditions
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName,
    const label& currentPatchID,
    const dictionary dict
)
  :
mesh_(mesh),
phaseName_(phaseName),
dictName_(dictName),
currentPatchID_(currentPatchID),
dict_(dict),
dynamicMesh_(isA<dynamicFvMesh>(mesh)),
rho_s_(meshLookupOrConstructScalar(mesh, "rho_s")),
omegaHeterogeneousRate_(meshLookupOrConstructScalar(mesh, "omegaHeterogeneousRate", dimensionedScalar("0", dimMass/pow3(dimLength)/dimTime, 0.0))),
rT(meshLookupOrConstructScalar(mesh, "rT", dimensionedScalar("1", dimLength, 1.0))),
physicsBasedErosionModel_(readScalar(dict_.lookup("physicsBasedErosionModel"))),
recession_(meshLookupOrConstructScalar(mesh, "recession", dimensionedScalar("0", dimLength, scalar(0.0)))),
initialPosition_
(
    meshLookupOrConstructVector(mesh,  "initialPosition",    dimensionedVector("0", dimLength, vector(0.0,0.0,0.0)))
),
materialDict_(
    IOobject
    (
        IOobject::groupName(dictName+"Properties", phaseName),
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    )
),
materialPropertiesDirectory
(
    fileName(materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand()
),
constantPropertiesDictionary
(
    IOobject
    (
        "constantProperties",
        materialPropertiesDirectory,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
rfFail(dimensionedScalar::lookupOrDefault("rfFail",constantPropertiesDictionary,dimensionSet(0, 1, 0, 0, 0, 0, 0),1e-6)),
debug_(materialDict_.lookupOrDefault<Switch>("debug","no"))
{
  if(dynamicMesh_) {
    forAll(initialPosition_.boundaryField()[currentPatchID_], faceI) {
      initialPosition_.boundaryFieldRef()[currentPatchID_][faceI] =  mesh_.Cf().boundaryField()[currentPatchID_][faceI];
    }
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::erosionModelBoundaryConditions::~erosionModelBoundaryConditions()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::erosionModelBoundaryConditions::update()
{
  if(dynamicMesh_) {
    if (debug_) {
      Info << "--- update motion --- Foam::erosionModelBoundaryConditions::update()" << endl;
    }
    updateMotionBC();
  }
}

void Foam::erosionModelBoundaryConditions::updateMotionBC()
{
  // OpenFOAM 4.x
  if(mesh_.objectRegistry::foundObject<volVectorField>("cellMotionU")) {
    volVectorField& cellMotionU_ = const_cast<volVectorField&>(meshLookupOrConstructVector(mesh_, "cellMotionU"));
    // Already updated
    const tmp<vectorField> nf_tmp = - mesh_.Sf().boundaryField()[currentPatchID_]
                                    / mesh_.magSf().boundaryField()[currentPatchID_];
    const vectorField& nf = nf_tmp();
    const tmp<volVectorField> grad_rho_s_tmp = fvc::grad(rho_s_);
    const volVectorField& grad_rho_s = grad_rho_s_tmp();
    const tmp<vectorField> grad_rho_s_int_tmp = grad_rho_s.boundaryField()[currentPatchID_].patchInternalField();
    const vectorField& grad_rho_s_int = grad_rho_s_int_tmp();
    const tmp<scalarField> rho_s_int_tmp = rho_s_.boundaryField()[currentPatchID_].patchInternalField();
    const scalarField& rho_s_int = rho_s_int_tmp();
    const tmp<scalarField> grad_rho_s_int_proj_tmp = grad_rho_s_int & nf;
    const scalarField& grad_rho_s_int_proj = grad_rho_s_int_proj_tmp();
    const tmp<volVectorField> grad_rT_tmp = fvc::grad(rT);
    const volVectorField& grad_rT = grad_rT_tmp();
    const tmp<vectorField> grad_rT_int_tmp = grad_rT.boundaryField()[currentPatchID_].patchInternalField();
    const vectorField& grad_rT_int = grad_rT_int_tmp();
    const tmp<scalarField> grad_rT_int_proj_tmp = grad_rT_int & nf;
    const scalarField& grad_rT_int_proj = grad_rT_int_proj_tmp();
    const tmp<scalarField> rTw_tmp = rT.boundaryField()[currentPatchID_].patchInternalField();
    const scalarField& rTw = rTw_tmp();
    const scalarField& invDxw = mesh_.deltaCoeffs().boundaryField()[currentPatchID_];
    const tmp<scalarField> omegaHeterogeneousRate_int_tmp = omegaHeterogeneousRate_.boundaryField()[currentPatchID_].patchInternalField();
    const scalarField& omegaHeterogeneousRate_int = omegaHeterogeneousRate_int_tmp();

    forAll(cellMotionU_.boundaryFieldRef()[currentPatchID_], faceI) {
      // Already updated
      const vector normal_ = nf[faceI];
      const scalar& rho_s_bf = rho_s_.boundaryField()[currentPatchID_][faceI];

      // Will be updated
      vector&  cellMotionU_BF = cellMotionU_.boundaryFieldRef()[currentPatchID_][faceI];

      if (physicsBasedErosionModel_ <= 0 || physicsBasedErosionModel_ >= 7) {
        FatalErrorInFunction << "Physics based erosion models (\"physicsBasedErosionModel\") are implemented from 1 to 6." << exit(FatalError);
      }
      if (physicsBasedErosionModel_ == 1) { // Rate based on min rf and grad(rf)
        if ((rTw[faceI] <= rfFail.value()) & (grad_rT_int_proj[faceI] > 0)) {
          cellMotionU_BF = (rfFail.value() - rTw[faceI]) / (grad_rT_int_proj[faceI] *  mesh_.time().deltaT().value()) * normal_;
        }
      }

      if (physicsBasedErosionModel_ == 2) { // Rate based on density gradient
        if ((rho_s_bf <= 50.0) & (grad_rho_s_int_proj[faceI] > 0.0)) {
          cellMotionU_BF = (50 - rho_s_bf) / (grad_rho_s_int_proj[faceI] *  mesh_.time().deltaT().value()) * normal_;
        }
      }

      if (physicsBasedErosionModel_ == 3) { // remove cell-by-cell
        if (rTw[faceI]  <= rfFail.value()) {
          cellMotionU_BF=  cellMotionU_BF * 0.98 + 0.02 * 1 / (invDxw[faceI] * mesh_.time().deltaT().value()) * normal_;
        } else {
          cellMotionU_BF = cellMotionU_BF * 0.98;
        }
      }

      if (physicsBasedErosionModel_ == 4) { // Test smoothing
        if ((rTw[faceI]  <= rfFail.value()) ) {
          cellMotionU_BF
            = max(1.05 * cellMotionU_BF & normal_, 1e-7 / rho_s_bf) * normal_;
        } else {
          cellMotionU_BF = cellMotionU_BF * 0.95;
        }
      }

      if (physicsBasedErosionModel_ == 5) { // Test omega
        if ((rTw[faceI]  <= rfFail.value()) ) {
          cellMotionU_BF = 0.95 * cellMotionU_BF + 0.05 * omegaHeterogeneousRate_int[faceI] / rho_s_int[faceI] * normal_;
        } else {
          cellMotionU_BF *= 0.95 ;
        }
      }

      if (physicsBasedErosionModel_ == 6) { // basic

        if ((rTw[faceI]  <= rfFail.value()) ) {
          cellMotionU_BF
            = (1 / invDxw[faceI]) * (1 / mesh_.time().deltaT().value()) * normal_;
        } else {
          cellMotionU_BF = 0.0 * normal_;
        }
      }
    }
  }

  forAll(mesh_.boundaryMesh()[currentPatchID_], faceI) {
    recession_.boundaryFieldRef()[currentPatchID_][faceI] = mag(initialPosition_.boundaryField()[currentPatchID_][faceI] -  mesh_.Cf().boundaryField()[currentPatchID_][faceI]);
  }
}


void Foam::erosionModelBoundaryConditions::write(Ostream& os) const
{
  os.writeKeyword("physicsBasedErosionModel") << physicsBasedErosionModel_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //

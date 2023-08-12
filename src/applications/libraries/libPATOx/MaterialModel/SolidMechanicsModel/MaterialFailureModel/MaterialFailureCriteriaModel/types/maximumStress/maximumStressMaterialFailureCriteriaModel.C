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

#include "maximumStressMaterialFailureCriteriaModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::maximumStressMaterialFailureCriteriaModel::maximumStressMaterialFailureCriteriaModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleMaterialFailureCriteriaModel(mesh, dictName),
failureCriteriaBC_(failureCriteriaBC()),
solidMechanicsModel_(meshLookupOrConstructModel<simpleSolidMechanicsModel>(mesh_,dictName,simpleSolidMechanicsModel::modelName)),
failureCriteria_(solidMechanicsModel_.createVolFieldIfNotFound<scalar>(solidMechanicsModel_,"failureCriteria", dimensionedScalar("0",dimless, 0),failureCriteriaBC_)),
mesh_extend_
(
    lookup_foam_extend_mesh(mesh,dictName)
),
sigma_(Foam_extend_::meshLookupOrConstructSymmTensor(mesh_extend_,"sigma",Foam_extend_::dimensionedSymmTensor("zero", Foam_extend_::dimForce/Foam_extend_::dimArea, Foam_extend_::symmTensor::zero))),
strength_tension_TTT_(readScalar(materialDict_.subDict("SolidMechanics").subDict("MaterialFailure").lookup("strength_tension_TTT"))),
strength_compression_TTT_(readScalar(materialDict_.subDict("SolidMechanics").subDict("MaterialFailure").lookup("strength_compression_TTT"))),
strength_tension_IP_(readScalar(materialDict_.subDict("SolidMechanics").subDict("MaterialFailure").lookup("strength_tension_IP"))),
strength_compression_IP_(readScalar(materialDict_.subDict("SolidMechanics").subDict("MaterialFailure").lookup("strength_compression_IP"))),
strength_shear_(readScalar(materialDict_.subDict("SolidMechanics").subDict("MaterialFailure").lookup("strength_shear")))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::maximumStressMaterialFailureCriteriaModel::~maximumStressMaterialFailureCriteriaModel()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::maximumStressMaterialFailureCriteriaModel::updateFailureCriteria()
{
  // for (const auto& cell : cells) {
  // Use extend max
  // We need to rotate sigma
  // RT*SIGMA*R

  // Algorithm to check if it has failed. Hardcoded
  // All_cells
  List<int> failing_mode(NMODES, 0);
  forAll(sigma_,cellI) {
    failing_mode[TENSION_TTT] = sigma_[cellI].xx() > strength_tension_TTT_ ? FAILS : SUCCEEDS;
    failing_mode[COMPRESSION_TTT] = abs(sigma_[cellI].xx()) > strength_compression_TTT_ ? FAILS : SUCCEEDS;
    failing_mode[TENSION_IP] = sigma_[cellI].yy() > strength_tension_IP_ ? FAILS : SUCCEEDS;
    failing_mode[COMPRESSION_IP] = abs(sigma_[cellI].yy()) > strength_tension_IP_ ? FAILS : SUCCEEDS;
    failing_mode[SHEAR_TTT_TTT] = abs(sigma_[cellI].xy()) > strength_shear_ ? FAILS : SUCCEEDS;
    failing_mode[SHEAR_IP_IP] = abs(sigma_[cellI].xz()) > strength_shear_ ? FAILS : SUCCEEDS;
    failing_mode[SHEAR_TTT_IP] = abs(sigma_[cellI].yz()) > strength_shear_ ? FAILS : SUCCEEDS;

    failureCriteria_[cellI] = max(failing_mode);
  }
  failureCriteria_.correctBoundaryConditions();
}

wordList Foam::maximumStressMaterialFailureCriteriaModel::failureCriteriaBC()
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

void Foam::maximumStressMaterialFailureCriteriaModel::cleanFields()
{
  forAll(failureCriteria_, cellI) {
    failureCriteria_[cellI] = SUCCEEDS;
  }
  forAll(failureCriteria_.boundaryField(), patchI) {
    forAll(failureCriteria_.boundaryField()[patchI], faceI) {
      failureCriteria_.boundaryFieldRef()[patchI][faceI] = SUCCEEDS;
    }
  }
}

// ************************************************************************* //

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

#include "FibrousMaterialTypeAVolumeAblationModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

//namespace Foam
//{
//defineTypeNameAndDebug(FibrousMaterialTypeAVolumeAblationModel, 0);
//}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FibrousMaterialTypeAVolumeAblationModel::FibrousMaterialTypeAVolumeAblationModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleVolumeAblationModel(mesh, dictName),
mesh_(mesh),
dictName_(dictName),
energyConservation_(simpleVolumeAblationModel::materialDict_.subDict("VolumeAblation").template lookupOrDefault<word>("energyConservation","none")),
specificSurface(meshLookupOrConstructScalar(mesh, "specificSurface", dimensionedScalar("0", pow(dimLength,-1), 1.0))),
siteDensity(meshLookupOrConstructScalar(mesh, "siteDensity", dimensionedScalar("0", dimMoles/pow(dimLength,2), 1.0))),
siteDensityIni(meshLookupOrConstructScalar(mesh, "siteDensityIni", dimensionedScalar("0", dimMoles/pow(dimLength,2), 1.0))),
rT(meshLookupOrConstructScalar(mesh, "rT", dimensionedScalar("1", dimLength, 1.0))),
rT0
(
    IOobject
    (
        "rT0",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    rT
),
eps_s_(meshLookupOrConstructScalar(mesh, "eps_s", dimensionedScalar("0", dimless, 0.0), "zeroGradient")),
eps_s0_
(
    IOobject
    (
        "eps_s0",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    eps_s_
),
m_dot_ablation("m_dot_ablation", dimensionSet(1, 0, -1, 0, 0, 0, 0), 0.0),
materialPropertiesDirectory(fileName(simpleVolumeAblationModel::materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand()),
constantPropertiesDictionary
(
    IOobject
    (
        "constantProperties",
        materialPropertiesDirectory,
        mesh.time().db().parent(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
rf0(dimensionedScalar::lookupOrDefault("rf0",constantPropertiesDictionary,dimensionSet(0, 1, 0, 0, 0, 0, 0),5e-6)),
rff(dimensionedScalar::lookupOrDefault("rff",constantPropertiesDictionary,dimensionSet(0, 1, 0, 0, 0, 0, 0),1e-6)),
rfFail(dimensionedScalar::lookupOrDefault("rfFail",constantPropertiesDictionary,dimensionSet(0, 1, 0, 0, 0, 0, 0),1e-6)),
fiberPhase(constantPropertiesDictionary.lookupOrDefault<int>("fiberPhase", 1)-1), // solidEps index starts at 0
matrixPhase(constantPropertiesDictionary.lookupOrDefault<int>("matrixPhase", -1)-1),
materialPropertiesModel_(meshLookupOrConstructModel<simpleMaterialPropertiesModel>(mesh_,dictName_,"MaterialProperties")),
solidRho_(materialPropertiesModel_.solidRho()),
solidEps_(materialPropertiesModel_.solidEps()),
solidEpsI_(materialPropertiesModel_.solidEpsI()),
solidRhoI_(materialPropertiesModel_.solidRhoI()),
rho_ext
(
    IOobject
    (
        "rho_ext",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    solidRho_[fiberPhase]
),
MaterialChemistryModel_(meshLookupOrConstructModel<simpleMaterialChemistryModel>(mesh_,dictName_,"MaterialChemistry")),
omegaHeterogeneousRate_(meshLookupOrConstructScalar(mesh, "omegaHeterogeneousRate", dimensionedScalar("0", dimMass/pow3(dimLength)/dimTime, 0.0))),
gasPropertiesModel_(meshLookupOrConstructModel<simpleGasPropertiesModel>(mesh_,dictName_,"GasProperties")),
rho_g_(gasPropertiesModel_.rho_g()),
M_Cs_("M_Cs", dimensionSet(1, 0, 0, 0, -1, 0, 0), 12.0112)
{
  // cylinder model
  if ((fiberPhase >= 0)&&(matrixPhase >= 0)) {
    eps_s0_ = solidEpsI_[fiberPhase] + solidEpsI_[matrixPhase];
  } else if ((fiberPhase >= 0)) {
    eps_s0_ = solidEpsI_[fiberPhase];
  } else {
    FatalErrorInFunction << "To run the volume ablation model, a microscopic model needs to be defined." << exit(FatalError);
  }

  eps_s_=eps_s0_;

  // initial total radius based on the assumption that the matrix perfectly coats cylindrical fibers
  forAll(rT0, cellI) {
    if ((matrixPhase>=0) && (solidEpsI_[matrixPhase][cellI] > 0)) {
      scalar solidEpsI_fiber = solidEpsI_[fiberPhase][cellI];
      scalar solidEpsI_matrix = solidEpsI_[matrixPhase][cellI];
      rT0[cellI] =
          std::sqrt
          (
              1. +
              solidEpsI_fiber / solidEpsI_matrix
          ) * rf0.value();
    } else {
      rT0[cellI] = rf0.value();
    }
  }
  forAll(rT0.boundaryField(), patchI) {
    forAll(rT0.boundaryField()[patchI], faceI) {
      if ((matrixPhase>=0) && (solidEpsI_[matrixPhase].boundaryField()[patchI][faceI] > 0)) {
        rT0.boundaryFieldRef()[patchI][faceI] =
            Foam::sqrt
            (
                1. +
                solidEpsI_[fiberPhase].boundaryField()[patchI][faceI] / solidEpsI_[matrixPhase].boundaryField()[patchI][faceI]
            ) * rf0.value();
      } else {
        rT0.boundaryFieldRef()[patchI][faceI] = rf0.value();
      }
    }
  }
  rT = rT0; // initialization

  // cylinder model (fiber with or without matrix)
  if ((matrixPhase >= 0)) {
    specificSurface =  rT * 2 * (solidEpsI_[fiberPhase] + solidEpsI_[matrixPhase]) / ( rT0 * rT0);
  } else {
    specificSurface =  rT * 2 * (solidEpsI_[fiberPhase]) / ( rT0 * rT0);
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FibrousMaterialTypeAVolumeAblationModel::~FibrousMaterialTypeAVolumeAblationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FibrousMaterialTypeAVolumeAblationModel::update()
{
  // Model for a type A fibrous material - may be generalized (see Lachaud2010_JSR)
  // update external radius density (matrix or fiber) for matrix-coated fibers being ablated
  forAll (rho_ext, i) {
    if (rT[i] <= rf0.value()) {
      rho_ext[i] = solidRhoI_[fiberPhase][i];
    } else if (matrixPhase >= 0) {
      rho_ext[i] = solidRho_[matrixPhase][i];
    } else {
      FatalErrorInFunction << "Error in the microscopic material model. Wrong initial fiber radius."
                           << "Add a matrix phase or increase radius." << exit(FatalError);
    }
  }

  forAll(rho_ext.boundaryField(), patchI) {
    forAll (rho_ext.boundaryField()[patchI], faceI) {
      if (rT.boundaryField()[patchI][faceI] <= rf0.value()) {
        rho_ext.boundaryFieldRef()[patchI][faceI] = solidRhoI_[fiberPhase].boundaryField()[patchI][faceI];
      } else if (matrixPhase >= 0) {
        rho_ext.boundaryFieldRef()[patchI][faceI] = solidRho_[matrixPhase].boundaryField()[patchI][faceI];
      } else {
        FatalErrorInFunction << "Error in the microscopic material model. Wrong initial fiber radius."
                             << "Add a matrix phase or increase radius." << exit(FatalError);
      }
    }
  }

  // Update solid mass removal by heterogenous reactions
  m_dot_ablation = 0. * m_dot_ablation;
  forAll(omegaHeterogeneousRate_, i) {
    m_dot_ablation.value() +=
        omegaHeterogeneousRate_[i] * mesh_.V()[i]; // total in-depth mass loss rate
  }

  // solid fraction evolution due to volume ablation
  if(this->dynamicMesh_) {
    solve
    (
        fvm::ddt(eps_s_) - fvm::div(mesh_.phi(), eps_s_)  // radius evolution with ALE correction
        + omegaHeterogeneousRate_ / rho_ext // mass loss rate
    );
  } else {
    solve(
        fvm::ddt(eps_s_)  // radius evolution
        + omegaHeterogeneousRate_ / rho_ext // mass loss rate
    );
  }
  eps_s_.max(0);

  // update geometrical properties
  rT = rT0 * sqrt(eps_s_/eps_s0_);
  specificSurface =  rT * 2 * eps_s0_ / (rT0*rT0);
  specificSurface.max(0);
  siteDensity =  siteDensityIni * (1. - exp(-pow4(rT / rff))); // smooth step function centered at rff (rT < rff -> 0, rT > rff -> 1)


  // If we had a material model (conductivity, cp, etc) as a function of phase fractions we would update the fractions
  // For now it is safer to not go into this level of detail when solving the energy conservation

  if (energyConservation_ == "isothermal") {

    forAll(rT, i) {
      if (rT[i] >= rf0.value()) {
        if ((matrixPhase >= 0)) {
          solidEps_[matrixPhase][i] =
              solidEpsI_[matrixPhase][i] *
              (rT[i] * rT[i] - rf0.value() * rf0.value()) /
              (rT0[i] * rT0[i]  - rf0.value() * rf0.value());
        }
        solidEps_[fiberPhase][i] = solidEpsI_[fiberPhase][i];
      } else {
        if ((matrixPhase >= 0)) {
          solidEps_[matrixPhase][i] = 0;
        }
        solidEps_[fiberPhase][i] =
            solidEpsI_[fiberPhase][i] *
            (rT[i] * rT[i]) /
            (rf0.value() * rf0.value());
      }
    }

    forAll(mesh_.boundaryMesh(), patchI) {
      const tmp<scalarField> rTint_tmp = rT.boundaryField()[patchI].patchInternalField();
      const scalarField& rTint = rTint_tmp();
      const tmp<scalarField> rT0int_tmp = rT0.boundaryField()[patchI].patchInternalField();
      const scalarField& rT0int = rT0int_tmp();

      forAll(rTint, i) {
        if (rTint[i] >= rf0.value()) {
          if ((matrixPhase >= 0)) {
            solidEps_[matrixPhase].boundaryFieldRef()[patchI][i] =
                solidEpsI_[matrixPhase].boundaryField()[patchI][i] *
                (rTint[i] * rTint[i] - rf0.value() * rf0.value()) /
                (rT0int[i] * rT0int[i]  - rf0.value() * rf0.value());
          }
          solidEps_[fiberPhase].boundaryFieldRef()[patchI][i] = solidEpsI_[fiberPhase].boundaryField()[patchI][i];
        } else {
          if ((matrixPhase >= 0)) {
            solidEps_[matrixPhase].boundaryFieldRef()[patchI][i] = 0. ;
          }
          solidEps_[fiberPhase].boundaryFieldRef()[patchI][i] =
              solidEpsI_[fiberPhase].boundaryField()[patchI][i] *
              (rTint[i] * rTint[i]) /
              (rf0.value() * rf0.value()) ;
        }
      }
    }
  } // end "isothermal"
}

void Foam::FibrousMaterialTypeAVolumeAblationModel::updateSolidCarbonMassFraction(volScalarField& YiCs)
{
  // loads apparent mass fraction of porous C(gr)
  // is divided by "M_Cs / rho_g" in MaterialChemistry solver to get the concentration
  YiCs = specificSurface * siteDensity * M_Cs_ / rho_g_ ;
}


// ************************************************************************* //

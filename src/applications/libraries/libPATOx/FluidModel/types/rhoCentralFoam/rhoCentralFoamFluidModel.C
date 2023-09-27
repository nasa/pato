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

#include "rhoCentralFoamFluidModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoCentralFoamFluidModel::rhoCentralFoamFluidModel
(
    Time& runTime
)
  :
basicFluidModel(runTime),
regionList(rp["fluid"]),
meshPtr
(
    dynamicFvMesh::New
    (
        IOobject
        (
            regionList[0], // only first fluid region
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    )
),
mesh(meshPtr()),
pThermo
(
    psiThermo::New(mesh)
),
thermo(pThermo()),
e(thermo.he()),
U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
),
rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
),
rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
),
rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
),
pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
),
neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
),
phi("phi", fvc::flux(rhoU)),
turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
),
p(thermo.p()),
T(thermo.T()),
psi(thermo.psi()),
mu(thermo.mu()),
v_zero("v_zero", dimVolume/dimTime, 0.0),
rho_pos(interpolate(rho, pos)),
rho_neg(interpolate(rho, neg)),
rhoU_pos(interpolate(rhoU, pos, U.name())),
rhoU_neg(interpolate(rhoU, neg, U.name())),
U_pos("U_pos", rhoU_pos/rho_pos), // ERROR: rho_pos == 0
U_neg("U_neg", rhoU_neg/rho_neg),
aphiv_pos("aphiv_pos",  U_pos & mesh.Sf()),
aphiv_neg("aphiv_neg", U_neg&mesh.Sf()),
phiv_pos("phiv_pos", U_pos & mesh.Sf()),
phiv_neg("phiv_neg", U_neg & mesh.Sf()),
rPsi("rPsi", 1.0/psi),
rPsi_pos(interpolate(rPsi, pos, T.name())),
rPsi_neg(interpolate(rPsi, neg, T.name())),
c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi)),
cSf_pos
(
    "cSf_pos",
    interpolate(c, pos, T.name())*mesh.magSf()
),
cSf_neg
(
    "cSf_neg",
    interpolate(c, neg, T.name())*mesh.magSf()
),
ap
(
    "ap",
    max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
),
am
(
    "am",
    min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
),
a_pos("a_pos", ap/(ap - am)),
a_neg("a_neg", 1.0-a_pos),
p_pos("p_pos", rho_pos*rPsi_pos),
p_neg("p_neg", rho_neg*rPsi_neg),
e_pos(interpolate(e, pos, T.name())),
e_neg(interpolate(e, neg, T.name())),
aSf("aSf", am*a_pos),
stitchingFluidPatchName
(
    const_cast<Time&>(mesh.time()).controlDict().lookup("stitchingFluidPatch")
),
maxDeltaTw_(basicFluidModel::maxDeltaTw_)

{
#include "rhoCentralFoam/createStitchFields.H"
#include "rhoCentralFoam/createFields.H"
#include "rhoCentralFoam/createFieldRefs.H"
#include "rhoCentralFoam/createTimeControls.H"
#include "rhoCentralFoam/createRDeltaT.H"

  turbulence->validate();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rhoCentralFoam/readFluxScheme.H"

// Courant numbers used to adjust the time-step
  CoNum = 0.0;
  meanCoNum = 0.0;

  Info<< "\nStarting time loop\n" << endl;

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rhoCentralFoamFluidModel::~rhoCentralFoamFluidModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::rhoCentralFoamFluidModel::updateBefore()
{
  // --- Directed interpolation of primitive fields onto faces

  rho_pos=(interpolate(rho, pos));
  rho_neg=(interpolate(rho, neg));

  rhoU_pos=(interpolate(rhoU, pos, U.name()));
  rhoU_neg=(interpolate(rhoU, neg, U.name()));

  rPsi=1.0/psi;
  rPsi_pos=interpolate(rPsi, pos, T.name());
  rPsi_neg=interpolate(rPsi, neg, T.name());

  e_pos=(interpolate(e, pos, T.name()));
  e_neg=(interpolate(e, neg, T.name()));

  U_pos=rhoU_pos/rho_pos;
  U_neg=rhoU_neg/rho_neg;

  p_pos=rho_pos*rPsi_pos;
  p_neg=rho_neg*rPsi_neg;

  phiv_pos=U_pos & mesh.Sf();
  phiv_neg=U_neg & mesh.Sf();

  c=sqrt(thermo.Cp()/thermo.Cv()*rPsi);
  cSf_pos= interpolate(c, pos, T.name())*mesh.magSf();
  cSf_neg= interpolate(c, neg, T.name())*mesh.magSf();

  ap=max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);
  am= min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);

  a_pos=ap/(ap - am);

  surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

  aSf=am*a_pos;

  if (fluxScheme == "Tadmor") {
    aSf = -0.5*amaxSf;
    a_pos = 0.5;
  }

  a_neg=1.0 - a_pos;

  phiv_pos *= a_pos;
  phiv_neg *= a_neg;

  aphiv_pos=phiv_pos - aSf;
  aphiv_neg=phiv_neg + aSf;

  // Reuse amaxSf for the maximum positive and negative fluxes
  // estimated by the central scheme
  amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

#include "rhoCentralFoam/centralCourantNo.H"
#include "readTimeControls.H"

  if (LTS) {
#include "rhoCentralFoam/setRDeltaT.H"
  } else {
#include "rhoCentralFoam/setDeltaT.H"
  }

}

void Foam::rhoCentralFoamFluidModel::updateAfter()
{

  // Do any mesh changes
  mesh.update();

  phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

  surfaceVectorField phiUp
  (
      (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
      + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf()
  );

  surfaceScalarField phiEp
  (
      "phiEp",
      aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
      + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
      + aSf*p_pos - aSf*p_neg
  );

  volScalarField muEff("muEff", turbulence->muEff());
  volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

  // --- Solve density
  solve(fvm::ddt(rho) + fvc::div(phi));

  // --- Solve momentum
  solve(fvm::ddt(rhoU) + fvc::div(phiUp));

  U.ref() =
      rhoU()
      /rho();
  U.correctBoundaryConditions();
  rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

  if (!inviscid) {
    solve
    (
        fvm::ddt(rho, U) - fvc::ddt(rho, U)
        - fvm::laplacian(muEff, U)
        - fvc::div(tauMC)
    );
    rhoU = rho*U;
  }

  // --- Solve energy
  surfaceScalarField sigmaDotU
  (
      "sigmaDotU",
      (
          fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
          + fvc::dotInterpolate(mesh.Sf(), tauMC)
      )
      & (a_pos*U_pos + a_neg*U_neg)
  );

  solve
  (
      fvm::ddt(rhoE)
      + fvc::div(phiEp)
      - fvc::div(sigmaDotU)
  );

  e = rhoE/rho - 0.5*magSqr(U);
  e.correctBoundaryConditions();
  thermo.correct();
  rhoE.boundaryFieldRef() ==
  rho.boundaryField()*
  (
      e.boundaryField() + 0.5*magSqr(U.boundaryField())
  );

  if (!inviscid) {
    solve
    (
        fvm::ddt(rho, e) - fvc::ddt(rho, e)
        - fvm::laplacian(turbulence->alphaEff(), e)
    );
    thermo.correct();
    rhoE = rho*(e + 0.5*magSqr(U));
  }

  p.ref() =
      rho()
      /psi();
  p.correctBoundaryConditions();
  rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();

  turbulence->correct();

  // Update multi region coupled criteria

#include "rhoCentralFoam/fluidStitchTolerance.H"

}

// ************************************************************************* //

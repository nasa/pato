bool closedVolume = p_rgh.needReference();
dimensionedScalar compressibility = fvc::domainIntegrate(psi);
bool compressible = (compressibility.value() > SMALL);

rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

surfaceScalarField phig("phig", -rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(rho*HbyA)
        + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
    )
    + phig
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
  fvScalarMatrix p_rghEqn
  (
      fvm::ddt(psi, p_rgh)
      + fvc::ddt(psi, rho)*gh
      + fvc::ddt(psi)*pRef
      + fvc::div(phiHbyA)
      - fvm::laplacian(rhorAUf, p_rgh)
      ==
      parcels.Srho()
      + surfaceFilm.Srho()
      + fvOptions(psi, p_rgh, rho.name())
  );

  p_rghEqn.solve();

  if (nonOrth == nNonOrthCorr) {
    phi = phiHbyA + p_rghEqn.flux();
    U = HbyA + rAU*fvc::reconstruct((p_rghEqn.flux() + phig)/rhorAUf);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
  }
}

p = p_rgh + rho*gh + pRef;

// Solve continuity
#include "fireFoam_rhoEqn.H"

// Update continuity errors
#include "compressibleContinuityErrors.H"

K = 0.5*magSqr(U);

// Update pressure time derivative if needed
if (thermo.dpdt())
{
  dpdt = fvc::ddt(p);
}

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume && compressible)
{
  p += (initialMass - fvc::domainIntegrate(thermo.rho()))
       /compressibility;
  rho = thermo.rho();
  p_rgh = p - rho*gh;
}

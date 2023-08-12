{
  label pRefCell = 0;
  scalar pRefValue = 0.0;
  setRefCell(p, simple, pRefCell, pRefValue);
  mesh.setFluxRequired(p.name());

  volScalarField rAU(1.0/UEqn.A());
  volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
  surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
  MRF.makeRelative(phiHbyA);
  adjustPhi(phiHbyA, U, p);

  tmp<volScalarField> rAtU(rAU);

  if (consistent) {
    rAtU = 1.0/(1.0/rAU - UEqn.H1());
    phiHbyA +=
    fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
  }

  tUEqn.clear();

  // Update the pressure BCs to ensure flux consistency
  constrainPressure(p, U, phiHbyA, rAtU(), MRF);

  // Non-orthogonal pressure corrector loop
  for (int oCorr=0; oCorr<=nNonOrthCorr; oCorr++) {
    bool finalIter = oCorr == nNonOrthCorr;

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    if (finalIter) {
      phi = phiHbyA - pEqn.flux();
    }
  }

#include "continuityErrors.H"

  // Explicitly relax pressure for momentum corrector
  p.relax();

  // Momentum corrector
  U = HbyA - rAtU()*fvc::grad(p);
  U.correctBoundaryConditions();
  fvOptions.correct(U);
}

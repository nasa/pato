/* Imported from OpenFOAM-7 CorrectPhi.C.
It needs pimple control to work. because PATO cannot use pimple controls we need
    to rewrite the function. */

correctUphiBCs(U, phi, true);

wordList pcorrTypes
(
    p.boundaryField().size(),
    zeroGradientFvPatchScalarField::typeName
);

// Set BCs of pcorr to fixed-value for patches at which p is fixed
forAll(p.boundaryField(), patchi)
{
  if (p.boundaryField()[patchi].fixesValue()) {
    pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
  }
}

volScalarField pcorr
(
    IOobject
    (
        "pcorr",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions(), 0),
    pcorrTypes
);

if (pcorr.needReference())
{
  fvc::makeRelative(phi, U);
  adjustPhi(phi, U, pcorr);
  fvc::makeAbsolute(phi, U);
}

mesh.setFluxRequired(pcorr.name());

for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
{
  bool finalIter = nonOrth == nNonOrthCorr;

  // Solve for pcorr such that the divergence of the corrected flux
  // matches the divU provided (from previous iteration, time-step...)
  fvScalarMatrix pcorrEqn
  (
      fvm::laplacian(dimensionedScalar("rAUf", dimTime, 1), pcorr)
      ==
      fvc::div(phi) - geometricZeroField()
  );

  pcorrEqn.setReference(0, 0);

  pcorrEqn.solve();

  if (finalIter) {
    phi -= pcorrEqn.flux();
  }
}

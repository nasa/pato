// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::div(phi, U)
    + MRF.DDt(U)
    + turb.divDevReff(U)
    ==
    fvOptions(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (momentumPredictor)
{
  solve(UEqn == -fvc::grad(p));

  fvOptions.correct(U);
}

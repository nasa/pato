if (divSigmaExpMethod == "standard")
{
  divSigmaExp = Foam_extend_::fvc::div
                (
                    mu*gradU.T() + lambda*(I*tr(gradU)) - (mu + lambda)*gradU,
                    "div(sigma)"
                );
} else if(divSigmaExpMethod == "surface")
{
  divSigmaExp = Foam_extend_::fvc::div
                (
                    muf*(mesh.Sf() & Foam_extend_::fvc::interpolate(gradU.T()))
                    + lambdaf*(mesh.Sf() & I*Foam_extend_::fvc::interpolate(tr(gradU)))
                    - (muf + lambdaf)*(mesh.Sf() & Foam_extend_::fvc::interpolate(gradU))
                );
} else if(divSigmaExpMethod == "decompose")
{
  snGradU = Foam_extend_::fvc::snGrad(U);

  surfaceTensorField shearGradU(((I - n*n) & Foam_extend_::fvc::interpolate(gradU)));

  divSigmaExp = Foam_extend_::fvc::div
                (
                    mesh.magSf()*
                    (
                        - (muf + lambdaf)*(snGradU & (I - n*n))
                        + lambdaf*tr(shearGradU & (I - n*n))*n
                        + muf*(shearGradU & n)
                    )
                );
} else if(divSigmaExpMethod == "expLaplacian")
{
  divSigmaExp =
      - Foam_extend_::fvc::laplacian(mu + lambda, U, "laplacian(DU,U)")
      + Foam_extend_::fvc::div(mu*gradU.T() + lambda*(I*tr(gradU)), "div(sigma)");
} else
{
  FatalError
      << "divSigmaExp method " << divSigmaExpMethod << " not found!" << exit(FatalError);
}

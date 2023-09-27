if(divSigmaExpMethod == "standard")
{
  //- calculating the full gradient has good convergence and no high freq oscillations
  divSigmaExp = Foam_extend_::fvc::div((C && Foam_extend_::symm(gradU)) - (K & gradU), "div(sigma)");
} else if(divSigmaExpMethod == "surface")
{
  //- this form seems to have the best convergence
  divSigmaExp = Foam_extend_::fvc::div
                (
                    mesh.magSf()*
                    (
                        (n & (Cf && Foam_extend_::fvc::interpolate(symm(gradU))))
                        - (n & (Kf & Foam_extend_::fvc::interpolate(gradU)))
                    )
                );
} else if(divSigmaExpMethod == "laplacian")
{
  //- can cause high freq oscillations and slow convergence
  divSigmaExp =
      Foam_extend_::fvc::div(C && symm(epsilon), "div(sigma)")
      - Foam_extend_::fvc::laplacian(K, U, "laplacian(K, U)");
} else
{
  FatalError
      << "divSigmaExp method " << divSigmaExpMethod << " not found!" << exit(FatalError);
}

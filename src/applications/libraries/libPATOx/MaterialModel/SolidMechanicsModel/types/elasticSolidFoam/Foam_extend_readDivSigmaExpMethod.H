//- how explicit component of sigma is to be calculated
word divSigmaExpMethod_tmp(mesh.solutionDict().subDict("solidMechanics").lookup("divSigmaExp"));
divSigmaExpMethod=divSigmaExpMethod_tmp;
Info << "Selecting divSigmaExp calculation method " << divSigmaExpMethod <<  endl;
if(divSigmaExpMethod != "standard" && divSigmaExpMethod != "surface" && divSigmaExpMethod != "decompose" && divSigmaExpMethod != "laplacian")
{
  FatalError << "divSigmaExp method " << divSigmaExpMethod << " not found!" << nl
             << "valid methods are:\nstandard\nsurface\ndecompose\nlaplacian"
             << exit(FatalError);
}

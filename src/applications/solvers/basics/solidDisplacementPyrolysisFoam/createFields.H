#include "readMechanicalProperties.H"
#include "readThermalProperties.H"

Info<< "Reading field D\n" << endl;
volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


autoPtr<volScalarField> Tptr(NULL);

if (thermalStress)
{
  Info<< "Reading field T\n" << endl;
  Tptr.reset
  (
      new volScalarField
      (
          IOobject
          (
              "T",
              runTime.timeName(),
              mesh,
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE
          ),
          mesh
      )
  );
}


Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mu*twoSymm(fvc::grad(D)) + lambda*(I*tr(fvc::grad(D)))
);

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);

if (compactNormalStress)
{
  divSigmaExp -= fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)");
} else
{
  divSigmaExp -= fvc::div((2*mu + lambda)*fvc::grad(D), "div(sigmaD)");
}

mesh.setFluxRequired(D.name());



//*** Pyrolysis ***//

// Opening the material property file ...

IOdictionary pyrolysisProperties
(
    IOobject
    (
        "pyrolysisProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar R(constant::physicoChemical::R); // ideal gas constant
dimensionedScalar A(pyrolysisProperties.lookup("A"));
dimensionedScalar Ep(pyrolysisProperties.lookup("Ep"));

Switch pyrolysisStress(pyrolysisProperties.lookup("pyrolysisStress"));
Switch dynamicMesh(mechanicalProperties.lookup("dynamicMesh"));


volScalarField Chi
(
    IOobject
    (
        "Chi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

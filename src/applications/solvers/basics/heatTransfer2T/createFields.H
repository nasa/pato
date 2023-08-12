/* Reading initial conditions */
Info<< "Reading field T_s\n" << endl;

volScalarField T_s
(
    IOobject
    (
        "T_s",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field T_g\n" << endl;

volScalarField T_g
(
    IOobject
    (
        "T_g",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

/* Reading material properties */
Info<< "Reading solidProperties\n" << endl;

IOdictionary solidProperties
(
    IOobject
    (
        "solidProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar eps
(
    solidProperties.lookup("eps")
);

dimensionedScalar rho_s0
(
    solidProperties.lookup("rho_s")
);

dimensionedScalar rho_g0
(
    solidProperties.lookup("rho_g")
);

volScalarField rho_s
(
    IOobject
    (
        "rho_s",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rho_s0",dimensionSet(1,-3,0,0,0,0,0),rho_s0.value())
);

volScalarField rho_g
(
    IOobject
    (
        "rho_g",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rho_g0",dimensionSet(1,-3,0,0,0,0,0),rho_g0.value())
);

dimensionedScalar cp_s
(
    solidProperties.lookup("cp_s")
);

dimensionedScalar cp_g
(
    solidProperties.lookup("cp_g")
);

dimensionedTensor k_s
(
    solidProperties.lookup("k_s")
);

dimensionedTensor k_g
(
    solidProperties.lookup("k_g")
);

dimensionedVector u_g
(
    solidProperties.lookup("u_g")
);

dimensionedScalar hgs
(
    solidProperties.lookup("hgs")
);

dimensionedScalar hsg
(
    solidProperties.lookup("hsg")
);

surfaceScalarField phi ("phi", u_g & mesh.Sf());

#include "createFvOptions.H"

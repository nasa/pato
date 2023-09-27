/* Reading initial conditions */
Info<< "Reading field p\n" << endl;

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

/* Reading physical properties */
Info<< "Reading physical properties\n" << endl;

IOdictionary constantProperties
(
    IOobject
    (
        "constantProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar mu
(
    constantProperties.lookup("mu")
);

dimensionedScalar M
(
    constantProperties.lookup("M")
);

dimensionedScalar R
(
    constantProperties.lookup("R")
);

dimensionedTensor K
(
    constantProperties.lookup("K")
);

dimensionedScalar eps
(
    constantProperties.lookup("eps")
);

dimensionedScalar rho_s0
(
    constantProperties.lookup("rho_s")
);

dimensionedScalar cp_s
(
    constantProperties.lookup("cp_s")
);

dimensionedScalar cp_g
(
    constantProperties.lookup("cp_g")
);

dimensionedTensor k_s
(
    constantProperties.lookup("k_s")
);

dimensionedTensor k_g
(
    constantProperties.lookup("k_g")
);

dimensionedScalar hgs
(
    constantProperties.lookup("hgs")
);

dimensionedScalar hsg
(
    constantProperties.lookup("hsg")
);

volTensorField Gamma
(
    IOobject
    (
        "Gamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    p*M/(mu*R*T_g)*K
);

volScalarField Beta
(
    IOobject
    (
        "Beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    eps*M/(R*T_g)
);

volVectorField u_g
(
    IOobject
    (
        "u_g",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("u_g", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero)
);
u_g.write();

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
    p*M/(R*T_g)
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

volVectorField q
(
    IOobject
    (
        "q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    u_g*p*M/(R*T_g)
);
q.write();

surfaceScalarField phi_g ("phi_g", linearInterpolate(q) & mesh.Sf());

#include "createFvOptions.H"

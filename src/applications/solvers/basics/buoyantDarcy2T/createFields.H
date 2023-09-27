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

Info<< "Reading field Ta\n" << endl;

volScalarField Ta
(
    IOobject
    (
        "Ta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Tg\n" << endl;

volScalarField Tg
(
    IOobject
    (
        "Tg",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "\nReading g" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
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

dimensionedScalar eps_g
(
    constantProperties.lookup("eps_g")
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

dimensionedTensor k_a
(
    constantProperties.lookup("k_a")
);

dimensionedTensor k_g
(
    constantProperties.lookup("k_g")
);

dimensionedScalar Hv0
(
    constantProperties.lookup("Hv0")
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
    p*M/(R*Tg)
);

rho_g.write();

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
    rho_g/mu*K
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
    eps_g*M/(R*Tg)
);

volTensorField Theta
(
    IOobject
    (
        "Theta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Gamma*M/(R*Tg)
);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero)
);
U.write();

volVectorField vG
(
    IOobject
    (
        "vG",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U/eps_g
);
vG.write();

volVectorField mDotG
(
    IOobject
    (
        "mDotG",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho_g*U
);

surfaceScalarField phi_g
(
    IOobject
    (
        "phi_g",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(mDotG) & mesh.Sf()
);

surfaceScalarField phiG("phiG", linearInterpolate((Theta & g)) & mesh.Sf());

#include "createFvOptions.H"

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

    dimensionedScalar rho_s
    (
        solidProperties.lookup("rho_s")
    );

    dimensionedScalar rho_g
    (
        solidProperties.lookup("rho_g")
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

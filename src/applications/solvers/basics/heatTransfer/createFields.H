/* Reading initial conditions */
Info<< "Reading field T\n" << endl;

volScalarField T
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

Info<< "Reading density \n" << endl;

dimensionedScalar rho
(
    solidProperties.lookup("rho")
);


dimensionedScalar cp
(
    solidProperties.lookup("cp")
);

dimensionedTensor k
(
    solidProperties.lookup("k")
);


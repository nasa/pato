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

dimensionedScalar eps
(
    constantProperties.lookup("eps")
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

dimensionedScalar T
(
    constantProperties.lookup("T")
);

dimensionedTensor K
(
    constantProperties.lookup("K")
);

dimensionedTensor Beta
(
    constantProperties.lookup("Beta")
);

volVectorField v
(
    IOobject
    (
        "v",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("v", dimensionSet(0, 1, -1, 0, 0, 0, 0), vector::zero)
);

v.write();


volScalarField magV
(
    IOobject
    (
        "magV",
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mag(v)
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
    mesh,
    dimensionedTensor("zero", dimensionSet(0,0,1,0,0,0,0), tensor::zero)
);

dimensionedScalar Eta (eps * M / (R * T));

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
    v*p*M/(R*T)
);

q.write();


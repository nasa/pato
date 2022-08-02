/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "SolidDisplacementModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SolidDisplacementModel::SolidDisplacementModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleSolidMechanicsModel(mesh, dictName),
mesh_(mesh),
dictName_(dictName),
gasPropertiesModel_(meshLookupOrConstructModel<simpleGasPropertiesModel>(mesh,dictName,"GasProperties")),
materialPropertiesModel_(meshLookupOrConstructModel<simpleMaterialPropertiesModel>(mesh,dictName,"MaterialProperties")),
pyrolysisModel_(meshLookupOrConstructModel<simplePyrolysisModel>(mesh,dictName,"Pyrolysis")),
p(meshLookupOrConstructScalar(mesh,"p")),
T(meshLookupOrConstructScalar(mesh,"Ta")),
M(gasPropertiesModel_.M()),
mu(gasPropertiesModel_.mu()),
eps_g(gasPropertiesModel_.eps_g()),
rho_g(gasPropertiesModel_.rho_g()),
R(constant::physicoChemical::R),
K(materialPropertiesModel_.K()),
piTotal(pyrolysisModel_.piTotal()),
Gamma
(
    IOobject
    (
        "Gamma",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    p * M / (mu * R * T) * K
),
Gamma_symm
(
    IOobject
    (
        "Gamma_symm",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    symm(Gamma)
),
Beta
(
    IOobject
    (
        "Beta",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    eps_g * M / (R * T)
),
U(meshLookupOrConstructVector(mesh, "U", dimensionedVector("U", dimLength/dimTime, vector(0,0,0)))),

/*****/

D(meshLookupOrConstructVector(mesh, "D", dimensionedVector("D", dimLength, vector(0,0,0)))),
nu(simpleSolidMechanicsModel::nu_),
E(simpleSolidMechanicsModel::E_),
mu_sM    //sM for solidMechanics
(
    IOobject
    (
        "mu_sM",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    E / (2 * (1 + nu))
),
lambda_sM
(
    IOobject
    (
        "lambda_sM",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    E * nu / ((1 + nu) * (1 - nu))    //if planeStress ; for now, only planeStress
),
threeK
(
    IOobject
    (
        "threeK",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    E / (1 - nu)  //if planeStress ; for now, only planeStress
),
alpha(simpleSolidMechanicsModel::alpha_),
threeKalpha
(
    IOobject
    (
        "threeKalpha",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    threeK * alpha
),
tau_(pyrolysisModel_.tau()),
rho(simpleSolidMechanicsModel::rho_), //otherwise, solidRho from materialPropeertiesModel, but is a list of densities
sigmaD
(
    IOobject
    (
        "sigmaD",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
        /* runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE */
    ),
    mu_sM * twoSymm(fvc::grad(D)) + lambda_sM * (I * tr(fvc::grad(D)))
),
sigma
(
    IOobject
    (
        "sigma",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    rho * sigmaD
),
divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
),
sigmaEq
(
    IOobject
    (
        "sigmaEq",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    sqrt((3.0/2.0)*magSqr(dev(sigma)))
),
cellCenter(meshLookupOrConstructVector(mesh, "cellCenter", dimensionedVector("cellCenter", dimLength, vector(0,0,0)))),
cellCenterX(meshLookupOrConstructScalar(mesh, "cellCenterX", dimensionedScalar("cellCenterX", dimLength, 0))),
cellCenterY(meshLookupOrConstructScalar(mesh, "cellCenterY", dimensionedScalar("cellCenterY", dimLength, 0))),
cellCenterZ(meshLookupOrConstructScalar(mesh, "cellCenterZ", dimensionedScalar("cellCenterZ", dimLength, 0))),
DisplacementX(meshLookupOrConstructScalar(mesh, "DisplacementX", dimensionedScalar("DisplacementX", dimLength, 0))),
DisplacementY(meshLookupOrConstructScalar(mesh, "DisplacementY", dimensionedScalar("DisplacementY", dimLength, 0))),
DisplacementZ(meshLookupOrConstructScalar(mesh, "DisplacementZ", dimensionedScalar("DisplacementZ", dimLength, 0)))

{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SolidDisplacementModel::~SolidDisplacementModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SolidDisplacementModel::update()
{
  // semi-implicit formulation of mass conservation (with Darcy's law -ie. ,
  // momentum conservation- substituted inside the mass conservation equation)

  Info<< "\nCalculating displacement field\n" << endl;

  beforeSolve();

  /*****/

  scalar initialResidual = 0;
// do
  {

    {
      fvVectorMatrix DEqn
      (
          fvm::d2dt2(D)
          ==
          fvm::laplacian(2*mu_sM + lambda_sM, D, "laplacian(DD,D)")
          + divSigmaExp
      );

      //if (thermalStress)
      {
        //const volScalarField& T = Tptr();
        DEqn += fvc::grad(threeKalpha*T);
      }

      //if (pyrolysisStress)
      {
        //Info << "pyro" << -threeK*Chi/6 << nl;
        DEqn -= fvc::grad(threeK*(1-tau_)/6);  // because (1-tau) = Chi
      }


      // DEqn.setComponentReference(1, 0, vector::X, 0);
      // DEqn.setComponentReference(1, 0, vector::Z, 0);

      initialResidual = DEqn.solve().max().initialResidual();

      /*     if (!compactNormalStress)
           {
               divSigmaExp = fvc::div(DEqn.flux());
           }
      */
    }

    {
//      volTensorField gradD(fvc::grad(D));
      sigmaD = mu_sM*twoSymm(fvc::grad(D)) + (lambda_sM*I)*tr(fvc::grad(D));

      // if (compactNormalStress)
      {
        divSigmaExp = fvc::div
                      (
                          sigmaD - (2*mu_sM + lambda_sM)*fvc::grad(D),
                          "div(sigmaD)"
                      );
      }
      //else
      /*    {
              divSigmaExp += fvc::div(sigmaD);
          } */
    }
  } // while (initialResidual > 1e-6); // while (initialResidual > convergenceTolerance && ++iCorr < nCorr);


  //T.correctBoundaryConditions();

  afterSolve();
  moveMesh();
}

void Foam::SolidDisplacementModel::beforeSolve()
{
  /*  const dictionary& stressControl = mesh.solutionDict().subDict("stressAnalysis");

    Switch compactNormalStress(stressControl.lookup("compactNormalStress"));

    int iCorr = 0;
    scalar initialResidual = 0;

    int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);
    convergenceTolerance = readScalar(stressControl.lookup("D"));
    stressControl.lookup("compactNormalStress") >> compactNormalStress;*/
}

void Foam::SolidDisplacementModel::afterSolve()
{

  //if (thermalStress)
  {
//            const volScalarField& T = Tptr();
    sigma = sigma - I*(rho*(threeKalpha*T));
  }

  //if (pyrolysisStress)
  {
    //const volScalarField& T = Tptr();
    sigma = sigma + I*(rho*(threeK*(1-tau_)/6));
  }

  sigmaEq = sqrt((3.0/2.0)*magSqr(dev(sigma)));

  Info<< "Max sigmaEq = " << max(sigmaEq).value()
      << endl;

}

void Foam::SolidDisplacementModel::moveMesh()
{
  Info << "Moving mesh using least squares interpolation" << endl;

// if (dynamicMesh)   //changer de nom
  {
    volPointInterpolation pointInterpolation(mesh_);

    // Create point mesh
    pointMesh pMesh(mesh_);

    wordList types
    (
        pMesh.boundary().size(),
        calculatedFvPatchVectorField::typeName
    );

    pointVectorField pointDU
    (
        IOobject
        (
            "pointDU",
            mesh_.time().timeName(),
            mesh_
        ),
        pMesh,
        dimensionedVector("zero", dimLength, vector::zero),
        types
    );


    tmp<volVectorField> deltaD = D - D.oldTime();
    //Info << "deltaD" << D - D.oldTime() << nl;
    //Info << "D" << D << nl;
    //Info << "D_old" << D.oldTime() << nl;
    pointInterpolation.interpolate(deltaD(), pointDU);
    //Info << "pointDU" <<pointDU<< nl;

    const tmp<vectorField> pointDUI_tmp =
        pointDU.internalField();
    const vectorField& pointDUI = pointDUI_tmp();

    //- Move mesh
    vectorField newPoints = mesh_.points();

    forAll (pointDUI, pointI) {
      newPoints[pointI] += pointDUI[pointI];
    }

    // Correct symmetryPlane points

    forAll(mesh_.boundaryMesh(), patchI) {
      if (isA<symmetryPolyPatch>(mesh_.boundaryMesh()[patchI])) {
        const labelList& meshPoints =
            mesh_.boundaryMesh()[patchI].meshPoints();

        vector avgN =
            gAverage(mesh_.boundaryMesh()[patchI].pointNormals());

        vector i(1, 0, 0);
        vector j(0, 1, 0);
        vector k(0, 0, 1);

        if (mag(avgN&i) > 0.95) {
          forAll(meshPoints, pI) {
            newPoints[meshPoints[pI]].x() = 0;
          }
        } else if (mag(avgN&j) > 0.95) {
          forAll(meshPoints, pI) {
            newPoints[meshPoints[pI]].y() = 0;
          }
        } else if (mag(avgN&k) > 0.95) {
          forAll(meshPoints, pI) {
            newPoints[meshPoints[pI]].z() = 0;
          }
        }
      }
    }

//#   include "calcUnusedNewPoints.H"

//    twoDPointCorrector twoDCorrector(mesh);
//    twoDCorrector.correctPoints(newPoints);

    dynamicFvMesh& mesh = (dynamicFvMesh&)(this->mesh_);
    mesh.movePoints(newPoints);

    cellCenter = mesh.C();

    cellCenterX = cellCenter.component(vector::X);
    cellCenterY = cellCenter.component(vector::Y);
    cellCenterZ = cellCenter.component(vector::Z);

    DisplacementX = D.component(vector::X);
    DisplacementY = D.component(vector::Y);
    DisplacementZ = D.component(vector::Z);

//    stressMesh.V00();
    mesh.moving(false);
    /*}
    else
    {
        FatalErrorIn(args.executable())
            << "Negative Jacobian"
            << exit(FatalError);
    } */
  }

}
// ************************************************************************* //

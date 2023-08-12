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

#include "elasticSolidFoamSolidMechanicsModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::elasticSolidFoamSolidMechanicsModel::elasticSolidFoamSolidMechanicsModel
(
    const Foam::fvMesh& mesh_foam,
    const Foam::word& regionName
)
  :
Foam::simpleSolidMechanicsModel(mesh_foam, regionName),
update_every_n_iter(materialDict_.subDict("SolidMechanics").lookupOrDefault<int>("update_every_n_iter",1)),
iter(0),
dynamicMesh(isA<dynamicFvMesh>(mesh_foam)),
mesh_extend
(
    lookup_foam_extend_mesh(mesh_foam,regionName)
),
mesh(mesh_extend),
runTime(const_cast<Foam_extend_::Time&>(mesh.time())),
g
(
    Foam_extend_::IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        Foam_extend_::IOobject::MUST_READ,
        Foam_extend_::IOobject::NO_WRITE
    )
),
U
(
    Foam_extend_::IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        Foam_extend_::IOobject::MUST_READ,
        Foam_extend_::IOobject::AUTO_WRITE
    ),
    mesh
),
gradU
(
    Foam_extend_::IOobject
    (
        "grad(U)",
        runTime.timeName(),
        mesh,
        Foam_extend_::IOobject::NO_READ,
        Foam_extend_::IOobject::NO_WRITE
    ),
    mesh,
    Foam_extend_::dimensionedTensor("zero", Foam_extend_::dimless, Foam_extend_::tensor::zero)
),
snGradU
(
    Foam_extend_::IOobject
    (
        "snGrad(U)",
        runTime.timeName(),
        mesh,
        Foam_extend_::IOobject::NO_READ,
        Foam_extend_::IOobject::NO_WRITE
    ),
    mesh,
    Foam_extend_::dimensionedVector("zero", Foam_extend_::dimless, Foam_extend_::vector::zero)
),
V
(
    Foam_extend_::IOobject
    (
        "V",
        runTime.timeName(),
        mesh,
        Foam_extend_::IOobject::READ_IF_PRESENT,
        Foam_extend_::IOobject::AUTO_WRITE
    ),
    Foam_extend_::fvc::ddt(U)
),
gradV(Foam_extend_::fvc::ddt(gradU)),
snGradV((snGradU - snGradU.oldTime())/runTime.deltaT()),
epsilon
(
    Foam_extend_::IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        Foam_extend_::IOobject::READ_IF_PRESENT,
        Foam_extend_::IOobject::AUTO_WRITE
    ),
    mesh,
    Foam_extend_::dimensionedSymmTensor("zero", Foam_extend_::dimless, Foam_extend_::symmTensor::zero)
),
sigma
(
    Foam_extend_::IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        Foam_extend_::IOobject::READ_IF_PRESENT,
        Foam_extend_::IOobject::AUTO_WRITE
    ),
    mesh,
    Foam_extend_::dimensionedSymmTensor("zero", Foam_extend_::dimForce/Foam_extend_::dimArea, Foam_extend_::symmTensor::zero)
),
divSigmaExp
(
    Foam_extend_::IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        Foam_extend_::IOobject::NO_READ,
        Foam_extend_::IOobject::NO_WRITE
    ),
    mesh,
    Foam_extend_::dimensionedVector("zero", Foam_extend_::dimForce/Foam_extend_::dimVolume, Foam_extend_::vector::zero)
),
rheology(sigma, U),
rho(rheology.rho()),
mu(rheology.mu()),
lambda(rheology.lambda()),
muf(rheology.muf()),
lambdaf(rheology.lambdaf()),
n(mesh.Sf()/mesh.magSf()),
aitkenDelta
(
    Foam_extend_::IOobject
    (
        "aitkenDelta",
        runTime.timeName(),
        mesh,
        Foam_extend_::IOobject::NO_READ,
        Foam_extend_::IOobject::NO_WRITE
    ),
    mesh,
    Foam_extend_::dimensionedVector("zero", Foam_extend_::dimLength, Foam_extend_::vector::zero)
),
aitkenInitialRes(1.0),
aitkenTheta(0.1)
{
# include "using_extend.H"
# include "Foam_extend_createHistory.H"
# include "Foam_extend_readDivSigmaExpMethod.H"
  solidInterfaceCorr = rheology.solidInterfaceActive();
  {
    if (solidInterfaceCorr) {
      //Info << "Creating solid interface correction" << endl;
      //solidInterfacePtr = new solidInterface(mesh, rheology);
      // constitutiveModel is now in charge of solidInterface
      solidInterfacePtr = &rheology.solInterface();
      //solidInterfacePtr->modifyProperties(muf, lambdaf);

      //- solidInterface needs muf and lambdaf to be used for divSigmaExp
      if (divSigmaExpMethod != "surface" && divSigmaExpMethod != "decompose") {
        FatalError
            << "divSigmaExp must be decompose or surface when "
            << "solidInterface is on"
            << exit(FatalError);
      }

      // check grad scheme
      if
      (
          word(mesh.schemesDict().gradSchemes().lookup("grad(U)"))
          != "leastSquaresSolidInterface"
      ) {
        Warning
            << "The grad(U) gradScheme should be "
            << "leastSquaresSolidInterface for the solidInterface "
            << " procedure to be correct"
            << endl;
      }
    }
  }
  modelInitialized();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::elasticSolidFoamSolidMechanicsModel::~elasticSolidFoamSolidMechanicsModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::elasticSolidFoamSolidMechanicsModel::update()
{
  iter++;
  if (iter%update_every_n_iter==0) {
    iter=0;
#   include "using_extend.H"
#   include "Foam_extend_readSolidMechanicsControls.H"

    int iCorr = 0;
    lduSolverPerformance solverPerf;
    scalar initialResidual = 1.0;
    scalar relativeResidual = 1.0;
//         lduMatrix::debug = 0;

    if (predictor) {
      Info<< "\nPredicting U, gradU and snGradU based on V,"
          << "gradV and snGradV\n" << endl;
      U += V*runTime.deltaT();
      gradU += gradV*runTime.deltaT();
      snGradU += snGradV*runTime.deltaT();
    }

    do {
      U.storePrevIter();

#   include "Foam_extend_calculateDivSigmaExp.H"

      // linear momentum equation
      fvVectorMatrix UEqn
      (
          rho*Foam_extend_::fvm::d2dt2(U)
          ==
          Foam_extend_::fvm::laplacian(2*muf + lambdaf, U, "laplacian(DU,U)")
          + divSigmaExp
          + rho*g
      );

      if (solidInterfaceCorr) {
        solidInterfacePtr->correct(UEqn);
      }

      solverPerf = UEqn.solve();

      if (iCorr == 0) {
        initialResidual = solverPerf.initialResidual();
        aitkenInitialRes = gMax(mag(U.internalField()));
      }

      if (aitkenRelax) {
#       include "Foam_extend_aitkenRelaxation.H"
      } else {
        U.relax();
      }

      gradU = Foam_extend_::fvc::grad(U);

#   include "Foam_extend_calculateRelativeResidual.H"

      if (iCorr % infoFrequency == 0) {
        Info<< "\tTime " << runTime.value()
            << ", Corrector " << iCorr
            << ", Solving for " << U.name()
            << " using " << solverPerf.solverName()
            << ", res = " << solverPerf.initialResidual()
            << ", rel res = " << relativeResidual;

        if (aitkenRelax) {
          Info<< ", aitken = " << aitkenTheta;
        }
        Info<< ", inner iters = " << solverPerf.nIterations() << endl;
      }
    } while
    (
        iCorr++ == 0
        ||
        (
            solverPerf.initialResidual() > convergenceTolerance
            && iCorr < nCorr
        )
    );

    Info<< nl << "Time " << runTime.value() << ", Solving for " << U.name()
        << ", Initial residual = " << initialResidual
        << ", Final residual = " << solverPerf.initialResidual()
        << ", Relative residual = " << relativeResidual
        << ", No outer iterations " << iCorr
        << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << endl;

    if (predictor) {
      V = Foam_extend_::fvc::ddt(U);
      gradV = Foam_extend_::fvc::ddt(gradU);
      snGradV = (snGradU - snGradU.oldTime())/runTime.deltaT();
    }

# include "Foam_extend_calculateEpsilonSigma.H"
# include "Foam_extend_writeFields.H"
# include "Foam_extend_writeHistory.H"

  }
}

// ************************************************************************* //

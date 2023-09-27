/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2021 PATO team
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

#include "ElasticThermalSolidMechanicsModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ElasticThermalSolidMechanicsModel::ElasticThermalSolidMechanicsModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleSolidMechanicsModel(mesh, regionName),
D(createVolField<vector>("D", dimensionedVector("D", dimLength, vector(0,0,0)))),
deltaD(createVolField<vector>("deltaD", dimensionedVector("deltaD", dimLength, vector(0,0,0)))),
gradD(createVolField<tensor>("gradD",fvc::grad(D))),
nu(createVolField<scalar>("nu",dimensionedScalar("0", dimless, scalar(0)))),
E(createVolField<scalar>("E",dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), scalar(0)))),
alpha(createVolField<scalar>("alpha",dimensionedScalar("0", dimensionSet(0,0,0,-1,0,0,0), scalar(0)))),
xi(createVolField<scalar>("xi",dimensionedScalar("0", dimensionSet(0,0,0,0,0,0,0), scalar(0)))),
mu_sM(createVolField<scalar>("mu_sM",E / (2 * (1 + nu)))),
lambda_sM(createVolField<scalar>("lambda_sM",nu*E/((1.0 + nu)*(1.0 - 2.0*nu)))),
threeK(createVolField<scalar>("threeK",E/(1 - 2*nu))),
threeKalpha(createVolField<scalar>("threeKalpha",threeK * alpha)),
threeKxi(createVolField<scalar>("threeKxi",threeK * xi / 3)),
sigma(createVolField<symmTensor>("sigma",mu_sM * twoSymm(gradD) + lambda_sM * (I * tr(gradD)))),
divSigmaExp(createVolField<vector>("divSigmaExp",fvc::div(sigma))),
sigmaEq(createVolField<scalar>("sigmaEq",sqrt((3.0/2.0)*magSqr(dev(sigma))))),
epsilon(createVolField<symmTensor>("epsilon",symm(gradD))),
epsilonEq(createVolField<scalar>("EpsilonEq",sqrt((2.0/3.0)*magSqr(dev(epsilon))))),
solidControl_(mesh.solutionDict().subDict("solidMechanics")),
nCorr_(readInt(solidControl_.lookup("nCorrectors"))),
convergenceTolerance_(readScalar(solidControl_.lookup("D"))),
planeStress_(simpleSolidMechanicsModel::planeStress_),
pyrolysisModel_(refModel<simplePyrolysisModel>()),
tau_(pyrolysisModel_.refVolField<scalar>("tau")),
energyModel_(refModel<simpleEnergyModel>()),
rho(energyModel_.refVolField<scalar>("rho_s")),
T(energyModel_.refVolField<scalar>("Ta")),
T0(createVolField<scalar>("T0",dimensionedScalar("0", dimensionSet(0,0,0,1,0,0,0),simpleSolidMechanicsModel::materialDict_.subDict("SolidMechanics").template lookupOrDefault<int>("T0",300)))),
initialized_(false)
{
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ElasticThermalSolidMechanicsModel::~ElasticThermalSolidMechanicsModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ElasticThermalSolidMechanicsModel::initialize()
{
  mu_sM = E / (2 * (1 + nu));
  lambda_sM = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));
  threeK = E/(1 - 2*nu);
  threeKalpha = threeK * alpha;
  threeKxi = threeK * xi / 3;
  sigma = mu_sM * twoSymm(gradD) + lambda_sM * (I * tr(gradD));
  divSigmaExp = fvc::div(sigma);
  sigmaEq = sqrt((3.0/2.0)*magSqr(dev(sigma)));
  epsilon = symm(gradD);
  epsilonEq = sqrt((2.0/3.0)*magSqr(dev(epsilon)));
}

void Foam::ElasticThermalSolidMechanicsModel::update()
{
  if (!initialized_) {
    initialize();
    initialized_=true;
  }

  // semi-implicit formulation of mass conservation (with Darcy's law -ie. ,
  // momentum conservation- substituted inside the mass conservation equation)

  Info<< "\nCalculating displacement field\n" << endl;

  beforeSolve();

  int iCorr = 0;
  scalar initialResidual = 1.0;
  scalar relResD = 1.0;

  do {
    D.storePrevIter();

    sigma = mu_sM * twoSymm(gradD) + lambda_sM * (I * tr(gradD));

    //Compact normal stress
    divSigmaExp =  fvc::div(sigma - (2*mu_sM + lambda_sM) * gradD);
    fvVectorMatrix DEqn
    (
        fvm::d2dt2(rho,D)
        ==
        fvm::laplacian(2*mu_sM + lambda_sM, D, "laplacian(DD,D)")
        + divSigmaExp
        - fvc::grad(threeKalpha*(T - T0))           // Thermal stress
        + fvc::grad(threeKxi*(1 - tau_))         // Pyrolysis stress
    );

    if (iCorr == 0) {
      initialResidual = DEqn.solve().max().initialResidual();
    } else {
      DEqn.solve();
    }

    D.relax();

    gradD = fvc::grad(D);

    // displacement
    scalar maxDD =
        gMax(mag(D.primitiveField() - D.oldTime().primitiveField()));
    relResD =
        gMax
        (
            mag
            (
                D.primitiveField() - D.prevIter().primitiveField()
            )
            /(maxDD + SMALL)
        );
    iCorr++;
  } while(iCorr == 0 || (relResD > convergenceTolerance_ && iCorr <nCorr_));

  Info<< "Solved for " << D.name()
      << " in " << iCorr << " iterations"
      << ", initial res = " << initialResidual
      << ", final rel res = " << relResD << nl
      << endl;

  sigma =  mu_sM * twoSymm(gradD) + lambda_sM * (I * tr(gradD));

  // Compact normal stress
  divSigmaExp = fvc::div(sigma - (2 * mu_sM + lambda_sM) *  gradD);

  afterSolve();
}

void Foam::ElasticThermalSolidMechanicsModel::beforeSolve()
{
  if (planeStress_) {
    lambda_sM = nu*E/((1.0 + nu)*(1.0 - nu));
    threeK = E/(1 - nu);
  } else {
    lambda_sM = nu*E/((1.0 + nu)*(1.0 -2*nu));
    threeK = E/(1 - 2*nu);
  }
  threeKalpha = threeK*alpha;
  threeKxi = threeK*xi/3;
}

void Foam::ElasticThermalSolidMechanicsModel::afterSolve()
{
  sigma = sigma + I * threeKxi*(1 - tau_);
  sigma = sigma - I * threeKalpha*(T - T0);
  epsilon = symm(gradD);
  sigmaEq = sqrt((3.0/2.0)*magSqr(dev(sigma)));
  epsilonEq = sqrt((2.0/3.0)*magSqr(dev(epsilon)));

  Info<< "Max sigmaEq = " << max(sigmaEq).value()
      << endl;
  Info<< "Max epsilonEq = " << max(epsilonEq).value()
      << endl;

  // Only needed for fluid solid interaction tutorial flapping console
  tmp<volVectorField> deltaD_tmp = D - D.oldTime();
  deltaD = deltaD_tmp();
}

// ************************************************************************* //

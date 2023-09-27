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

#include "DarcyForchheimerLawMassModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DarcyForchheimerLawMassModel::DarcyForchheimerLawMassModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleMassModel(mesh, regionName),
Beta(createVolField<tensor>("Beta",dimensionedTensor("0",dimensionSet(0,-1,0,0,0,0,0),tensor(0,0,0,0,0,0,0,0,0)))),
g(createUniformField<vector>("g",dimensionedVector("0", dimLength/(dimTime*dimTime), vector::zero))),
K(createVolField<tensor>("K",dimensionedTensor("0",dimLength*dimLength,tensor::zero))),
mDotG(createVolField<vector>("mDotG",dimensionedVector("0", dimMass/pow(dimLength,2)/dimTime, vector::zero))),
mDotGFace(createSurfaceField<vector>("mDotGface",linearInterpolate(mDotG))),
mDotGw(createVolField<scalar>("mDotGw",dimensionedScalar("0", dimMass/pow(dimLength,2)/dimTime, 0.0))),
p(createVolField<scalar>("p")),
U(createVolField<vector>("U",dimensionedVector("0",dimLength/dimTime,vector::zero))),
v(createVolField<vector>("vG",dimensionedVector("0", dimLength/dimTime, vector::zero))),
R(constant::physicoChemical::R),
energyModel_(refModel<simpleEnergyModel>()),
T(energyModel_.refVolField<scalar>("Ta")),
gasPropertiesModel_(refModel<simpleGasPropertiesModel>()),
eps_g(gasPropertiesModel_.refVolField<scalar>("eps_g")),
M(gasPropertiesModel_.refVolField<scalar>("M_g")),
mu(gasPropertiesModel_.refVolField<scalar>("mu_g")),
rho_g(gasPropertiesModel_.refVolField<scalar>("rho_g")),
pyrolysisModel_(refModel<simplePyrolysisModel>()),
piTotal(pyrolysisModel_.refVolField<scalar>("piTotal")),
Gamma(createVolField<tensor>("Gamma", dimensionedTensor("zero", dimensionSet(0,0,1,0,0,0,0), tensor::zero))),
Gamma_symm(createVolField<symmTensor>("Gamma_symm", symm(Gamma))),
Eta(createVolField<scalar>("Eta", eps_g * M / (R * T))),
Gamma_G(createVolField<tensor>("Gamma_G", Gamma * p * M / (R * T))),
phiG(createSurfaceField<scalar>("phiG", linearInterpolate(Gamma_G & g) & mesh.Sf()))
{
  U = eps_g*v;
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DarcyForchheimerLawMassModel::~DarcyForchheimerLawMassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DarcyForchheimerLawMassModel::update()
{
  // semi-implicit formulation of mass conservation (with Darcy's law -ie. ,
  // momentum conservation- substituted inside the mass conservation equation)

  beforeSolve();

  // compute implicitly the new pressure
  if(!simpleMassModel::dynamicMesh_) {
    solve
    (
        fvm::ddt(Eta, p)                                 // storage
        - fvm::laplacian(Gamma, p)                         // convection
        - piTotal                                          // source : pyrolysis
        + fvc::div(phiG)                                   // convection : buoyancy part
    );
  } else {
    solve
    (
        fvm::ddt(Eta, p)                                   // storage
        - fvm::div(fvc::interpolate(Eta)*mesh_.phi(), p)   // mesh motion correction (ALE)
        - fvm::laplacian(Gamma, p)                         // convection
        - piTotal                                          // source : pyrolysis
        + fvc::div(phiG)                                   // convection : buoyancy part
    );
  }

  const_cast<volScalarField&>(T).correctBoundaryConditions();
  afterSolve();
}

void Foam::DarcyForchheimerLawMassModel::beforeSolve()
{
  if(simpleMassModel::debug_) {
    Info<<"\t update Eta"<<endl;
  }
  // update explicit quantities
  Eta.ref() = eps_g()*M()/(R * T());
  Eta.boundaryFieldRef() =
      eps_g.boundaryField()*M.boundaryField()
      /(R.value()*T.boundaryField());

  if(simpleMassModel::debug_) {
    Info<<"\t update Gamma"<<endl;
  }
  Gamma.replace
  (
      0,
      K.component(Tensor<double>::XX)
      /(
          mu + Beta.component(Tensor<double>::XX)*p*M/(R*T)*mag(U)
          *K.component(Tensor<double>::XX)
      )
  );
  Gamma.replace
  (
      1,
      K.component(Tensor<double>::XY)
      /(
          mu + Beta.component(Tensor<double>::XY)*p*M/(R*T)*mag(U)
          *K.component(Tensor<double>::XY)
      )
  );
  Gamma.replace
  (
      2,
      K.component(Tensor<double>::XZ)
      /(
          mu + Beta.component(Tensor<double>::XZ)*p*M/(R*T)*mag(U)
          *K.component(Tensor<double>::XZ)
      )
  );
  Gamma.replace
  (
      3,
      K.component(Tensor<double>::YX)
      /(
          mu + Beta.component(Tensor<double>::YX)*p*M/(R*T)*mag(U)
          *K.component(Tensor<double>::YX)
      )
  );
  Gamma.replace
  (
      4,
      K.component(Tensor<double>::YY)
      /(
          mu + Beta.component(Tensor<double>::YY)*p*M/(R*T)*mag(U)
          *K.component(Tensor<double>::YY)
      )
  );
  Gamma.replace
  (
      5,
      K.component(Tensor<double>::YZ)
      /(
          mu + Beta.component(Tensor<double>::YZ)*p*M/(R*T)*mag(U)
          *K.component(Tensor<double>::YZ)
      )
  );
  Gamma.replace
  (
      6,
      K.component(Tensor<double>::ZX)
      /(
          mu + Beta.component(Tensor<double>::ZX)*p*M/(R*T)*mag(U)
          *K.component(Tensor<double>::ZX)
      )
  );
  Gamma.replace
  (
      7,
      K.component(Tensor<double>::ZY)
      /(
          mu + Beta.component(Tensor<double>::ZY)*p*M/(R*T)*mag(U)
          *K.component(Tensor<double>::ZZ)
      )
  );
  Gamma.replace
  (
      8,
      K.component(Tensor<double>::ZZ)
      /(
          mu + Beta.component(Tensor<double>::ZZ)*p*M/(R*T)*mag(U)
          *K.component(Tensor<double>::ZZ)
      )
  );
  Gamma.field() *= p()*M()/(R.value()*T());
  Gamma.boundaryFieldRef() *=
      p.boundaryField()*M.boundaryField()/(R.value()*T.boundaryField());
  Gamma_symm=symm(Gamma);

  if(simpleMassModel::debug_) {
    Info<< "\t Update Gamma_G" << endl;
  }
  Gamma_G.ref() = Gamma()*p()*M()/(R*T());
  Gamma_G.boundaryFieldRef() =
      Gamma.boundaryFieldRef() * p.boundaryField() * M.boundaryField()
      / (R.value() * T.boundaryField());

  if(simpleMassModel::debug_) {
    Info<<"\t update p"<<endl;
  }
}

void Foam::DarcyForchheimerLawMassModel::afterSolve()
{
  // print min and max pressures in the terminal during the run
  Info << "p_max = "
       << max(p).value()
       << " / p_min = "
       << min(p).value()
       << endl;

  if(simpleMassModel::debug_) {
    Info<<"\t update rho_g"<<endl;
  }
  // gas density (Ideal gas law)
  rho_g = p*M/(R*T);

  if(simpleMassModel::debug_) {
    Info<<"\t update v"<<endl;
  }
  // It's not really Gamma that is updated but the velocity part of Gamma
  Gamma.replace
  (
      0,
      K.component(Tensor<double>::XX)
      /(
          mu + Beta.component(Tensor<double>::XX)*rho_g*mag(U)
          *K.component(Tensor<double>::XX)
      )
  );
  Gamma.replace
  (
      1,
      K.component(Tensor<double>::XY)
      /(
          mu + Beta.component(Tensor<double>::XY)*rho_g*mag(U)
          *K.component(Tensor<double>::XY)
      )
  );
  Gamma.replace
  (
      2,
      K.component(Tensor<double>::XZ)
      /(
          mu + Beta.component(Tensor<double>::XZ)*rho_g*mag(U)
          *K.component(Tensor<double>::XZ)
      )
  );
  Gamma.replace
  (
      3,
      K.component(Tensor<double>::YX)
      /(
          mu + Beta.component(Tensor<double>::YX)*rho_g*mag(U)
          *K.component(Tensor<double>::YX)
      )
  );
  Gamma.replace
  (
      4,
      K.component(Tensor<double>::YY)
      /(
          mu + Beta.component(Tensor<double>::YY)*rho_g*mag(U)
          *K.component(Tensor<double>::YY)
      )
  );
  Gamma.replace
  (
      5,
      K.component(Tensor<double>::YZ)
      /(
          mu + Beta.component(Tensor<double>::YZ)*rho_g*mag(U)
          *K.component(Tensor<double>::YZ)
      )
  );
  Gamma.replace
  (
      6,
      K.component(Tensor<double>::ZX)
      /(
          mu + Beta.component(Tensor<double>::ZX)*rho_g*mag(U)
          *K.component(Tensor<double>::ZX)
      )
  );
  Gamma.replace
  (
      7,
      K.component(Tensor<double>::ZY)
      /(
          mu + Beta.component(Tensor<double>::ZY)*rho_g*mag(U)
          *K.component(Tensor<double>::ZZ)
      )
  );
  Gamma.replace
  (
      8,
      K.component(Tensor<double>::ZZ)
      /(
          mu + Beta.component(Tensor<double>::ZZ)*rho_g*mag(U)
          *K.component(Tensor<double>::ZZ)
      )
  );
  Gamma.field() *= p()*M()/(R.value()*T());
  Gamma.boundaryFieldRef() *=
      p.boundaryField()*M.boundaryField()/(R.value()*T.boundaryField());

  /* average gas velocity (Darcy-Forchheimer's law)
      from the new pressure gradient (cell centered).
      This formulation is not optimized.*/
  v = - 1./eps_g*(((Gamma*R*T/(p*M)) & fvc::grad(p))
                  - (((Gamma*R*T/(p*M))*rho_g) & g));

  // velocity used in the fluid region
  U = eps_g*v;

  if(simpleMassModel::debug_) {
    Info<<"\t update mDotG"<<endl;
  }
  // pyrolysis gas mass flow rate (cell centered)
  mDotG = eps_g*rho_g*v;

  if(simpleMassModel::debug_) {
    Info<<"\t update mDotGFace"<<endl;
  }

  // pyrolysis gas mass flow rate (interpolated on the faces)
  mDotGFace =
      - fvc::interpolate(rho_g/mu)
      * fvc::interpolate(K) & mesh_.Sf()/mesh_.magSf()
      * fvc::snGrad(p);

  forAll(mDotGw.boundaryField(), patchI) {
    forAll(mDotGw.boundaryField()[patchI], faceI) {
      const vector nf =
          - mesh_.Sf().boundaryField()[patchI][faceI]
          / mesh_.magSf().boundaryField()[patchI][faceI];
      const scalar mDotGFace_BF = mDotGFace.boundaryField()[patchI][faceI]&(-nf);
      mDotGw.boundaryFieldRef()[patchI][faceI] = mDotGFace_BF;
    }
  }

  if(simpleMassModel::debug_) {
    Info<<"\t Update phiG" << endl;
  }
  // Gravitational mass flow rate
  phiG = linearInterpolate(Gamma_G & g) & mesh_.Sf();

}

// ************************************************************************* //

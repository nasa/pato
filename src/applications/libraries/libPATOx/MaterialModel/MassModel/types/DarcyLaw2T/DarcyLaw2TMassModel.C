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

#include "DarcyLaw2TMassModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DarcyLaw2TMassModel::DarcyLaw2TMassModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleMassModel(mesh, regionName),
p(createVolField<scalar>("p")),
K(createVolField<tensor>("K",dimensionedTensor("0",dimLength*dimLength,tensor::zero))),
v(createVolField<vector>("vG",dimensionedVector("0", dimLength/dimTime, vector::zero))),
U(createVolField<vector>("U",dimensionedVector("0", dimLength/dimTime, vector::zero))),
g(createUniformField<vector>("g",dimensionedVector("0", dimLength/(dimTime*dimTime), vector::zero))),
gasPropertiesModel(refModel<simpleGasPropertiesModel>()),
M(gasPropertiesModel.refVolField<scalar>("M_g")),
mu(gasPropertiesModel.refVolField<scalar>("mu_g")),
eps_g(gasPropertiesModel.refVolField<scalar>("eps_g")),
rho_g(gasPropertiesModel.refVolField<scalar>("rho_g")),
pyrolysisModel(refModel<simplePyrolysisModel>()),
piTotal(pyrolysisModel.refVolField<scalar>("piTotal")),
energyModel(refModel<simpleEnergyModel>()),
Tg_(energyModel.refVolField<scalar>("Tg")),
R(constant::physicoChemical::R),
Gamma(createVolField<tensor>("Gamma", p * M / (mu * R * Tg_) * K)),
Gamma_symm(createVolField<symmTensor>("Gamma_symm", symm(Gamma))),
Eta(createVolField<scalar>("Eta", eps_g * M / (R * Tg_))),
Gamma_G(createVolField<tensor>("Gamma_G", Gamma * p * M / (R * Tg_))),
mDotG(createVolField<vector>("mDotG",dimensionedVector("0", dimMass/pow(dimLength,2)/dimTime, vector::zero))),
mDotGFace(createSurfaceField<vector>("mDotGface",linearInterpolate(mDotG))),
mDotGw(createVolField<scalar>("mDotGw",dimensionedScalar("0", dimMass/pow(dimLength,2)/dimTime, 0.0))),
phiG(createSurfaceField<scalar>("phiG", linearInterpolate(Gamma_G & g) & mesh.Sf()))
{
  U = eps_g*v;
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DarcyLaw2TMassModel::~DarcyLaw2TMassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DarcyLaw2TMassModel::update()
{
  // semi-implicit formulation of mass conservation (with Darcy's law -ie. ,
  // momentum conservation- substituted inside the mass conservation equation)

  beforeSolve();

  // compute implicitly the new pressure
  if(!simpleMassModel::dynamicMesh_) {
    solve
    (
        fvm::ddt(Eta, p)                                  // storage
        - fvm::laplacian(Gamma, p)                         // convection
        - piTotal                                         // source : pyrolysis
        + fvc::div(phiG)                                   // convection : buoyancy part
    );
  } else {
    solve
    (
        fvm::ddt(Eta, p)                                  // storage
        - fvm::div(fvc::interpolate(Eta)*mesh_.phi(), p)  // mesh motion correction (ALE)
        - fvm::laplacian(Gamma, p)                         // convection
        - piTotal                                          // source : pyrolysis
        + fvc::div(phiG)                                   // convection : buoyancy part
    );
  }

  const_cast<volScalarField&>(Tg_).correctBoundaryConditions();

  afterSolve();
}

void Foam::DarcyLaw2TMassModel::beforeSolve()
{
  if(simpleMassModel::debug_) {
    Info<<"\t update Beta"<<endl;
  }
  // update explicit quantities
  Eta.ref()  = eps_g() * M() / (R * Tg_());
  Eta.boundaryFieldRef()  = eps_g.boundaryField() * M.boundaryField() / (R.value() * Tg_.boundaryField());

  if(simpleMassModel::debug_) {
    Info<<"\t update Gamma"<<endl;
  }
  Gamma.ref() = p() * M() / (mu() * R * Tg_()) * K();
  Gamma.boundaryFieldRef() = p.boundaryField() * M.boundaryField() / (mu.boundaryField() * R.value() * Tg_.boundaryField()) * K.boundaryField();
  Gamma_symm=symm(Gamma);

  if(simpleMassModel::debug_) {
    Info<< "\t Update Gamma_G" << endl;
  }
  Gamma_G.ref() = Gamma.ref() * p() * M() / (R * Tg_());
  Gamma_G.boundaryFieldRef() =
      Gamma.boundaryFieldRef() * p.boundaryField() * M.boundaryField()
      / (R.value() * Tg_.boundaryField());

  if(simpleMassModel::debug_) {
    Info<<"\t update p"<<endl;
  }
}

void Foam::DarcyLaw2TMassModel::afterSolve()
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
  rho_g = p * M / (R * Tg_);

  if(simpleMassModel::debug_) {
    Info<<"\t update v"<<endl;
  }
  // average gas velocity (Darcy's law) from the new pressure gradient (cell centered)
  v = - 1./eps_g*1./mu*((K & fvc::grad(p)) - rho_g*(K & g));

  // velocity used in the fluid region
  U = eps_g*v;

  if(simpleMassModel::debug_) {
    Info<<"\t update mDotG"<<endl;
  }
  // pyrolysis gas mass flow rate (cell centered)
  mDotG = eps_g * rho_g * v;

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

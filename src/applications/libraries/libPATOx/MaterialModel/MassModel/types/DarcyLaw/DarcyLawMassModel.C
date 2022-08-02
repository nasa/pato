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

#include "DarcyLawMassModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DarcyLawMassModel::DarcyLawMassModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleMassModel(mesh, dictName),
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
v(simpleMassModel::vG_),
U(meshLookupOrConstructVector(mesh, "U", dimensionedVector("U", dimLength/dimTime, vector(0,0,0)))),
mDotG(simpleMassModel::mDotG_),
mDotGFace(simpleMassModel::mDotGFace_)
{
  U = eps_g*v;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DarcyLawMassModel::~DarcyLawMassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DarcyLawMassModel::update()
{
  // semi-implicit formulation of mass conservation (with Darcy's law -ie. ,
  // momentum conservation- substituted inside the mass conservation equation)

  beforeSolve();

  // compute implicitly the new pressure
  if(!simpleMassModel::dynamicMesh_) {
    solve
    (
        fvm::ddt(Beta, p)                                  // storage
        - fvm::laplacian(Gamma, p)                         // convection
        - piTotal                                          // source : pyrolysis
    );
  } else {
    solve
    (
        fvm::ddt(Beta, p)                                  // storage
        - fvm::div(fvc::interpolate(Beta)*mesh_.phi(), p)   // mesh motion correction (ALE)
        - fvm::laplacian(Gamma, p)                         // convection
        - piTotal                                          // source : pyrolysis
    );
  }

  T.correctBoundaryConditions();

  afterSolve();
}

void Foam::DarcyLawMassModel::beforeSolve()
{
  if(simpleMassModel::debug_) {
    Info<<"\t update Beta"<<endl;
  }
  // update explicit quantities
  Beta.ref()  = eps_g() * M() / (R * T());
  Beta.boundaryFieldRef()  = eps_g.boundaryField() * M.boundaryField() / (R.value() * T.boundaryField());

  if(simpleMassModel::debug_) {
    Info<<"\t update Gamma"<<endl;
  }
  Gamma.ref() = p() * M() / (mu() * R * T()) * K();
  Gamma.boundaryFieldRef() = p.boundaryField() * M.boundaryField() / (mu.boundaryField() * R.value() * T.boundaryField()) * K.boundaryField();
  Gamma_symm=symm(Gamma);

  if(simpleMassModel::debug_) {
    Info<<"\t update p"<<endl;
  }
}

void Foam::DarcyLawMassModel::afterSolve()
{
  // print min and max pressures in the terminal during the run
  Info << "p_max = "
       << max(p).value()
       << " / p_min = "
       << min(p).value()
       << endl;

  if(simpleMassModel::debug_) {
    Info<<"\t update v"<<endl;
  }
  // average gas velocity (Darcy's law) from the new pressure gradient (cell centered)
  v = - 1. / eps_g * 1. / mu * (K & fvc::grad(p));

  // velocity used in the fluid region
  U = eps_g*v;

  if(simpleMassModel::debug_) {
    Info<<"\t update rho_g"<<endl;
  }
  // gas density (Ideal gas law)
  rho_g = p * M / (R * T);

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
      - fvc::interpolate(rho_g / mu) *
      fvc::interpolate(K) & mesh_.Sf() / mesh_.magSf() *
      fvc::snGrad(p);
}

// ************************************************************************* //

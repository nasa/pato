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

#include "Pyrolysis2TEnergyModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Pyrolysis2TEnergyModel::Pyrolysis2TEnergyModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleEnergyModel(mesh, regionName),
fvOptions_(fv::options::New(mesh_)),
cp_s(createVolField<scalar>("cp",dimensionedScalar("0", dimensionSet(0,2,-2,-1,0,0,0), scalar(0)))),
Hv_(createVolField<scalar>("Hv",dimensionedScalar("0",dimensionSet(1,-1,-3,-1,0,0,0),0))),
k_a(createVolField<tensor>("k",dimensionedTensor("0",dimensionSet(1, 1, -3, -1, 0, 0, 0),I))),
pyrolysisFlux_(createVolField<scalar>("pyrolysisFlux",dimensionedScalar("0",dimensionSet(1,-1,-3,0,0,0,0),0))),
rho_s(createVolField<scalar>("rho_s",dimensionedScalar("0",dimMass/dimVolume,0))),
Tg(createVolField<scalar>("Tg")),
Ts(createVolField<scalar>("Ta")),
Dp_(createDimScalarProp("Dp",true,dimensionedScalar("1",dimLength,1))),
Hv0_(createDimScalarProp("Hv0",true,dimensionedScalar("0",dimensionSet(1,-1,-3,-1,0,0,0),0))),
HvType_(createWordProp("HvType",true,"constant")),
R(Foam::constant::physicoChemical::R),
massModel_(refModel<simpleMassModel>()),
g(massModel_.refUniformField<vector>("g")),
K(massModel_.refVolField<tensor>("K")),
p(massModel_.refVolField<scalar>("p")),
gasPropertiesModel_(refModel<simpleGasPropertiesModel>()),
cp_g(gasPropertiesModel_.refVolField<scalar>("cp_g")),
eps_g(gasPropertiesModel_.refVolField<scalar>("eps_g")),
h_g(gasPropertiesModel_.refVolField<scalar>("h_g")),
k_g(gasPropertiesModel_.refVolField<scalar>("k_g")),
M(gasPropertiesModel_.refVolField<scalar>("M_g")),
mu(gasPropertiesModel_.refVolField<scalar>("mu_g")),
rho_g(gasPropertiesModel_.refVolField<scalar>("rho_g")),
vg(massModel_.refVolField<vector>("vG")),
epsgRhog(createVolField<scalar>("epsgRhog",(eps_g*rho_g))),
GammaHg(createVolField<tensor>("GammaHg",((h_g * p * M) / (mu * R * Tg)) * K)),
Gamma_GHg(createVolField<tensor>("Gamma_GHg",GammaHg*p*M/(R*Tg))),
k_s(createVolField<tensor>("k_s",(k_a - I * k_g))),
phiHg(createSurfaceField<scalar>("phiHg", linearInterpolate(Gamma_GHg & g) & mesh.Sf())),
phi_g(createSurfaceField<scalar>("phi_g", linearInterpolate(eps_g*rho_g*vg) & mesh_.Sf()))
{
  if (!fvOptions_.optionList::size()) {
    Info << getTabLevel() << "No Energy finite volume options present" << endl;
  }
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Pyrolysis2TEnergyModel::~Pyrolysis2TEnergyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Pyrolysis2TEnergyModel::update()
{

  beforeSolve();

  // Mesh motion ?
  if(this->dynamicMesh_) {
    // Global energy balance
    // Fluid
    fvScalarMatrix TgEqn
    (
        eps_g*rho_g*cp_g*(fvm::ddt(Tg) - fvm::div(mesh_.phi(), Tg)) // storage - implicit in T, explicit in rhoCp - with mesh motion correction (ALE)
        + h_g*(fvc::ddt(epsgRhog) - fvc::div(mesh_.phi(), epsgRhog))  // gas storage - explicit
        - fvc::ddt(eps_g, p)                                          // Pressure energy
        - fvc::laplacian(GammaHg, p)                                  // convection - explicit
        // - ( (eps_g * vg) & fvc::grad(p) )				// Pressure work
        - fvm::laplacian(k_g, Tg)					// Gas conductivity
        + fvc::div(phiHg)                                             // convection (gravitational part) - explicit
        + fvm::Sp(Hv_, Tg) - Hv_ * Ts				        // Exchange between solid and fluid
        ==
        fvOptions_(rho_g*cp_g, Tg)
    );
    fvOptions_.constrain(TgEqn);
    TgEqn.solve();
    fvOptions_.correct(Tg);

    // Solid
    fvScalarMatrix TsEqn
    (
        rho_s*cp_s*(fvm::ddt(Ts) - fvm::div(mesh_.phi(), Ts)) // storage - implicit in T, explicit in rhoCp - with mesh_ motion correction (ALE)
        + pyrolysisFlux_                                                    // storage (second part of the derivative - mass variation) - solid mass loss by pyrolsyis - explicit
        - fvm::laplacian(k_s, Ts)                                           // conduction - implicit in T, explicit in k
        + Hv_*(Ts - Tg)                                                   // Exchange between solid and fluid
        ==
        fvOptions_(rho_s*cp_s, Ts)
    );
    fvOptions_.constrain(TsEqn);
    TsEqn.solve();
    fvOptions_.correct(Ts);
  } else {
    // Global energy balance
    // Fluid
    fvScalarMatrix TgEqn
    (
        eps_g*rho_g*cp_g*fvm::ddt(Tg) // storage - implicit in T, explicit in rhoCp
        + h_g*fvc::ddt(epsgRhog)        // gas storage - explicit
        - fvc::ddt(eps_g, p)            // Pressure energy
        - fvc::laplacian(GammaHg, p)    // convection - explicit
        // - ((eps_g*vg) & fvc::grad(p))   // Pressure work
        - fvm::laplacian(k_g, Tg)       // Gas conductivity
        + fvc::div(phiHg)               // convection (gravitational part) - explicit
        + fvm::Sp(Hv_, Tg) - Hv_*Ts     // Exchange between solid and fluid
        ==
        fvOptions_(rho_g*cp_g, Tg)      // Optional energy source term
    );
    fvOptions_.constrain(TgEqn);
    TgEqn.solve();
    fvOptions_.correct(Tg);

    // Solid
    fvScalarMatrix TsEqn
    (
        rho_s*cp_s*(fvm::ddt(Ts)) // storage - implicit in T, explicit in rhoCp
        + pyrolysisFlux_                        // storage (second part of the derivative - mass variation) - solid mass loss by pyrolsyis - explicit
        - fvm::laplacian(k_s, Ts)               // conduction - implicit in T, explicit in k
        + Hv_*(Ts - Tg)                       // Exchange between solid and fluid
        ==
        fvOptions_(rho_s*cp_s, Ts)            // Optional energy source term
    );
    fvOptions_.constrain(TsEqn);
    TsEqn.solve();
    fvOptions_.correct(Ts);
  }


  afterSolve();
}

void Foam::Pyrolysis2TEnergyModel::beforeSolve()
{
  epsgRhog.ref() = (eps_g()*rho_g());
  epsgRhog.boundaryFieldRef() = (eps_g.boundaryField()*rho_g.boundaryField());
  GammaHg.ref() = ((h_g()*p()*M())/(mu()*R*Tg()))*K();
  GammaHg.boundaryFieldRef() =
      ((h_g.boundaryField()*p.boundaryField()*M.boundaryField())
       /(mu.boundaryField()*R.value()*Tg.boundaryField()))*K.boundaryField();

  Gamma_GHg.ref() = GammaHg.ref()*p()*M()/(R*Tg());
  Gamma_GHg.boundaryFieldRef() =
      GammaHg.boundaryFieldRef()*p.boundaryField()*M.boundaryField()
      /(R.value()*Tg.boundaryField());

  k_s.ref() = k_a() - I*k_g();
  k_s.boundaryFieldRef() = k_a.boundaryField() - I*k_g.boundaryField();

  // Calculation of heat transfer coefficient if not constant
  if (HvType_ == "Wakao") {
    const volScalarField Re(Dp_*eps_g*mag(vg)*rho_g/mu);
    const volScalarField Pr(cp_g*mu/k_g);
    const volScalarField Nu(2 + 1.1*pow(Pr,1./3)*pow(Re,0.6));
    Hv_ = (6*(1 - eps_g)/Dp_)*k_g*Nu/Dp_;
  } else if (HvType_ == "constant") {
    Hv_ = Hv0_;
  } else {
    FatalErrorInFunction << HvType_ << "is not implemented."
                         << exit(FatalError);
  }
}

void Foam::Pyrolysis2TEnergyModel::afterSolve()
{
  if(this->debug_) {
    forAll(Ts, cellI) {
      if (Ts[cellI]<0) {
        const wordList fieldNames =
            mesh_.objectRegistry::sortedNames("volScalarField");
        forAll(fieldNames, nameI) {
          const volScalarField field =
              const_cast<volScalarField&>
              (
                  mesh_.objectRegistry::lookupObject<volScalarField>
                  (
                      fieldNames[nameI]
                  )
              );
          field.write();
        }
        FatalErrorInFunction
            << "The solid temperature is negative in cell "
            << cellI << "." << exit(FatalError);
      }
    }
    forAll(Ts.boundaryField(), patchI) {
      forAll(Ts.boundaryField()[patchI], faceI) {
        if (Ts.boundaryField()[patchI][faceI]<0) {
          const wordList fieldNames =
              mesh_.objectRegistry::sortedNames("volScalarField");
          forAll(fieldNames, nameI) {
            const volScalarField field =
                const_cast<volScalarField&>
                (mesh_.objectRegistry::lookupObject<volScalarField>
                 (
                     fieldNames[nameI]
                 )
                );
            field.write();
          }
          FatalErrorInFunction
              << "The solid temperature is negative in patch "
              << patchI << " and face " << faceI << "."
              << exit(FatalError);
        }
      }
    }
  }

  if(this->debug_) {
    forAll(Tg, cellI) {
      if (Tg[cellI]<0) {
        const wordList fieldNames =
            mesh_.objectRegistry::sortedNames("volScalarField");
        forAll(fieldNames, nameI) {
          const volScalarField field =
              const_cast<volScalarField&>
              (
                  mesh_.objectRegistry::lookupObject<volScalarField>
                  (
                      fieldNames[nameI]
                  )
              );
          field.write();
        }
        FatalErrorInFunction
            << "The fluid temperature is negative in cell "
            << cellI << "." << exit(FatalError);
      }
    }
    forAll(Tg.boundaryField(), patchI) {
      forAll(Tg.boundaryField()[patchI], faceI) {
        if (Tg.boundaryField()[patchI][faceI]<0) {
          const wordList fieldNames =
              mesh_.objectRegistry::sortedNames("volScalarField");
          forAll(fieldNames, nameI) {
            const volScalarField field =
                const_cast<volScalarField&>
                (
                    mesh_.objectRegistry::lookupObject<volScalarField>
                    (
                        fieldNames[nameI]
                    )
                );
            field.write();
          }

          FatalErrorInFunction
              << "The fluid temperature is negative in patch "
              << patchI << " and face " << faceI << "."
              << exit(FatalError);
        }
      }
    }
  }
}

// ************************************************************************* //

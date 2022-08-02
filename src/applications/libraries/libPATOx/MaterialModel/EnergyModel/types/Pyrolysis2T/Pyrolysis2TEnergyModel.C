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
    along with OpenFOAM.  If Pyrolysist, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Pyrolysis2TEnergyModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Pyrolysis2TEnergyModel::Pyrolysis2TEnergyModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleEnergyModel(mesh, dictName),
mesh_(mesh),
dictName_(dictName),
MassModel_(meshLookupOrConstructModel<simpleMassModel>(mesh,dictName,"Mass")),
materialPropertiesModel_(meshLookupOrConstructModel<simpleMaterialPropertiesModel>(mesh,dictName,"MaterialProperties")),

materialPropertiesDirectory(fileName(simpleEnergyModel::materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand()),
constantPropertiesDictionary
(
    IOobject
    (
        materialPropertiesDirectory+"/constantProperties",
        mesh.time().db().parent(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
gasPropertiesModel_(meshLookupOrConstructModel<simpleGasPropertiesModel>(mesh,dictName,"GasProperties")),
pyrolysisModel_(meshLookupOrConstructModel<simplePyrolysisModel>(mesh,dictName,"Pyrolysis")),
h_g(gasPropertiesModel_.h_g()),
p(meshLookupOrConstructScalar(mesh,"p")),
M(gasPropertiesModel_.M()),
mu(gasPropertiesModel_.mu()),

R(Foam::constant::physicoChemical::R),
Tg(meshLookupOrConstructScalar(mesh,"Tg")),
Ts(meshLookupOrConstructScalar(mesh,"Ta")),
K(materialPropertiesModel_.K()),

eps_g(gasPropertiesModel_.eps_g()),
rho_g(gasPropertiesModel_.rho_g()),
rho_s(materialPropertiesModel_.rho_s()),
cp_g(gasPropertiesModel_.cp_g()),
cp_s(materialPropertiesModel_.cp()),
pyrolysisFlux_(materialPropertiesModel_.pyrolysisFlux()),
vg(MassModel_.vG()),
k_g(gasPropertiesModel_.k_g()),
k_a(materialPropertiesModel_.k()),
k_s
(
    IOobject
    (
        "k_s",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    (k_a - I * k_g)
),
epsgRhog
(
    IOobject
    (
        "epsgRhog",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),

    (eps_g*rho_g)
),
phi_g
(
    IOobject
    (
        "phi_g",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    linearInterpolate(eps_g*rho_g*vg) & mesh_.Sf()
),
GammaHg
(
    IOobject
    (
        "GammaHg",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    ((h_g * p * M) / (mu * R * Tg)) * K
),
Hv0_(constantPropertiesDictionary.lookup("Hv0"))

{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Pyrolysis2TEnergyModel::~Pyrolysis2TEnergyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Pyrolysis2TEnergyModel::update()
{

  beforeSolve();  

  phi_g = linearInterpolate(eps_g*rho_g*vg) & mesh_.Sf();

// Solid

  if(this->dynamicMesh_) {
    // global energy balance
    solve
    (
        (1.-eps_g) * rho_s * cp_s * (fvm::ddt(Ts) - fvm::div(mesh_.phi(), Ts))       	// storage - implicit in T, explicit in rhoCp - with mesh_ motion correction (ALE)
        + pyrolysisFlux_                                              			// storage (second part of the derivative - mass variation) - solid mass loss by pyrolsyis - explicit
        - fvm::laplacian(k_s, Ts)                                        		// conduction - implicit in T, explicit in k
        + fvm::Sp(Hv0_, Ts) - Hv0_ * Tg							// Exchange between solid and fluid
    );
  } else {
    // global energy balance
    solve
    (
        (1.-eps_g) * rho_s * cp_s * (fvm::ddt(Ts))					// storage - implicit in T, explicit in rhoCp - with mesh motion correction (ALE)
        + pyrolysisFlux_                                            			// storage (second part of the derivative - mass variation) - solid mass loss by pyrolsyis - explicit
        - fvm::laplacian(k_s, Ts)                                     			// conduction - implicit in T, explicit in k
        + fvm::Sp(Hv0_, Ts) - Hv0_ * Tg							   // Exchange between solid and fluid
    );

  }

// Fluid

  if(this->dynamicMesh_) {
    // global energy balance
    solve
    (
        eps_g * rho_g * cp_g * (fvm::ddt(Tg) - fvm::div(mesh_.phi(), Tg))      	// storage - implicit in T, explicit in rhoCp - with mesh motion correction (ALE)
	+ h_g * (fvc::ddt(epsgRhog) - fvc::div(mesh_.phi(), epsgRhog))     	// gas storage - explicit
	- fvc::ddt(eps_g, p)							// Pressure energy
         - fvc::laplacian(GammaHg, p)                                  		// convection - explicit
//	+ cp_g * fvm::div(phi_g, Tg)						// convection - implicit
//	- ( (eps_g * vg) & fvc::grad(p) )					// Pressure work
	- fvm::laplacian(k_g, Tg)						// Gas conductivity
        + fvm::Sp(Hv0_, Tg) - Hv0_ * Ts						// Exchange between solid and fluid
    );
  } else {
    // global energy balance
    solve
    (
        eps_g * rho_g * cp_g * fvm::ddt(Tg)                                 	// storage - implicit in T, explicit in rhoCp - with mesh motion correction (ALE)
	+ h_g * fvc::ddt(epsgRhog)                                     		// gas storage - explicit
	- fvc::ddt(eps_g, p)							// Pressure energy
        - fvc::laplacian(GammaHg, p)                                  		// convection - explicit
	//+ cp_g * fvm::div(phi_g, Tg)						// convection - implicit
//	- ( (eps_g * vg) & fvc::grad(p) )					// Pressure work
	- fvm::laplacian(k_g, Tg)						// Gas conductivity
        + fvm::Sp(Hv0_, Tg) - Hv0_ * Ts						// Exchange between solid and fluid
    );

  }

  afterSolve();
}

void Foam::Pyrolysis2TEnergyModel::beforeSolve()
{
  epsgRhog.ref() =  (eps_g()*rho_g());
  epsgRhog.boundaryFieldRef() =  (eps_g.boundaryField()*rho_g.boundaryField());

  GammaHg.ref() =  ((h_g() * p() * M()) / (mu() * R * Tg())) * K();
  GammaHg.boundaryFieldRef() =  ((h_g.boundaryField() * p.boundaryField() * M.boundaryField()) / (mu.boundaryField() * R.value() * Tg.boundaryField())) * K.boundaryField();
}

void Foam::Pyrolysis2TEnergyModel::afterSolve()
{
  if(this->debug_) {

    forAll(Ts, cellI) {
      if (Ts[cellI]<0) {
        const wordList fieldNames = mesh_.objectRegistry::sortedNames("volScalarField");
        forAll(fieldNames, nameI) {
          const volScalarField field = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(fieldNames[nameI]));
          field.write();
        }
        FatalErrorInFunction << "The solid temperature is negative in cell " << cellI << "." << exit(FatalError);
      }
    }
    forAll(Ts.boundaryField(), patchI) {
      forAll(Ts.boundaryField()[patchI], faceI) {
        if (Ts.boundaryField()[patchI][faceI]<0) {
          const wordList fieldNames = mesh_.objectRegistry::sortedNames("volScalarField");
          forAll(fieldNames, nameI) {
            const volScalarField field = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(fieldNames[nameI]));
            field.write();
          }
          FatalErrorInFunction << "The solid temperature is negative in patch " << patchI << " and face " << faceI << "." << exit(FatalError);
        }
      }
    }
  }

  if(this->debug_) {

    forAll(Tg, cellI) {
      if (Tg[cellI]<0) {
        const wordList fieldNames = mesh_.objectRegistry::sortedNames("volScalarField");
        forAll(fieldNames, nameI) {
          const volScalarField field = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(fieldNames[nameI]));
          field.write();
        }
        FatalErrorInFunction << "The fluid temperature is negative in cell " << cellI << "." << exit(FatalError);

      }
    }
    forAll(Tg.boundaryField(), patchI) {
      forAll(Tg.boundaryField()[patchI], faceI) {
        if (Tg.boundaryField()[patchI][faceI]<0) {
          const wordList fieldNames = mesh_.objectRegistry::sortedNames("volScalarField");
          forAll(fieldNames, nameI) {
            const volScalarField field = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(fieldNames[nameI]));
            field.write();
          }

          FatalErrorInFunction << "The fluid temperature is negative in patch " << patchI << " and face " << faceI << "." << exit(FatalError);

        }
      }
    }
  }
}

// ************************************************************************* //

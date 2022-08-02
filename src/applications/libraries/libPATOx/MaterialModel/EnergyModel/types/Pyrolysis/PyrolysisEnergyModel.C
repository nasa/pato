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

#include "PyrolysisEnergyModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PyrolysisEnergyModel::PyrolysisEnergyModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleEnergyModel(mesh, dictName),
mesh_(mesh),
dictName_(dictName),
materialPropertiesModel_(meshLookupOrConstructModel<simpleMaterialPropertiesModel>(mesh,dictName,"MaterialProperties")),
gasPropertiesModel_(meshLookupOrConstructModel<simpleGasPropertiesModel>(mesh,dictName,"GasProperties")),
pyrolysisModel_(meshLookupOrConstructModel<simplePyrolysisModel>(mesh,dictName,"Pyrolysis")),
h_g(gasPropertiesModel_.h_g()),
p(meshLookupOrConstructScalar(mesh,"p")),
M(gasPropertiesModel_.M()),
mu(gasPropertiesModel_.mu()),
R(Foam::constant::physicoChemical::R),
T(meshLookupOrConstructScalar(mesh,"Ta")),
K(materialPropertiesModel_.K()),
eps_g(gasPropertiesModel_.eps_g()),
rho_g(gasPropertiesModel_.rho_g()),
rho_s(materialPropertiesModel_.rho_s()),
cp(materialPropertiesModel_.cp()),
pyrolysisFlux_(materialPropertiesModel_.pyrolysisFlux()),
epsgRhogEg
(
    IOobject
    (
        "epsgRhogEg",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    (eps_g*rho_g*h_g) - (eps_g*p)
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
    ((h_g * p * M) / (mu * R * T)) * K
),
k(materialPropertiesModel_.k())
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PyrolysisEnergyModel::~PyrolysisEnergyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PyrolysisEnergyModel::update()
{
  beforeSolve();

  if(this->dynamicMesh_) {
    // global energy balance
    solve
    (
        rho_s * cp * (fvm::ddt(T) - fvm::div(mesh_.phi(), T))         // storage - implicit in T, explicit in rhoCp - with mesh_ motion correction (ALE)
        + pyrolysisFlux_                                              // storage (second part of the derivative - mass variation) - solid mass loss by pyrolsyis - explicit
        + fvc::ddt(epsgRhogEg) -  fvc::div(mesh_.phi(), epsgRhogEg)   // gas storage - explicit
        - fvm::laplacian(k, T)                                        // conduction - implicit in T, explicit in k
        - fvc::laplacian(GammaHg, p)                                  // convection - explicit
    );
  } else {
    // global energy balance
    solve
    (
        rho_s * cp * (fvm::ddt(T))                                 // storage - implicit in T, explicit in rhoCp - with mesh motion correction (ALE)
        + pyrolysisFlux_                                           // storage (second part of the derivative - mass variation) - solid mass loss by pyrolsyis - explicit
        + fvc::ddt(epsgRhogEg)                                     // gas storage - explicit
        - fvm::laplacian(k, T)                                     // conduction - implicit in T, explicit in k
        - fvc::laplacian(GammaHg, p)                               // convection - explicit
    );

  }

  afterSolve();
}

void Foam::PyrolysisEnergyModel::beforeSolve()
{
  epsgRhogEg.ref() =  (eps_g()*rho_g()*h_g()) - (eps_g()*p());
  epsgRhogEg.boundaryFieldRef() =  (eps_g.boundaryField()*rho_g.boundaryField()*h_g.boundaryField()) - (eps_g.boundaryField()*p.boundaryField());

  GammaHg.ref() =  ((h_g() * p() * M()) / (mu() * R * T())) * K();
  GammaHg.boundaryFieldRef() =  ((h_g.boundaryField() * p.boundaryField() * M.boundaryField()) / (mu.boundaryField() * R.value() * T.boundaryField())) * K.boundaryField();
}

void Foam::PyrolysisEnergyModel::afterSolve()
{
  if(this->debug_) {
    // Find if temperature < 0
    int cellID=0;
    int print_fields=0;
    forAll(T, cellI) {
      if (T[cellI]<0) {
        cellID=cellI+1;
        print_fields++;
        break;
      }
    }
    // cellID per processor
    List<int> cellIDs(Pstream::nProcs());
    for(int i = 0; i<Pstream::nProcs(); i++) {
      cellIDs[i]=0;
    }
    cellIDs[Pstream::myProcNo()]=cellID;
    if (Pstream::parRun()) {
      reduce(print_fields,sumOp<int>());
      reduce(cellIDs,sumOp<List<int>>());
    }
    // if temperature < 0, print the fields and raise an error.
    if (print_fields>0) {
      const wordList fieldNames = mesh_.objectRegistry::sortedNames("volScalarField");
      forAll(fieldNames, nameI) {
        const volScalarField field = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(fieldNames[nameI]));
        field.write();
      }
      List<int> cellIDs_;
      List<int> ranks_;
      forAll(cellIDs, i) {
        if (cellIDs[i]>0) {
          cellIDs_.append(cellIDs[i]-1);
          ranks_.append(i);
        }
      }
      if (Pstream::parRun()) {
        if (Pstream::myProcNo()==ranks_[0]) {
          FatalErrorInFunction << "The temperature is negative in cell " << cellIDs_[0] << " from rank " << ranks_[0] << "." << exit(FatalError);
        }
      } else {
        FatalErrorInFunction << "The temperature is negative in cell " << cellID << "." << exit(FatalError);
      }
    }
    int patchID=0;
    int faceID=0;
    print_fields=0;
    forAll(T.boundaryField(), patchI) {
      forAll(T.boundaryField()[patchI], faceI) {
        if (T.boundaryField()[patchI][faceI]<0) {
          patchID=patchI+1;
          faceID=faceI+1;
          print_fields++;
          break;
        }
      }
    }
    List<int> patchIDs(Pstream::nProcs());
    List<int> faceIDs(Pstream::nProcs());
    for(int i = 0; i<Pstream::nProcs(); i++) {
      patchIDs[i]=0;
      faceIDs[i]=0;
    }
    patchIDs[Pstream::myProcNo()]=patchID;
    faceIDs[Pstream::myProcNo()]=faceID;
    if (Pstream::parRun()) {
      reduce(print_fields,sumOp<int>());
      reduce(patchIDs,sumOp<List<int>>());
      reduce(faceIDs,sumOp<List<int>>());
    }
    if (print_fields>0) {
      List<int> patchIDs_;
      List<int> faceIDs_;
      List<int> ranks_;
      forAll(patchIDs, i) {
        if (patchIDs[i]>0) {
          patchIDs_.append(patchIDs[i]-1);
          ranks_.append(i);
        }
        if (faceIDs[i]>0) {
          faceIDs_.append(faceIDs[i]-1);
        }
      }
      const wordList fieldNames = mesh_.objectRegistry::sortedNames("volScalarField");
      forAll(fieldNames, nameI) {
        const volScalarField field = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(fieldNames[nameI]));
        field.write();
      }
      if (Pstream::parRun()) {
        if (Pstream::myProcNo()==ranks_[0]) {
          FatalErrorInFunction << "The temperature is negative in patch " << patchIDs_[0] << " and face " << faceIDs_[0] << " from rank " << ranks_[0] << "." << exit(FatalError);
        }
      } else {
        FatalErrorInFunction << "The temperature is negative in patch " << patchID << " and face " << faceID << "." << exit(FatalError);
      }
    }
  }
}

// ************************************************************************* //

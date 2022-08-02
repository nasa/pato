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

#include "massIOModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massIOModel::massIOModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
noIOModel(mesh, dictName),
recession(meshLookupOrConstructScalar(mesh, "recession", dimensionedScalar("0", dimLength, scalar(0.0)))),
rho_s(meshLookupOrConstructScalar(mesh, "rho_s")),
mDotCw(meshLookupOrConstructScalar(mesh, "mDotCw", dimensionedScalar("0", dimMass/pow(dimLength,2)/dimTime, scalar(0.0)))),
MassModel(meshLookupOrConstructModel<simpleMassModel>(mesh,dictName,"Mass")),
mDotGFace(MassModel.mDotGFace()),
PyrolysisModel(meshLookupOrConstructModel<simplePyrolysisModel>(mesh,dictName,"Pyrolysis")),
rho_v(PyrolysisModel.rho_v()),
rho_c(PyrolysisModel.rho_c()),
topPatchName(simpleIOModel::materialDict_.subDict("IO").lookup("topPatchName")),
topPatchID(noIOModel::mesh_.boundaryMesh().findPatchID(topPatchName)),
massOutputFile(simpleIOModel::materialDict_.subDict("IO").lookup("massOutputFile")),
os_mass(massOutputFile),
massLossOutputFile(simpleIOModel::materialDict_.subDict("IO").lookup("massLossOutputFile")),
os_massLoss(massLossOutputFile),
rho_s0(0.0),
V0(0.0)
{
  if (topPatchID<0) {
    FatalError << topPatchName << " patch not found." << exit(FatalError);
  }
  if (mesh.boundaryMesh()[topPatchID].size()!=1) {
    FatalError << "massIOModel works only in 1D. The top patch size must equal to 1." << exit(FatalError);
  }
  initialPosition = mesh_.Cf().boundaryField()[topPatchID][0];

  os_mass << "// t(s)" << "	" << "m_dot_g_surf(kg/m²/s)" << " " << "m_dot_c(kg/m²/s)"
          << " "  << "< 0.98 Virgin" << " " << "> 1.02 Char" << " " << "Total recession (m)" << nl;
  os_mass << "0" << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << nl;

  os_massLoss << "// t(s) Mean density(kg/m³) V/V0 m/m0" << endl;

  forAll(rho_s, i) {
    V0+=mesh.V()[i];
    rho_s0+=mesh.V()[i]*rho_s[i];
  }
  rho_s0=rho_s0/V0;

  scalar VtimeName = 0;
  scalar rho_s_mean = 0;
  forAll(rho_s, i) {
    VtimeName+=mesh.V()[i];
    rho_s_mean+=mesh.V()[i]*rho_s[i];
  }
  rho_s_mean = rho_s_mean / VtimeName;
  os_massLoss << mesh.time().timeName()  << " " << rho_s_mean << " " << 1 << " " << 1 << endl;

  scalar rho_v_ = rho_v[0];
  forAll(rho_v, cellI) {
    if(rho_v_!=rho_v[cellI]) {
      FatalErrorInFunction << "rho_v is not constant." << exit(FatalError);
    }
  }
  forAll(rho_v.boundaryField(), patchI) {
    forAll(rho_v.boundaryField()[patchI], faceI) {
      if(rho_v_!=rho_v.boundaryField()[patchI][faceI]) {
        FatalErrorInFunction << "rho_v is not constant." << exit(FatalError);
      }
    }
  }

  scalar rho_c_ = rho_c[0];
  forAll(rho_c, cellI) {
    if(rho_c_!=rho_c[cellI]) {
      FatalErrorInFunction << "rho_c is not constant." << exit(FatalError);
    }
  }
  forAll(rho_c.boundaryField(), patchI) {
    forAll(rho_c.boundaryField()[patchI], faceI) {
      if(rho_c_!=rho_c.boundaryField()[patchI][faceI]) {
        FatalErrorInFunction << "rho_c is not constant." << exit(FatalError);
      }
    }
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::massIOModel::~massIOModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::massIOModel::update()
{
  noIOModel::writeOutput();

  if (noIOModel::mesh_.time().outputTime() || noIOModel::mesh_.time().timeName() == "0") {
    int rhoSize = rho_s.size() - 1;

    scalar rho_c_2 = rho_c[0] + 0.02 * (rho_v[0] - rho_c[0]);
    scalar rho_v_98 = rho_c[0] + 0.98 * (rho_v[0] - rho_c[0]);

    const vector nf =
        - mesh_.Sf().boundaryField()[topPatchID][0]
        / mesh_.magSf().boundaryField()[topPatchID][0];

    const scalar mDotGw = mDotGFace.boundaryField()[topPatchID][0]&(-nf);

    if ((rho_s[rhoSize] >= rho_v_98)) {
      os_mass << noIOModel::mesh_.time().timeName() << " " << mDotGw << " " << mDotCw.boundaryField()[topPatchID][0]
              << " " << 0 << " " << 0 << " " << recession.boundaryField()[topPatchID][0] << endl;
    } else {
      vector char_front = noIOModel::mesh_.C()[rhoSize];
      vector virgin_front = noIOModel::mesh_.C()[rhoSize];
      int index_rho = 0;
      while (rho_s[index_rho] > rho_v_98) {
        index_rho++;
      }
      if (rho_s[index_rho] == rho_v_98) {
        virgin_front = noIOModel::mesh_.C()[index_rho];
      } else {
        scalar weight_front = (rho_s[index_rho] - rho_v_98) / (rho_s[index_rho] - rho_s[index_rho - 1]);
        virgin_front = (1 - weight_front) * noIOModel::mesh_.C()[index_rho] + weight_front * noIOModel::mesh_.C()[index_rho - 1];
      }
      if ((rho_s[rhoSize] < rho_c_2)) {
        index_rho--;
        while (rho_s[index_rho] > rho_c_2) {
          index_rho++;
        }
        if (rho_s[index_rho] == rho_c_2) {
          char_front = noIOModel::mesh_.C()[index_rho];
        } else {
          scalar weight_front = (rho_s[index_rho] - rho_c_2) / (rho_s[index_rho] - rho_s[index_rho - 1]);
          char_front = (1 - weight_front) * noIOModel::mesh_.C()[index_rho] + weight_front * noIOModel::mesh_.C()[index_rho - 1];
        }
      }

      os_mass << noIOModel::mesh_.time().timeName() << " " << mDotGw << " " << mDotCw.boundaryField()[topPatchID][0]
              << " " << mag(initialPosition-virgin_front) << " " << mag(initialPosition-char_front) << " " << recession.boundaryField()[topPatchID][0] << endl;
    }

    // Mass loss
    scalar VtimeName = 0;
    scalar rho_s_mean = 0;
    forAll(rho_s, i) {
      VtimeName+=noIOModel::mesh_.V()[i];
      rho_s_mean+=noIOModel::mesh_.V()[i]*rho_s[i];
    }
    rho_s_mean = rho_s_mean / VtimeName;
    os_massLoss <<  noIOModel::mesh_.time().timeName()  << " " << rho_s_mean << " " << VtimeName/V0 << " " << rho_s_mean/rho_s0 << endl;
  }
}


// ************************************************************************* //

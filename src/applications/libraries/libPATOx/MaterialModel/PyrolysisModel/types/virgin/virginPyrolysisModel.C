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
    along with OpenFOAM.  If LinearArrheniust, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "virginPyrolysisModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::virginPyrolysisModel::virginPyrolysisModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simplePyrolysisModel(mesh, regionName),
nSolidPhases_(createScalarProp("nSolidPhases","yes",1)),
piTotal_(createVolField<scalar>("piTotal",dimensionedScalar("0",dimMass/dimVolume/dimTime,0))),
tau(createVolField<scalar>("tau",dimensionedScalar("1", dimless, 1.0),wordList(nPatches_,"zeroGradient"))),
rho_v_(createVolField<scalar>("rho_v",dimensionedScalar("0",dimMass/dimVolume,0))),
rho_c_(createVolField<scalar>("rho_c",dimensionedScalar("0",dimMass/dimVolume,0))),
oneVolScalarField(createVolField<scalar>("one",dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0))),
sizeXsi(0)
{
  if (this->debug_) {
    Info << getTabLevel() << "debug: start --- Foam::virginPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
    Info << getTabLevel() << "debug: solidEps_, solidEpsI_, solidRho_, solidRhoI_ --- Foam::virginPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  // Initialize solidEps_, solidEpsI_, solidRho_, solidRhoI_
  // The solid phases are defined from value 1 (0 is attributed to the gas).
  for(int i = 0; i < nSolidPhases_; i++) {
    word name_epsI = "epsI_s["+Foam::name(i+1)+"]";
    word name_eps = "eps_s["+Foam::name(i+1)+"]";
    dimensionedScalar solidEpsI_value = createDimScalarProp("epsI[" + std::to_string(i+1) + "]",true,dimensionedScalar("0",dimless,0.0));
    solidEpsI_.append(createVolField<scalar>(name_epsI.c_str(),solidEpsI_value));
    solidEps_.append(createVolField<scalar>(name_eps.c_str(),solidEpsI_value));
    word name_rhoI = "rhoI_s["+Foam::name(i+1)+"]";
    word name_rho = "rho_s["+Foam::name(i+1)+"]";
    dimensionedScalar solidRhoI_value = createDimScalarProp("rhoI[" + std::to_string(i+1) + "]",true,dimensionedScalar("0",dimMass/dimVolume,0.0));
    solidRhoI_.append(createVolField<scalar>(name_rhoI.c_str(),solidRhoI_value));
    solidRho_.append(createVolField<scalar>(name_rho.c_str(),solidRhoI_value));
  }
  if (this->debug_) {
    Info << getTabLevel() << "debug: nSolidPhase and nPyroReac --- Foam::virginPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  // Init nSolidPhase and nPyroReac
  if (nSolidPhases_ < 1) {
    FatalErrorInFunction
        << "The number of solid phases 'nSolidPhases' must be > 0. "
        "If not, just use a CFD solver!" << nl
        << exit(FatalError);
  }
  if (this->debug_) {
    Info << getTabLevel() << "debug: nPyroReac --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  nPyroReac.resize(nSolidPhases_);
  forAll(nPyroReac, solidI) {
    const word namePyro = "nPyroReac[" + Foam::name(solidI + 1) + "]";
    nPyroReac.set(solidI, createScalarProp(namePyro,"yes")); // default value = 0
    sizeXsi += nPyroReac[solidI];
  }
  if (this->debug_) {
    Info << getTabLevel() << "debug: hashXsi --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  hashXsi.setSize(sizeXsi, 2);
  int nindex = 0;
  for (int solidI = 0; solidI < nSolidPhases_; solidI++) {
    for (int reacI = 0; reacI < nPyroReac[solidI]; reacI++) {
      hashXsi[nindex + reacI][0] = solidI + 1;
      hashXsi[nindex + reacI][1] = reacI + 1;
    }
    nindex += nPyroReac[solidI];
  }
  // Init Fp
  for(int i = 0; i < sizeXsi; i++) {
    word phase_reac="[" + std::to_string(hashXsi[i][0]) + "][" + std::to_string(hashXsi[i][1]) + "]";
    Fp.append(createDimScalarProp("F"+phase_reac,"yes"));
  }
  initialize();
  if (this->debug_) {
    Info << getTabLevel() << "debug: rho_v_ and rho_c --- Foam::virginPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  // this model is initialized
  modelInitialized();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::virginPyrolysisModel::~virginPyrolysisModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::virginPyrolysisModel::update()
{}

void Foam::virginPyrolysisModel::initialize()
{
  // Init rho_v and rho_c
  rho_v_ *= 0;
  forAll(solidRho_, phaseI) {
    rho_v_ += solidEpsI_[phaseI]*solidRhoI_[phaseI];
  }
  rho_c_ *= 0;
  for (int i = 0; i < nSolidPhases_; i++) {
    rho_c_ += solidRhoI_[i] * solidEpsI_[i];
  }
  for(int i = 0; i < sizeXsi; i++) {
    rho_c_ -=
        Fp[i] *
        solidEpsI_[hashXsi[i][0]-1] *
        solidRhoI_[hashXsi[i][0]-1];
  }
}

// ************************************************************************* //

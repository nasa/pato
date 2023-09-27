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
    along with OpenFOAM.  If Tabulated2Tt, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Tabulated2TGasPropertiesModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Tabulated2TGasPropertiesModel::Tabulated2TGasPropertiesModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleGasPropertiesModel(mesh, dictName),
cp_g(createVolField<scalar>("cp_g",dimensionedScalar("zero", dimensionSet(0, 2, -2, -1, 0, 0, 0), 0.0))),
eps_g(createVolField<scalar>("eps_g",dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0))),
h_g(createVolField<scalar>("h_g",dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0))),
k_g(createVolField<scalar>("k_g",dimensionedScalar("zero", dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0))),
M(createVolField<scalar>("M_g",dimensionedScalar("zero", dimensionSet(1, 0, 0, 0, -1, 0, 0), 0.02))),
mu(createVolField<scalar>("mu_g",dimensionedScalar("mu", dimensionSet(1, -1, -1, 0, 0, 0, 0), 2e-5))),
rho_g(createVolField<scalar>("rho_g",dimensionedScalar("zero", dimensionSet(1, -3, 0, 0, 0, 0, 0), 0.0))),
eps_g_c_(createDimScalarProp("eps_g_c")),
eps_g_v_(createDimScalarProp("eps_g_v")),
energyModel_(refModel<simpleEnergyModel>()),
Tg(energyModel_.refVolField<scalar>("Tg")),
massModel_(refModel<simpleMassModel>()),
p(massModel_.refVolField<scalar>("p")),
pyrolysisModel_(meshLookupOrConstructModel<simplePyrolysisModel>(mesh,dictName,"Pyrolysis")),
tau_(pyrolysisModel_.refVolField<scalar>("tau")),
gasPropertiesFile_(fileName(simpleGasPropertiesModel::materialDict_.subDict("GasProperties").lookup("GasPropertiesFile")).expand())
{
  // read and store the gas-property table into the RAM for faster access
  Info << getTabLevel() << "Reading the gas property table" << nl;

  IFstream gasPropertyInputFile(
      gasPropertiesFile_);  // opens an input file
  IFstream gasPropertyInputFileTemp(
      gasPropertiesFile_);  // opens an input file

  if (gasPropertyInputFile.good() == false) {
    FatalErrorInFunction << "gasProperties file not found" << nl
                         << exit(FatalError);
  }


  Info << getTabLevel() << "The gas-property table must have 7 columns: p, T, M, cp, h, mu, k_g"
       << nl;
  int columnTableGP = 7;

  int rawTableGP = 0;
  scalar tempGP, rawTableFracGP, rawTableIntGP;
  int i_rawGP = 0;

  while (true) {
    gasPropertyInputFileTemp >> tempGP;
    if (gasPropertyInputFileTemp.eof() == 1)
      break;
    i_rawGP++;
  }

  rawTableFracGP =
      modf
      (
          static_cast<scalar>(i_rawGP) /
          static_cast<scalar>(columnTableGP),
          &rawTableIntGP
      );

  if (rawTableFracGP != 0) {

    FatalErrorInFunction << "... and it does not seem to have 7 columns." << nl
                         << rawTableIntGP << " lines and "
                         << rawTableFracGP * columnTableGP << " 'scalar' have been read."
                         << nl
                         << exit(FatalError);

  } else {
    rawTableGP = i_rawGP / columnTableGP;
    Info << simpleModel::getTabLevel() << "The gas-property table has been read ("
         << rawTableGP << " lines)." << nl;
  }

  RectangularMatrix<scalar> TableGP(rawTableGP, columnTableGP);
  for (int x = 0; x < rawTableGP; x++) {
    for (int i = 0; i < columnTableGP; i++) {
      gasPropertyInputFile >> TableGP[x][i];
    }
  }

  // "gasPT" is an object of the class Tabulated2TGasPropertiesObject, which handles the gas property lookup and interpolations when called
  gasPT_ptr.reset(new Tabulated2TGasPropertiesObject(TableGP));
  update();
  modelInitialized();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Tabulated2TGasPropertiesModel::~Tabulated2TGasPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Tabulated2TGasPropertiesModel::update()
{
  if (this->debug_) {
    Info << "--- update gas properties --- Foam::Tabulated2TGasPropertiesModel::update()" << endl;
  }
  Tabulated2TGasPropertiesObject& gasPT = gasPT_ptr();
  // Gas thermodynamics and transport properties

  if (this->debug_) {
    Info << "--- M_g --- Foam::Tabulated2TGasPropertiesModel::update()" << endl;
    Info << "---  min,max p,T " << min(p) << " " << max(p) << " " << min(Tg) << " " << max(Tg) << "--- Foam::Tabulated2TGasPropertiesModel::update()" << endl;
  }
  // Read data from a table
  forAll (M, i) {
    M[i] = gasPT.M(p[i], Tg[i]);
  }

  if (this->debug_) {
    Info << "--- cp_g --- Foam::Tabulated2TGasPropertiesModel::update()" << endl;
  }
  forAll (cp_g, i) {
    cp_g[i] = gasPT.cp(p[i], Tg[i]);
  }


  if (this->debug_) {
    Info << "--- h_g --- Foam::Tabulated2TGasPropertiesModel::update()" << endl;
  }
  forAll (h_g, i) {
    h_g[i] = gasPT.h(p[i], Tg[i]);
  }

  if (this->debug_) {
    Info << "--- mu_g --- Foam::Tabulated2TGasPropertiesModel::update()" << endl;
  }
  forAll (mu, i) {
    mu[i] = gasPT.mu(p[i], Tg[i]);
  }

  if (this->debug_) {

    Info << "--- mu_g --- Foam::Tabulated2TGasPropertiesModel::update()" << endl;
  }
  /*forAll (rho_g, i) {
    rho_g[i] = gasPT.rho(p[i], Tg[i]);
  }*/
  forAll (k_g, i) {
    k_g[i] = gasPT.k(p[i], Tg[i]);
  }

  if (this->debug_) {

    Info << "--- boundaries --- Foam::Tabulated2TGasPropertiesModel::update()" << endl;
  }
  forAll(mesh_.boundaryMesh(), patchI) {
    if (!isA<emptyPolyPatch>(mesh_.boundaryMesh()[patchI])) {
      forAll(mesh_.boundaryMesh()[patchI], faceI) {
        M.boundaryFieldRef()[patchI][faceI] =
            gasPT.M
            (
                p.boundaryField()[patchI][faceI],
                Tg.boundaryField()[patchI][faceI]
            );
        cp_g.boundaryFieldRef()[patchI][faceI] =
            gasPT.cp
            (
                p.boundaryField()[patchI][faceI],
                Tg.boundaryField()[patchI][faceI]
            );
        h_g.boundaryFieldRef()[patchI][faceI] =
            gasPT.h
            (
                p.boundaryField()[patchI][faceI],
                Tg.boundaryField()[patchI][faceI]
            );
        mu.boundaryFieldRef()[patchI][faceI] =
            gasPT.mu
            (
                p.boundaryField()[patchI][faceI],
                Tg.boundaryField()[patchI][faceI]
            );

        /*rho_g.boundaryFieldRef()[patchI][faceI] =
            gasPT.rho
            (
                p.boundaryField()[patchI][faceI],
                Tg.boundaryField()[patchI][faceI]
            );*/
        k_g.boundaryFieldRef()[patchI][faceI] =
            gasPT.k
            (
                p.boundaryField()[patchI][faceI],
                Tg.boundaryField()[patchI][faceI]
            );

      }
    }
  }

  if (this->debug_) {
    Info << "--- rho_g eps_g --- Foam::Tabulated2TGasPropertiesModel::update()" << endl;
  }
  const dimensionedScalar R = constant::physicoChemical::R;
  rho_g  = p * M / (R * Tg); // gas density (perfect gas law)
  eps_g = eps_g_c_ + (eps_g_v_ - eps_g_c_) * tau_; // empirical porosity

  if (this->debug_) {
    Info << "--- end --- Foam::Tabulated2TGasPropertiesModel::update()" << endl;
  }
}

// ************************************************************************* //

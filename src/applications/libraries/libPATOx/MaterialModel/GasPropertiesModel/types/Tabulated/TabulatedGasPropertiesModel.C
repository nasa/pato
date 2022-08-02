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
    along with OpenFOAM.  If Tabulatedt, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "TabulatedGasPropertiesModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TabulatedGasPropertiesModel::TabulatedGasPropertiesModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleGasPropertiesModel(mesh, dictName),
M(simpleGasPropertiesModel::M_),
mu(simpleGasPropertiesModel::mu_),
h_g(simpleGasPropertiesModel::h_g_),
eps_g(simpleGasPropertiesModel::eps_g_),
rho_g(simpleGasPropertiesModel::rho_g_),
Tg(meshLookupOrConstructScalar(mesh,"Ta")),
p(meshLookupOrConstructScalar(mesh,"p")),
gasPropertiesFile_(fileName(simpleGasPropertiesModel::materialDict_.subDict("GasProperties").lookup("GasPropertiesFile")).expand()),
materialPropertiesDirectory(fileName(simpleGasPropertiesModel::materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand()),
constantPropertiesDictionary
(
    IOobject
    (
        "constantProperties",
        materialPropertiesDirectory,
        mesh.time().db().parent(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
eps_g_c_(simpleGasPropertiesModel::materialDict_.subDict("GasProperties").found("eps_g_c")?simpleGasPropertiesModel::materialDict_.subDict("GasProperties").lookup("eps_g_c"):constantPropertiesDictionary.lookup("eps_g_c")),
eps_g_v_(simpleGasPropertiesModel::materialDict_.subDict("GasProperties").found("eps_g_v")?simpleGasPropertiesModel::materialDict_.subDict("GasProperties").lookup("eps_g_v"):constantPropertiesDictionary.lookup("eps_g_v")),
pyrolysisModel_(meshLookupOrConstructModel<simplePyrolysisModel>(mesh,dictName,"Pyrolysis")),
tau_(pyrolysisModel_.tau())
{
  // read and store the gas-property table into the RAM for faster access
  Info << "Reading the gas property table" << nl;

  IFstream gasPropertyInputFile(
      gasPropertiesFile_);  // opens an input file
  IFstream gasPropertyInputFileTemp(
      gasPropertiesFile_);  // opens an input file

  if (gasPropertyInputFile.good() == false) {
    FatalErrorInFunction << "gasProperties file not found" << nl
                         << exit(FatalError);
  }

  Info << "The gas-property table must have 5 columns: p, T, M, h, nu"
       << nl;
  int columnTableGP = 5;
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
    FatalErrorInFunction << "... and it does not seem to have 5 columns." << nl
                         << rawTableIntGP << " lines and "
                         << rawTableFracGP * columnTableGP << " 'scalar' have been read."
                         << nl
                         << exit(FatalError);

  } else {
    rawTableGP = i_rawGP / columnTableGP;
    Info << "The gas-property table has been read (" << rawTableGP << " lines)."
         << nl;
  }

  RectangularMatrix<scalar> TableGP(rawTableGP, columnTableGP);
  for (int x = 0; x < rawTableGP; x++) {
    for (int i = 0; i < columnTableGP; i++) {
      gasPropertyInputFile >> TableGP[x][i];
    }
  }

  // "gasPT" is an object of the class TabulatedGasPropertiesObject, which handles the gas property lookup and interpolations when called
  gasPT_ptr.reset(new TabulatedGasPropertiesObject(TableGP));
  update();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TabulatedGasPropertiesModel::~TabulatedGasPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::TabulatedGasPropertiesModel::update()
{
  if (this->debug_) {
    Info << "--- update gas properties --- Foam::TabulatedGasPropertiesModel::update()" << endl;
  }
  TabulatedGasPropertiesObject& gasPT = gasPT_ptr();
  // Gas thermodynamics and transport properties

  if (this->debug_) {
    Info << "--- M_g --- Foam::TabulatedGasPropertiesModel::update()" << endl;
    Info << "---  min,max p,T " << min(p) << " " << max(p) << " " << min(Tg) << " " << max(Tg) << "--- Foam::TabulatedGasPropertiesModel::update()" << endl;
  }
  // Read data from a table
  forAll (M, i) {
    M[i] = gasPT.M(p[i], Tg[i]);
  }

  if (this->debug_) {
    Info << "--- h_g --- Foam::TabulatedGasPropertiesModel::update()" << endl;
  }
  forAll (h_g, i) {
    h_g[i] = gasPT.h(p[i], Tg[i]);
  }

  if (this->debug_) {
    Info << "--- mu_g --- Foam::TabulatedGasPropertiesModel::update()" << endl;
  }
  forAll (mu, i) {
    mu[i] = gasPT.mu(p[i], Tg[i]);
  }

  if (this->debug_) {
    Info << "--- boundaries --- Foam::TabulatedGasPropertiesModel::update()" << endl;
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
      }
    }
  }
  if (this->debug_) {
    Info << "--- rho_g eps_g --- Foam::TabulatedGasPropertiesModel::update()" << endl;
  }
  const dimensionedScalar R = constant::physicoChemical::R;
  rho_g  = p * M / (R * Tg); // gas density (perfect gas law)
  eps_g = eps_g_c_ + (eps_g_v_ - eps_g_c_) * tau_; // empirical porosity

  if (this->debug_) {
    Info << "--- end --- Foam::TabulatedGasPropertiesModel::update()" << endl;
  }
}

// ************************************************************************* //

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

#include "EquilibriumGasPropertiesModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EquilibriumGasPropertiesModel::EquilibriumGasPropertiesModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleGasPropertiesModel(mesh, dictName),
verifyMaterialChemistry_(verifyMaterialChemistry()),
M(simpleGasPropertiesModel::M_),
mu(simpleGasPropertiesModel::mu_),
h_g(simpleGasPropertiesModel::h_g_),
eps_g(simpleGasPropertiesModel::eps_g_),
rho_g(simpleGasPropertiesModel::rho_g_),
Tg(meshLookupOrConstructScalar(mesh,"Ta")),
p(meshLookupOrConstructScalar(mesh,"p")),
Tg_old
(
    IOobject
    (
        "Ta_old",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    Tg
),
p_old
(
    IOobject
    (
        "p_old",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    p
),
T_threshold_(simpleGasPropertiesModel::materialDict_.subDict("MaterialProperties").template lookupOrDefault<scalar>("T_threshold",10)),
p_threshold_(simpleGasPropertiesModel::materialDict_.subDict("MaterialProperties").template lookupOrDefault<scalar>("p_threshold",10)),
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
tau_(pyrolysisModel_.tau()),
MaterialChemistryModel_(meshLookupOrConstructModel<simpleMaterialChemistryModel>(mesh,dictName,"MaterialChemistry")),
mix_(MaterialChemistryModel_.mixture()),
pTp(new double [2]),
massFractions_(MaterialChemistryModel_.massFractions()),
moleFractions_(MaterialChemistryModel_.moleFractions()),
p_Zx(new double [moleFractions_.size()])
{
  Tg_old == dimensionedScalar("0",dimTemperature,0);
  Tg_old.boundaryFieldRef() == 0;
  p_old == dimensionedScalar("0",dimensionSet(1,-1,-2,0,0,0,0),0);
  p_old.boundaryFieldRef() == 0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::EquilibriumGasPropertiesModel::~EquilibriumGasPropertiesModel()
{
  delete[] pTp;
  delete[] p_Zx;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::EquilibriumGasPropertiesModel::update()
{

  if (simpleGasPropertiesModel::debug_) {
    Info << "\t update M, hg, mu"<< endl;
  }
  forAll(M, cellI) {
    if ((mag(Tg[cellI] - Tg_old[cellI]) > T_threshold_) || (mag(p[cellI] - p_old[cellI]) > p_threshold_)) {
      Tg_old[cellI] = Tg[cellI];
      p_old[cellI] = p[cellI];
      scalar T_var = Tg[cellI];
      scalar p_var = p[cellI];

      forAll(moleFractions_, elemI) {
        p_Zx[elemI] = moleFractions_[elemI][cellI];
      }

      pTp[1] = T_var;
      pTp[0] = p_var;
      mix_().setState(p_Zx, pTp, 2);

      M[cellI]  = mix_().mixtureMw();
      h_g[cellI] = mix_().mixtureHMass();
      mu[cellI] = mix_().viscosity();
    }
  }

  forAll(mesh_.boundaryMesh(), patchI) {
    forAll(M.boundaryField()[patchI], faceI) {

      if ((mag(Tg.boundaryField()[patchI][faceI] - Tg_old.boundaryField()[patchI][faceI]) > T_threshold_)
          || (mag(p.boundaryField()[patchI][faceI] - p_old.boundaryField()[patchI][faceI]) > p_threshold_)) {
        Tg_old.boundaryFieldRef()[patchI][faceI] = Tg.boundaryField()[patchI][faceI];
        p_old.boundaryFieldRef()[patchI][faceI] = p.boundaryField()[patchI][faceI];
        scalar T_var = Tg.boundaryField()[patchI][faceI];
        scalar p_var = p.boundaryField()[patchI][faceI];
        forAll(moleFractions_, elemI) {
          p_Zx[elemI] = moleFractions_[elemI].boundaryField()[patchI][faceI];
        }

        pTp[1] = T_var;
        pTp[0] = p_var;

        mix_().setState(p_Zx, pTp, 2);

        M.boundaryFieldRef()[patchI][faceI]  = mix_().mixtureMw();
        h_g.boundaryFieldRef()[patchI][faceI] = mix_().mixtureHMass();
        mu.boundaryFieldRef()[patchI][faceI] = mix_().viscosity();
      }
    }
  }

  const dimensionedScalar R = constant::physicoChemical::R;
  rho_g  = p * M / (R * Tg); // gas density (perfect gas law)
  eps_g = eps_g_c_ + (eps_g_v_ - eps_g_c_) * tau_; // empirical porosity
}

Switch Foam::EquilibriumGasPropertiesModel::verifyMaterialChemistry()
{
  if (!simpleGasPropertiesModel::materialDict_.isDict("MaterialChemistry")) {
    FatalErrorInFunction << "MaterialChemistry not found in \"constant/" << simpleGasPropertiesModel::materialDict_.path().name() << "/" << simpleGasPropertiesModel::materialDict_.name() << "\"" << exit(FatalError);
  }

  word chemType_(simpleGasPropertiesModel::materialDict_.subDict("MaterialChemistry").lookup("MaterialChemistryType"));
  wordList validChemType_;
  validChemType_.append("ConstantEquilibrium");
  validChemType_.append("EquilibriumElement");
  bool valid_ = false;
  forAll(validChemType_, typeI) {
    if (chemType_ == validChemType_[typeI]) {
      valid_=true;
    }
  }
  if (!valid_) {
    FatalErrorInFunction << "MaterialChemistryType is not compatible with the GasPropertiesType." << nl << "Only MaterialChemistryType valid: " <<  validChemType_ << exit(FatalError);
  }

  return "yes";
}

// ************************************************************************* //

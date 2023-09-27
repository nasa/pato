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
    along with OpenFOAM.  If FiniteRatet, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FiniteRateGasPropertiesModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FiniteRateGasPropertiesModel::FiniteRateGasPropertiesModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleGasPropertiesModel(mesh, regionName),
verifyMaterialChemistry_(verifyMaterialChemistry()),
Ediff(createVolField<scalar>("Ediff",dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0))),
eps_g(createVolField<scalar>("eps_g",dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0))),
M(createVolField<scalar>("M_g",dimensionedScalar("0", dimensionSet(1, 0, 0, 0, -1, 0, 0), 0.02))),
mu(createVolField<scalar>("mu_g",dimensionedScalar("0", dimensionSet(1, -1, -1, 0, 0, 0, 0), 2e-5))),
h_g(createVolField<scalar>("h_g",dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0))),
rho_g(createVolField<scalar>("rho_g",dimensionedScalar("zero", dimensionSet(1, -3, 0, 0, 0, 0, 0), 0.0))),
eps_g_c_(createDimScalarProp("eps_g_c")),
eps_g_v_(createDimScalarProp("eps_g_v")),
eta0_(createDimScalarProp("eta0")),
p_threshold_(createScalarProp("p_threshold","yes",10)),
T_threshold_(createScalarProp("T_threshold","yes",10)),
R(constant::physicoChemical::R),
energyModel_(refModel<simpleEnergyModel>()),
Tg(energyModel_.refVolField<scalar>("Ta")),
massModel_(refModel<simpleMassModel>()),
p(massModel_.refVolField<scalar>("p")),
pyrolysisModel_(refModel<simplePyrolysisModel>()),
tau_(pyrolysisModel_.refVolField<scalar>("tau")),
materialChemistryModel_(refModel<simpleMaterialChemistryModel>()),
mix_(const_cast<autoPtr<Mutation::Mixture>&>(materialChemistryModel_.mixture())),
ns_mix(mix_().nSpecies()),
pTp(new double [2]),
massFractions_(materialChemistryModel_.massFractions()),
speciesNames_(materialChemistryModel_.speciesNames()),
speciesIndexMutation_(materialChemistryModel_.speciesIndexMutation()),
Tg_old(createVolField<scalar>("Ta_old",Tg)),
p_old(createVolField<scalar>("p_old",p)),
initialized_(false)
{
  Dm_.resize(ns_mix);
  forAll(Dm_, specI) {
    Dm_.set(specI,createVolField<scalar>("Dm[" + speciesNames_[specI] + "]",dimensionedScalar("0",dimensionSet(0,2,-1,0,0), 0.0)));
  }
  hSpeciesField_.resize(ns_mix);
  forAll(hSpeciesField_, specI) {
    hSpeciesField_.set(specI,createVolField<scalar>("hSpeciesField[" + std::to_string(specI) +"](" + speciesNames_[specI] + ")",dimensionedScalar("0",dimensionSet(0,2,-2,0,0), 0.0)));
  }
  p_Y = new double[ns_mix]; // species mass fraction
  p_Dm = new double [ns_mix]; // species diffusion coefficients
  p_hSpecies = new double [ns_mix]; // Species enthalpy
  for(int i = 0; i < ns_mix; i++) {
    p_Y[i]=0;
    p_Dm[i]=0;
    p_hSpecies[i]=0;
  }
  update();
  modelInitialized();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FiniteRateGasPropertiesModel::~FiniteRateGasPropertiesModel()
{
  delete[] pTp;
  delete[] p_Y;
  delete[] p_Dm;
  delete[] p_hSpecies;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FiniteRateGasPropertiesModel::update()

{
  if (simpleGasPropertiesModel::debug_) {
    Info << "update M, hg, mu, Dm, hSpecies, Ediff"<< endl;
  }

  if (!initialized_) {
    std::cout.setstate(std::ios_base::failbit);
  }

  forAll(M, cellI) {
    scalar T_var = Tg[cellI];
    scalar p_var = p[cellI];
    forAll(massFractions_, specI) {
      p_Y[speciesIndexMutation_[specI]] = massFractions_[specI][cellI];
    }

    pTp[1] = T_var;
    pTp[0] = p_var;
    mix_->setState(p_Y, pTp, 2);

    M[cellI]  = mix_->mixtureMw();
    h_g[cellI] = mix_->mixtureHMass();
    mu[cellI] = mix_->viscosity();

    mix_->averageDiffusionCoeffs(p_Dm);

    forAll(massFractions_, specI) {
      Dm_[specI][cellI]  = p_Dm[speciesIndexMutation_[specI]];
    }

    mix_->speciesHOverRT(T_var, p_hSpecies);
    forAll(massFractions_, specI) {
      hSpeciesField_[specI][cellI] =
          p_hSpecies[speciesIndexMutation_[specI]] / mix_->speciesMw(speciesIndexMutation_[specI]) *
          R.value() *
          T_var;
    }
  }

  forAll(mesh_.boundaryMesh(), patchI) {
    forAll(M.boundaryField()[patchI], faceI) {
      scalar T_var = Tg.boundaryField()[patchI][faceI];
      scalar p_var = p.boundaryField()[patchI][faceI];
      forAll(massFractions_, specI) {
        p_Y[speciesIndexMutation_[specI]] = massFractions_[specI].boundaryField()[patchI][faceI];
      }

      pTp[1] = T_var;
      pTp[0] = p_var;

      mix_->setState(p_Y, pTp, 2);

      M.boundaryFieldRef()[patchI][faceI]  = mix_->mixtureMw();
      h_g.boundaryFieldRef()[patchI][faceI] = mix_->mixtureHMass();
      mu.boundaryFieldRef()[patchI][faceI] = mix_->viscosity();

      mix_->averageDiffusionCoeffs(p_Dm);
      for (int i = 0; i < ns_mix; i++) {
        Dm_[i].boundaryFieldRef()[patchI][faceI]  = p_Dm[speciesIndexMutation_[i]];
      }
      mix_->speciesHOverRT(T_var, p_hSpecies);
      for (int i = 0; i < ns_mix; i++) {
        hSpeciesField_[i].boundaryFieldRef()[patchI][faceI] =
            p_hSpecies[speciesIndexMutation_[i]] / mix_->speciesMw(speciesIndexMutation_[i]) *
            R.value() *
            T_var;
      }
    }
  }

  if (!initialized_) {
    std::cout.clear();
    initialized_=true;
  }

  rho_g  = p * M / (R * Tg); // gas density (perfect gas law)
  eps_g = eps_g_c_ + (eps_g_v_ - eps_g_c_) * tau_; // empirical porosity

  // Compute the flux transported by diffusion of the molecules
  Ediff = 0*Ediff;
  forAll(massFractions_, specI) {
    if (speciesNames_[specI]!="C(gr)") {
      const volScalarField& Dmi = Dm_[specI];
      const volScalarField& hSpeciesFieldi = hSpeciesField_[specI];
      // Flux transported by molecular diffusion
      Ediff += fvc::laplacian(Dmi / eta0_ * eps_g * rho_g * hSpeciesFieldi, massFractions_[specI]);
    }
  }
}


Switch Foam::FiniteRateGasPropertiesModel::verifyMaterialChemistry()
{
  if (!simpleGasPropertiesModel::materialDict_.isDict("MaterialChemistry")) {
    FatalErrorInFunction << "MaterialChemistry not found in \"constant/" << simpleGasPropertiesModel::materialDict_.path().name() << "/" << simpleGasPropertiesModel::materialDict_.name() << "\"" << exit(FatalError);
  }

  word chemType_(simpleGasPropertiesModel::materialDict_.subDict("MaterialChemistry").lookup("MaterialChemistryType"));
  wordList validChemType_;
  validChemType_.append("SpeciesConservation");
  validChemType_.append("OnlyFiniteRate");

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

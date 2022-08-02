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

#include "basicMutation.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
defineTypeNameAndDebug(basicMutation, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicMutation::basicMutation
(
    const fvMesh& mesh,
    const word& dictName
)
  :
basicUserThermo(mesh, dictName)
{
  IOdictionary thermoDict_
  (
      IOobject
      (
          dictName,
          mesh.time().constant(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE,
          false
      )
  );

  mixtureName_ = thermoDict_.lookupOrDefault<word>("mutationMixture", "air13");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicMutation::~basicMutation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::basicMutation::getElementsComposition(scalarList& Xe, scalarList& Ye, scalarList& rhoe)
{
  const scalar* const pXs   = pMutationMix_->X(); //new scalar[this->nSpecies()];
  scalar* pXe   = new scalar[this->nElements()];
  scalar* pYe   = new scalar[this->nElements()];

  pMutationMix_->convert<Mutation::Thermodynamics::ConversionType::X_TO_XE>(pXs, pXe);
  pMutationMix_->convert<Mutation::Thermodynamics::ConversionType::XE_TO_YE>(pXe, pYe);

  scalar rho;
  this->getRho(rho);

  forAll(rhoe, ie) {
    Xe[ie] = pXe[ie];
    Ye[ie] = pYe[ie];
    rhoe[ie] = rho * Ye[ie];
  }

  // delete pXs;
  delete[] pXe;
  delete[] pYe;
}

// ====================================================== //

void Foam::basicMutation::getSpeciesComposition(scalarList& Xs, scalarList& Ys, scalarList& rhos)
{
  const scalar* const pXs   = pMutationMix_->X(); //new scalar[this->nSpecies()];
  const scalar* const pYs   = pMutationMix_->Y(); //new scalar[this->nSpecies()];
  scalar* pRhos = new scalar[this->nSpecies()];

  pMutationMix_->densities(pRhos);

  forAll(Xs, is) {
    Xs[is]   = pXs[is];
    Ys[is]   = pYs[is];
    rhos[is] = pRhos[is];
  }

  //delete pXs;
  //delete pYs;
  delete[] pRhos;
}

// ====================================================== //

void Foam::basicMutation::getMixtureR(scalar& mixR)
{
  scalar mixtureMM;
  this->getMixtureMolarMass(mixtureMM);

  mixR = this->RU() / mixtureMM;
}

// ====================================================== //

void Foam::basicMutation::getMixtureMolarMass(scalar& mixMM)
{
  mixMM = pMutationMix_->mixtureMw();
}

// ====================================================== //

void Foam::basicMutation::getE(scalar& e)
{
  e = pMutationMix_->mixtureEnergyMass();
}

// ====================================================== //

void Foam::basicMutation::getH(scalar& h)
{
  h = pMutationMix_->mixtureHMass();
}

// ====================================================== //

void Foam::basicMutation::getT(scalar& T)
{
  T = pMutationMix_->T();
}

// ====================================================== //

void Foam::basicMutation::getP(scalar& p)
{
  p = pMutationMix_->P();
}

// ====================================================== //

void Foam::basicMutation::getRho(scalar& rho)
{
  rho = pMutationMix_->density();
}

// ====================================================== //

void Foam::basicMutation::getPsi(scalar& psi)
{
  scalar rho;
  this->getRho(rho);

  scalar p;
  this->getP(p);

  psi = rho / p;
}

// ====================================================== //

void Foam::basicMutation::getMu(scalar& mu)
{
  mu = pMutationMix_->viscosity();
}

// ====================================================== //

void Foam::basicMutation::getNu(scalar& nu)
{
  scalar mu;
  this->getMu(mu);

  scalar rho;
  this->getRho(rho);

  nu = mu / rho;
}

// ====================================================== //

void Foam::basicMutation::getSigma(scalar& sigma)
{
  sigma = pMutationMix_->electricConductivity();
}

// ====================================================== //

void Foam::basicMutation::getThermalDiffRatios(scalarList& tDiffR)
{
  scalar* pThermalDiffRatios = new scalar[this->nSpecies()];

  this->pMutationMix_->thermalDiffusionRatios(pThermalDiffRatios);

  forAll(tDiffR, is) {
    tDiffR[is] = pThermalDiffRatios[is];
  }

  tDiffR[0] = this->pMutationMix_->electronThermalDiffusionRatio();

  delete[] pThermalDiffRatios;
}

// ************************************************************************* //

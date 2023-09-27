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

#include "PTEquilMutation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PTEquilMutation::PTEquilMutation
(
    const fvMesh& mesh,
    const word& dictName
)
  :
basicMutation(mesh, dictName)
{
  stateModelName_ = "Equil";

  configureMixture();
  createCompositionFields(mesh);
  update();
  e_.write();
  isInConstructor_ = 0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PTEquilMutation::~PTEquilMutation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PTEquilMutation::update()
{
  scalarField& TCells = this->T_.primitiveFieldRef();
  scalarField& pCells = this->p_.primitiveFieldRef();
  scalarField& eCells     = this->e_.primitiveFieldRef();
  scalarField& rhoCells   = this->rho_.primitiveFieldRef();

  if (isInConstructor_) {
    forAll(TCells, iCell) {
      //- ATM: assumes that here we set the state knowing p, T (at initialization)
      setState(pCells[iCell], TCells[iCell], 1);

      getRho(rhoCells[iCell]);
      getE(eCells[iCell]);
    }

    volScalarField::Boundary& pBf     = this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf     = this->T_.boundaryFieldRef();
    volScalarField::Boundary& eBf     = this->e_.boundaryFieldRef();
    volScalarField::Boundary& rhoBf   = this->rho_.boundaryFieldRef();

    forAll(TBf, iPatch) {
      fvPatchScalarField& patchP     = pBf[iPatch];
      fvPatchScalarField& patchT     = TBf[iPatch];
      fvPatchScalarField& patchE     = eBf[iPatch];
      fvPatchScalarField& patchRho   = rhoBf[iPatch];

      forAll(patchT, iFace) {
        setState(patchP[iFace], patchT[iFace], 1);

        getRho(patchRho[iFace]);
        getE(patchE[iFace]);
      }
    }
  }

  scalarField& psiCells   = this->psi_.primitiveFieldRef();
  scalarField& muCells    = this->mu_.primitiveFieldRef();
  scalarField& nuCells    = this->nu_.primitiveFieldRef();
  scalarField& alphaCells = this->alpha_.primitiveFieldRef();
  scalarField& sigmaCells = this->sigma_.primitiveFieldRef();
  scalarField& gammaCells = this->gamma_.primitiveFieldRef();
  scalarField& kappaCells = this->kappa_.primitiveFieldRef();

  PtrList<volScalarField>& YeCells   = this->Ye_;
  PtrList<volScalarField>& YsCells   = this->Ys_;
  PtrList<volScalarField>& XeCells   = this->Xe_;
  PtrList<volScalarField>& XsCells   = this->Xs_;
  PtrList<volScalarField>& rhoeCells = this->rhoe_;
  PtrList<volScalarField>& rhosCells = this->rhos_;

  forAll(TCells, iCell) {
    currentRho_ = rhoCells[iCell];
    currentE_   = eCells[iCell];

    //- ATM: assumes that here we set the state knowing rho, e
    setState(currentRho_, currentE_, 0);

    this->updateCurrent();

    forAll(this->currentXe_, ie) {
      XeCells[ie][iCell]   = currentXe_[ie];
      YeCells[ie][iCell]   = currentYe_[ie];
      rhoeCells[ie][iCell] = currentRhoe_[ie];
    }

    forAll(this->currentXs_, is) {
      XsCells[is][iCell]   = currentXs_[is];
      YsCells[is][iCell]   = currentYs_[is];
      rhosCells[is][iCell] = currentRhos_[is];
    }

    TCells[iCell]     = currentT_;
    pCells[iCell]     = currentP_;
    psiCells[iCell]   = currentPsi_;
    muCells[iCell]    = currentMu_;
    nuCells[iCell]    = currentNu_;
    gammaCells[iCell] = currentGamma_;
    kappaCells[iCell] = currentKappa_;
    alphaCells[iCell] = currentAlpha_;
    sigmaCells[iCell] = currentSigma_;
  }

  this->updateBoundaries();
}

// ====================================================== //

void Foam::PTEquilMutation::updateBoundaries()
{
  volScalarField::Boundary& pBf     = this->p_.boundaryFieldRef();
  volScalarField::Boundary& TBf     = this->T_.boundaryFieldRef();
  volScalarField::Boundary& eBf     = this->e_.boundaryFieldRef();
  volScalarField::Boundary& psiBf   = this->psi_.boundaryFieldRef();
  volScalarField::Boundary& rhoBf   = this->rho_.boundaryFieldRef();
  volScalarField::Boundary& muBf    = this->mu_.boundaryFieldRef();
  volScalarField::Boundary& nuBf    = this->nu_.boundaryFieldRef();
  volScalarField::Boundary& alphaBf = this->alpha_.boundaryFieldRef();
  volScalarField::Boundary& gammaBf = this->gamma_.boundaryFieldRef();
  volScalarField::Boundary& sigmaBf = this->sigma_.boundaryFieldRef();
  volScalarField::Boundary& kappaBf = this->kappa_.boundaryFieldRef();

  forAll(this->T_.boundaryField(), iPatch) {
    fvPatchScalarField& patchP     = pBf[iPatch];
    fvPatchScalarField& patchT     = TBf[iPatch];
    fvPatchScalarField& patchE     = eBf[iPatch];
    fvPatchScalarField& patchPsi   = psiBf[iPatch];
    fvPatchScalarField& patchRho   = rhoBf[iPatch];
    fvPatchScalarField& patchMu    = muBf[iPatch];
    fvPatchScalarField& patchNu    = nuBf[iPatch];
    fvPatchScalarField& patchAlpha = alphaBf[iPatch];
    fvPatchScalarField& patchGamma = gammaBf[iPatch];
    fvPatchScalarField& patchKappa = kappaBf[iPatch];
    fvPatchScalarField& patchSigma = sigmaBf[iPatch];

    if (patchT.fixesValue()||isInConstructor_) {
      forAll(patchT, iFace) {
        currentT_ = patchT[iFace];
        currentP_ = patchP[iFace];

        setState(currentP_, currentT_, 1);

        this->updateCurrent();

        forAll(this->currentXe_, ie) {
          this->Xe_[ie].boundaryFieldRef()[iPatch][iFace]   = currentXe_[ie];
          this->Ye_[ie].boundaryFieldRef()[iPatch][iFace]   = currentYe_[ie];
          this->rhoe_[ie].boundaryFieldRef()[iPatch][iFace] = currentRhoe_[ie];
        }

        forAll(this->currentXs_, is) {
          this->Xs_[is].boundaryFieldRef()[iPatch][iFace]   = currentXs_[is];
          this->Ys_[is].boundaryFieldRef()[iPatch][iFace]   = currentYs_[is];
          this->rhos_[is].boundaryFieldRef()[iPatch][iFace] = currentRhos_[is];
        }

        patchT[iFace]     = currentT_;
        patchP[iFace]     = currentP_;
        patchE[iFace]     = currentE_;
        patchRho[iFace]   = currentRho_;
        patchPsi[iFace]   = currentPsi_;
        patchMu[iFace]    = currentMu_;
        patchNu[iFace]    = currentNu_;
        patchGamma[iFace] = currentGamma_;
        patchKappa[iFace] = currentKappa_;
        patchAlpha[iFace] = currentAlpha_;
        patchSigma[iFace] = currentSigma_;
      }
    } else {
      forAll(patchT, iFace) {
        currentE_   = patchE[iFace];
        currentRho_ = patchRho[iFace];

        setState(currentRho_, currentE_, 0);

        this->updateCurrent();

        forAll(this->currentXe_, ie) {
          this->Xe_[ie].boundaryFieldRef()[iPatch][iFace] = currentXe_[ie];
          this->Ye_[ie].boundaryFieldRef()[iPatch][iFace] = currentYe_[ie];
          this->rhoe_[ie].boundaryFieldRef()[iPatch][iFace] = currentRhoe_[ie];
        }

        forAll(this->currentXs_, is) {
          this->Xs_[is].boundaryFieldRef()[iPatch][iFace] = currentXs_[is];
          this->Ys_[is].boundaryFieldRef()[iPatch][iFace] = currentYs_[is];
          this->rhos_[is].boundaryFieldRef()[iPatch][iFace] = currentRhos_[is];
        }

        patchT[iFace]     = currentT_;
        patchP[iFace]     = currentP_;
        patchE[iFace]     = currentE_;
        patchRho[iFace]   = currentRho_;
        patchPsi[iFace]   = currentPsi_;
        patchMu[iFace]    = currentMu_;
        patchNu[iFace]    = currentNu_;
        patchGamma[iFace] = currentGamma_;
        patchKappa[iFace] = currentKappa_;
        patchAlpha[iFace] = currentAlpha_;
        patchSigma[iFace] = currentSigma_;
      }
    }
  }
}

// ====================================================== //

void Foam::PTEquilMutation::updateCurrent()
{
  getElementsComposition(currentXe_, currentYe_, currentRhoe_);
  getSpeciesComposition(currentXs_, currentYs_, currentRhos_);
  getT(currentT_);
  getP(currentP_);
  getE(currentE_);
  getRho(currentRho_);
  getPsi(currentPsi_);
  getMu(currentMu_);
  getNu(currentNu_);
  getGamma(currentGamma_);
  getKappa(currentKappa_);
  getAlpha(currentAlpha_);
  getSigma(currentSigma_);
}

// ====================================================== //

void Foam::PTEquilMutation::setState(const scalar& var1, const scalar& var2, const int vars = 0)
{
  switch (vars) {
  case 0: {
    //- Set the state using rho and e
    //	param @var1 : mass density
    //	param @var2 : specific energy per unit mass

    scalar rho  = var1;
    scalar rhoE = var1 * var2;

    this->pMutationMix_->setState(&rho, &rhoE, 0);

    break;
  }

  case 1: {
    //- Set the state using p and T
    //	param @var1 : pressure
    //	param @var2 : temperature

    scalar p = var1;
    scalar T = var2;

    this->pMutationMix_->setState(&p, &T, 1);

    break;
  }

  default:
    FatalErrorInFunction
        << "Unknown combination of parameters for setState"
        << nl
        << exit(FatalError);
    // ATM: fatal error TO-DO
  }
}

// ====================================================== //

void Foam::PTEquilMutation::getGamma(scalar& gamma)
{
  gamma = this->pMutationMix_->mixtureEquilibriumGamma();
}

// ====================================================== //

void Foam::PTEquilMutation::getKappa(scalar& kappa)
{
  kappa = this->pMutationMix_->equilibriumThermalConductivity();
}

// ====================================================== //

void Foam::PTEquilMutation::getAlpha(scalar& alpha)
{
  scalar kappa;
  this->getKappa(kappa);


  scalar Cv = this->pMutationMix_->mixtureEquilibriumCvMass();
  alpha = kappa / Cv;
}

// ************************************************************************* //

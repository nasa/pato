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

#include "PTXEquilMutation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PTXEquilMutation::PTXEquilMutation
(
    const fvMesh& mesh,
    const word& dictName
)
  :
basicMutation(mesh, dictName)
{
  stateModelName_ = "Equil";

  configureMixture();

  this->diffMassFluxes_.resize(nElements_);

  this->currentFluxesP_.resize(nElements_+1);
  this->currentFluxesT_.resize(nElements_+1);
  this->currentFluxesX_.resize((nElements_+1) * nElements_);

  this->elementsMassDiffFluxesGradP_.resize(nElements_+1);
  this->elementsMassDiffFluxesGradT_.resize(nElements_+1);
  this->elementsMassDiffFluxesGradX_.resize(nElements_ * nElements_);
  this->elementsEnergyDiffFluxesGradX_.resize(nElements_);

  createCompositionFields(mesh);
  createElementFluxFields(mesh);
  update();

  //- ATM: the following is needed because Mutation++ has both models in 'EquilStateModel'
  //- ATM: over-writing number of mass conservation equations for elements
  this->nMassEqns_ = this->nElements_;

  isInConstructor_ = 0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PTXEquilMutation::~PTXEquilMutation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PTXEquilMutation::update()
{
  scalarField& TCells			= this->T_.primitiveFieldRef();
  scalarField& pCells 		= this->p_.primitiveFieldRef();
  scalarField& eCells     = this->e_.primitiveFieldRef();
  scalarField& rhoCells   = this->rho_.primitiveFieldRef();

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
  PtrList<volVectorField>& diffMassFluxesCells = this->diffMassFluxes_;

  vectorField& qDiffCells = this->qDiff_();

  PtrList<volScalarField>& massFluxPCells = this->elementsMassDiffFluxesGradP_;
  PtrList<volScalarField>& massFluxTCells = this->elementsMassDiffFluxesGradT_;
  PtrList<volScalarField>& massFluxXCells = this->elementsMassDiffFluxesGradX_;

  scalarField& energyFluxPCells = this->elementsEnergyDiffFluxesGradP_();
  scalarField& energyFluxTCells = this->elementsEnergyDiffFluxesGradT_();
  PtrList<volScalarField>& energyFluxXCells = this->elementsEnergyDiffFluxesGradX_;

  vectorField gradPCells = (fvc::grad(this->p_));
  vectorField gradTCells = (fvc::grad(this->T_));
  List<vectorField> gradXeCells(this->nElements());

  forAll(gradXeCells, ie) {
    gradXeCells[ie] = (fvc::grad(this->Xe_[ie]));
  }

  scalarList tmpRhoe(this->nElements());
  scalarList tmpMassFluxP(this->nElements());
  scalarList tmpMassFluxT(this->nElements());
  scalarList tmpMassFluxXe((this->nElements())*(this->nElements()));
  scalarList tmpEnergyFluxXe(this->nElements());
  List<vector> tmpGradXe(this->nElements());
  List<vector> tmpDiffMassFluxes(this->nElements());


  if (isInConstructor_) {
    scalarList tmpYe(this->nElements());

    forAll(TCells, iCell) {
      forAll(tmpYe, ie) {
        tmpYe[ie] = YeCells[ie][iCell];
      }

      //- ATM: assumes that here we set the state knowing p, T, Ye (at initialization)
      setState(pCells[iCell], TCells[iCell], tmpYe, 1);

      getRho(rhoCells[iCell]);
      getE(eCells[iCell]);

      forAll(tmpYe, ie) {
        rhoeCells[ie][iCell] = tmpYe[ie] * rhoCells[iCell];
      }
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
        forAll(tmpYe,ie) {
          tmpYe[ie] = this->Ye_[ie].boundaryFieldRef()[iPatch][iFace];
        }

        setState(patchP[iFace], patchT[iFace], tmpYe, 1);

        getRho(patchRho[iFace]);
        getE(patchE[iFace]);

        forAll(tmpYe, ie) {
          this->rhoe_[ie].boundaryFieldRef()[iPatch][iFace] = tmpYe[ie] * patchRho[iFace];
        }
      }
    }
  }

  forAll(TCells, iCell) {
    currentRho_ = rhoCells[iCell];
    currentE_   = eCells[iCell];
    forAll(this->currentYe_, ie) {
      currentYe_[ie] = rhoeCells[ie][iCell] / rhoCells[iCell];
    }

    //- ATM: assumes that here we set the state knowing rho, e, rhoe (after initialization)
    setState(currentRho_, currentE_, currentYe_, 0);

    this->updateCurrent();

    forAll(this->currentXe_, ie) {
      XeCells[ie][iCell]   = currentXe_[ie];
      YeCells[ie][iCell]   = currentYe_[ie];

      massFluxPCells[ie][iCell] = currentFluxesP_[ie];
      massFluxTCells[ie][iCell] = currentFluxesT_[ie];

      energyFluxXCells[ie][iCell] = currentFluxesX_[(ie+1)*this->nElements() - 1];

      forAll(this->currentXe_, ie2) {
        massFluxXCells[ie*this->nElements() + ie2][iCell] = currentFluxesX_[ie*(this->nElements()+1) + ie2];
      }

      tmpRhoe[ie] = rhoeCells[ie][iCell];
      tmpGradXe[ie] = gradXeCells[ie][iCell];
      tmpMassFluxP[ie] = massFluxPCells[ie][iCell];
      tmpMassFluxT[ie] = massFluxTCells[ie][iCell];
      tmpEnergyFluxXe[ie] = energyFluxXCells[ie][iCell];

      forAll(this->currentXe_, ie2) {
        tmpMassFluxXe[ie*this->nElements() + ie2] = massFluxXCells[ie*this->nElements() + ie2][iCell];
      }
    }

    energyFluxPCells[iCell] = currentFluxesP_[this->nElements()];
    energyFluxTCells[iCell] = currentFluxesT_[this->nElements()];

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

    getDiffMassFluxes(tmpMassFluxP, tmpMassFluxT, tmpMassFluxXe, gradPCells[iCell], gradTCells[iCell], tmpGradXe, tmpDiffMassFluxes);

    forAll(tmpDiffMassFluxes, ie) {
      diffMassFluxesCells[ie][iCell] = tmpDiffMassFluxes[ie];
    }

    getDiffHeatFlux(energyFluxPCells[iCell], energyFluxTCells[iCell], tmpEnergyFluxXe, gradPCells[iCell], gradTCells[iCell], tmpGradXe, qDiffCells[iCell]);
  }

  this->updateBoundaries();
}

// ====================================================== //

void Foam::PTXEquilMutation::updateBoundaries()
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

  scalarList tmpRhoe(this->nElements());
  scalarList tmpMassFluxP(this->nElements());
  scalarList tmpMassFluxT(this->nElements());
  scalarList tmpMassFluxXe((this->nElements())*(this->nElements()));
  scalarList tmpEnergyFluxXe(this->nElements());
  List<vector> tmpGradXe(this->nElements());
  List<vector> tmpDiffMassFluxes(this->nElements());

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

        forAll(this->currentYe_,ie) {
          currentYe_[ie] = this->Ye_[ie].boundaryFieldRef()[iPatch][iFace];
        }

        setState(currentP_, currentT_, currentYe_, 1);

        this->updateCurrent();

        forAll(this->currentXe_, ie) {
          this->Xe_[ie].boundaryFieldRef()[iPatch][iFace]   = currentXe_[ie];
          this->Ye_[ie].boundaryFieldRef()[iPatch][iFace]   = currentYe_[ie];

          tmpMassFluxP[ie] = this->elementsMassDiffFluxesGradP_[ie].boundaryFieldRef()[iPatch][iFace];
          tmpMassFluxT[ie] = this->elementsMassDiffFluxesGradT_[ie].boundaryFieldRef()[iPatch][iFace];
          this->elementsMassDiffFluxesGradP_[ie].boundaryFieldRef()[iPatch][iFace] = currentFluxesP_[ie];
          this->elementsMassDiffFluxesGradT_[ie].boundaryFieldRef()[iPatch][iFace] = currentFluxesT_[ie];

          this->elementsEnergyDiffFluxesGradX_[ie].boundaryFieldRef()[iPatch][iFace] = currentFluxesX_[(ie+1)*this->nElements() - 1];

          forAll(this->currentXe_, ie2) {
            this->elementsMassDiffFluxesGradX_[ie*this->nElements() + ie2].boundaryFieldRef()[iPatch][iFace] = currentFluxesX_[ie*(this->nElements()+1) + ie2];

            tmpMassFluxXe[ie*this->nElements() + ie2] = this->elementsMassDiffFluxesGradX_[ie*this->nElements() + ie2].boundaryFieldRef()[iPatch][iFace];
          }

          tmpRhoe[ie] = this->rhoe_[ie].boundaryFieldRef()[iPatch][iFace];
          tmpEnergyFluxXe[ie] = this->elementsEnergyDiffFluxesGradX_[ie].boundaryFieldRef()[iPatch][iFace];
        }

        this->elementsEnergyDiffFluxesGradP_->boundaryFieldRef()[iPatch][iFace] = currentFluxesP_[this->nElements()];
        this->elementsEnergyDiffFluxesGradT_->boundaryFieldRef()[iPatch][iFace] = currentFluxesT_[this->nElements()];

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

        vector tmpGradP = (fvc::grad(this->p_))->boundaryField()[iPatch][iFace];
        vector tmpGradT = (fvc::grad(this->T_))->boundaryField()[iPatch][iFace];

        getDiffMassFluxes(tmpMassFluxP, tmpMassFluxT, tmpMassFluxXe, tmpGradP, tmpGradT, tmpGradXe, tmpDiffMassFluxes);

        forAll(tmpDiffMassFluxes, ie) {
          this->diffMassFluxes_[ie].boundaryFieldRef()[iPatch][iFace] = tmpDiffMassFluxes[ie];
        }

        getDiffHeatFlux(this->elementsEnergyDiffFluxesGradP_->boundaryFieldRef()[iPatch][iFace], this->elementsEnergyDiffFluxesGradT_->boundaryFieldRef()[iPatch][iFace], tmpEnergyFluxXe, tmpGradP, tmpGradT, tmpGradXe, this->qDiff_->boundaryFieldRef()[iPatch][iFace]);
      }
    } else {
      forAll(patchT, iFace) {
        forAll(this->currentYe_,ie) {
          currentYe_[ie] = this->Ye_[ie].boundaryFieldRef()[iPatch][iFace];
        }

        currentE_   = patchE[iFace];
        currentRho_ = patchRho[iFace];

        setState(currentRho_, currentE_, currentYe_, 0);

        this->updateCurrent();

        forAll(this->currentXe_, ie) {
          this->Xe_[ie].boundaryFieldRef()[iPatch][iFace]   = currentXe_[ie];
          this->Ye_[ie].boundaryFieldRef()[iPatch][iFace]   = currentYe_[ie];

          tmpMassFluxP[ie] = this->elementsMassDiffFluxesGradP_[ie].boundaryFieldRef()[iPatch][iFace];
          tmpMassFluxT[ie] = this->elementsMassDiffFluxesGradT_[ie].boundaryFieldRef()[iPatch][iFace];
          this->elementsMassDiffFluxesGradP_[ie].boundaryFieldRef()[iPatch][iFace] = currentFluxesP_[ie];
          this->elementsMassDiffFluxesGradT_[ie].boundaryFieldRef()[iPatch][iFace] = currentFluxesT_[ie];

          this->elementsEnergyDiffFluxesGradX_[ie].boundaryFieldRef()[iPatch][iFace] = currentFluxesX_[(ie+1)*this->nElements() - 1];

          forAll(this->currentXe_, ie2) {
            this->elementsMassDiffFluxesGradX_[ie*this->nElements() + ie2].boundaryFieldRef()[iPatch][iFace] = currentFluxesX_[ie*(this->nElements()+1) + ie2];

            tmpMassFluxXe[ie*this->nElements() + ie2] = this->elementsMassDiffFluxesGradX_[ie*this->nElements() + ie2].boundaryFieldRef()[iPatch][iFace];
          }

          tmpRhoe[ie] = this->rhoe_[ie].boundaryFieldRef()[iPatch][iFace];
          tmpEnergyFluxXe[ie] = this->elementsEnergyDiffFluxesGradX_[ie].boundaryFieldRef()[iPatch][iFace];
        }

        this->elementsEnergyDiffFluxesGradP_->boundaryFieldRef()[iPatch][iFace] = currentFluxesP_[this->nElements()];
        this->elementsEnergyDiffFluxesGradT_->boundaryFieldRef()[iPatch][iFace] = currentFluxesT_[this->nElements()];

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

        vector tmpGradP = (fvc::grad(this->p_))->boundaryField()[iPatch][iFace];
        vector tmpGradT = (fvc::grad(this->T_))->boundaryField()[iPatch][iFace];

        getDiffMassFluxes(tmpMassFluxP, tmpMassFluxT, tmpMassFluxXe, tmpGradP, tmpGradT, tmpGradXe, tmpDiffMassFluxes);

        forAll(tmpDiffMassFluxes, ie) {
          this->diffMassFluxes_[ie].boundaryFieldRef()[iPatch][iFace] = tmpDiffMassFluxes[ie];
        }

        getDiffHeatFlux(this->elementsEnergyDiffFluxesGradP_->boundaryFieldRef()[iPatch][iFace], this->elementsEnergyDiffFluxesGradT_->boundaryFieldRef()[iPatch][iFace], tmpEnergyFluxXe, tmpGradP, tmpGradT, tmpGradXe, this->qDiff_->boundaryFieldRef()[iPatch][iFace]);
      }
    }
  }
}

// ====================================================== //

void Foam::PTXEquilMutation::updateCurrent()
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
  getElementFluxes(currentFluxesP_, currentFluxesT_, currentFluxesX_);
}

// ====================================================== //

void Foam::PTXEquilMutation::setState(const scalar& var1, const scalar& var2, const scalarList& Ye, const int vars = 0)
{
  scalar* pVars = new scalar[2];
  scalar* pXe   = new scalar[this->nElements()];
  scalar* pYe   = new scalar[this->nElements()];

  forAll(Ye, ie) {
    pYe[ie] = Ye[ie];
  }

  this->pMutationMix_->convert<Mutation::Thermodynamics::ConversionType::YE_TO_XE>(pYe, pXe);

  switch (vars) {
  case 0: {
    //- Set the state using rho and e and Xe

    int iter = 0;
    scalar tol = 1e-10;
    scalar maxIter = 1e6;
    scalar p(currentP_);
    scalar T(currentT_);

    this->setState(p, T, Ye, 1);
    this->updateCurrent();

    scalar Rgas;
    getMixtureR(Rgas);

    scalar Cv = Rgas / (currentGamma_ - 1);
    scalar rho  = var1;
    scalar e    = var2;
    scalar rhoE = rho * e;

    scalar deltaRhoE = rho * (e - currentE_);
    scalar deltaRho  = rho - currentRho_;
    scalar dRhoEdT = rho * Cv; 										//- ATM: this is not accurate, but probably works
    scalar dRhodT  = - rho * Cv / e;							//- ATM: this is not accurate, but probably works
    scalar dRhoEdP = 1 / (currentGamma_ - 1);			//- ATM: this is not accurate, but probably works
    scalar dRhodP  = 1 / (Rgas * (T - e / Cv));		//- ATM: this is not accurate, but probably works
    scalar deltaT = 0.;
    scalar deltaP = 0.;

    while (std::max(std::abs(deltaRhoE / rhoE), std::abs(deltaRho / rho)) > tol) {
      deltaT = deltaRhoE / dRhoEdT + deltaRho / dRhodT;
      deltaP = deltaRhoE / dRhoEdP + deltaRho / dRhodP;

      while (T + deltaT < 0)	{
        deltaT *= 0.5;
      }
      while (p + deltaP < 0)	{
        deltaP *= 0.5;
      }

      T += deltaT;
      p += deltaP;

      this->setState(p, T, Ye, 1);
      this->updateCurrent();

      getMixtureR(Rgas);
      Cv = Rgas / (currentGamma_ - 1);

      deltaRhoE = rho * (e - currentE_);
      deltaRho  = rho - currentRho_;     //rho - (p / (Rgas * T));
      dRhoEdT = rho * Cv;
      dRhodT  = - rho * Cv / e;
      dRhoEdP = 1 / (currentGamma_ - 1);
      dRhodP  = 1 / (Rgas * T - Rgas * e / Cv);

      if (iter++ > maxIter) {
        FatalErrorInFunction
            << "Maximum number of iterations exceeded in TE function (computing p and T for given rho and e)."
            << abort(FatalError);
      }
    }

    this->currentP_ = p;
    this->currentT_ = T;

    pVars[0] = currentP_;

    pVars[1] = currentT_;

    this->pMutationMix_->setState(pXe, pVars, 2);

    break;
  }

  case 1: {
    //- Set the state using p and T and Ye

    //- var1 = p
    pVars[0] = var1;

    //- var2 = T
    pVars[1] = var2;

    this->pMutationMix_->setState(pXe, pVars, 2);

    break;
  }

  default:
    FatalErrorInFunction
        << "Unknown combination of parameters for setState"
        << nl
        << exit(FatalError);
    // ATM: fatal error TO-DO
  }

  delete[] pVars;
  delete[] pXe;
  delete[] pYe;
}

// ====================================================== //

void Foam::PTXEquilMutation::getGamma(scalar& gamma)
{
  gamma = this->pMutationMix_->mixtureEquilibriumGamma();
}

// ====================================================== //

void Foam::PTXEquilMutation::getKappa(scalar& kappa)
{
  kappa = this->pMutationMix_->equilibriumThermalConductivity();

  //- ATM: this may be wrong -- includes Soret conductivity TO-DO FOR-JM
  //			 Energy diffusion due to gradients is added in the
  //			 source terms, too. Need to check that it is consistent.
  //			 May need to remove Soret conductivity if we add all terms.
}

// ====================================================== //

void Foam::PTXEquilMutation::getAlpha(scalar& alpha)
{
  scalar kappa;
  this->getKappa(kappa);

  scalar Cv = this->pMutationMix_->mixtureEquilibriumCvMass();

  alpha = kappa / Cv;
}

// ====================================================== //

void Foam::PTXEquilMutation::getElementFluxes(scalarList& fluxP, scalarList& fluxT, scalarList& fluxX)
{
  int size = this->nElements_+1;

  scalar* pFluxP = new scalar[size];
  scalar* pFluxT = new scalar[size];
  scalar* pFluxX = new scalar[size * this->nElements_];

  pMutationMix_->equilDiffFluxFacsP(pFluxP);
  pMutationMix_->equilDiffFluxFacsT(pFluxT);
  pMutationMix_->equilDiffFluxFacsZ(pFluxX);

  forAll(fluxP, iF) {
    fluxP[iF] = pFluxP[iF];
    fluxT[iF] = pFluxT[iF];
  }

  forAll(fluxX, iFX) {
    fluxX[iFX] = pFluxX[iFX];
  }

  delete[] pFluxP;
  delete[] pFluxT;
  delete[] pFluxX;
}

// ====================================================== //

void Foam::PTXEquilMutation::getDiffMassFluxes(const scalarList& fluxP, const scalarList& fluxT, const scalarList& fluxX, const vector& gradP, const vector& gradT, const List<vector>& gradXe, List<vector>& diffMF)
{
  forAll(fluxP, ie) {
    diffMF[ie] = fluxP[ie] * gradP + fluxT[ie] * gradT;
    forAll(fluxP, ie2) {
      diffMF[ie] += fluxX[ie*this->nElements() + ie2] * gradXe[ie2]; 			//- ATM: check the order of the terms TO-DO
    }
  }
}

// ====================================================== //

void Foam::PTXEquilMutation::getDiffHeatFlux(const scalar& eFluxP, const scalar& eFluxT, const scalarList& eFluxXe, vector& gradP, vector& gradT, List<vector>& gradXe, vector& qDiff)
{
  qDiff = eFluxP * gradP + eFluxT * gradT;
  forAll(eFluxXe, ie) {
    qDiff += eFluxXe[ie] * gradXe[ie];
  }
}

// ************************************************************************* //

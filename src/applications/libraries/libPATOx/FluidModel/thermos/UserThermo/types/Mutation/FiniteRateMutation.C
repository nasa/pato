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

#include "FiniteRateMutation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FiniteRateMutation::FiniteRateMutation
(
    const fvMesh& mesh,
    const word& dictName
)
  :
basicMutation(mesh, dictName)
{
  stateModelName_ = "ChemNonEq1T";

  configureMixture();

  this->currentOmegas_.resize(nSpecies_);

  this->diffVs_.resize(nSpecies_);
  this->omegas_.resize(nSpecies_);

  createCompositionFields(mesh);
  createFiniteRateFields(mesh);
  update();

  isInConstructor_ = 0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FiniteRateMutation::~FiniteRateMutation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FiniteRateMutation::update()
{
  scalarField& TCells 		= this->T_.primitiveFieldRef();
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
  PtrList<volScalarField>& omegasCells = this->omegas_;
  PtrList<volVectorField>& diffVsCells = this->diffVs_;

  vectorField& qDiffCells = this->qDiff_();

  vectorField gradPCells = (fvc::grad(this->p_));
  vectorField gradTCells = (fvc::grad(this->T_));
  List<vectorField> gradXsCells(this->nSpecies());

  forAll(gradXsCells, is) {
    gradXsCells[is] = (fvc::grad(this->Xs_[is]));
  }

  scalarList tmpRhos(this->nSpecies());
  List<vector> tmpDiffVs(this->nSpecies());
  List<vector> tmpGradXs(this->nSpecies());


  if (isInConstructor_) {
    scalarList tmpYs(this->nSpecies());

    forAll(TCells, iCell) {
      forAll(tmpYs, is) {
        tmpYs[is] = YsCells[is][iCell];
      }

      //- ATM: assumes that here we set the state knowing p, T (at initialization)
      setState(pCells[iCell], TCells[iCell], tmpYs, 1);

      getRho(rhoCells[iCell]);
      getE(eCells[iCell]);

      forAll(tmpYs, is) {
        rhosCells[is][iCell] = tmpYs[is] * rhoCells[iCell];
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
        forAll(tmpYs, is) {
          tmpYs[is] = this->Ys_[is].boundaryFieldRef()[iPatch][iFace];
        }

        setState(patchP[iFace], patchT[iFace], tmpYs, 1);

        getRho(patchRho[iFace]);
        getE(patchE[iFace]);

        forAll(tmpYs, is) {
          this->rhos_[is].boundaryFieldRef()[iPatch][iFace] = tmpYs[is] * patchRho[iFace];
        }
      }
    }
  }

  forAll(TCells, iCell) {
    currentRho_ = rhoCells[iCell];
    currentE_ 	= eCells[iCell];

    forAll(this->currentYs_, is) {
      currentYs_[is] = rhosCells[is][iCell] / rhoCells[iCell];
    }

    //- ATM: assumes that here we set the state knowing p, T, Ys (at initialization)
    setState(currentRho_, currentE_, currentYs_, 0);

    this->updateCurrent();

    forAll(this->currentXe_, ie) {
      XeCells[ie][iCell]   = currentXe_[ie];
      YeCells[ie][iCell]   = currentYe_[ie];
      rhoeCells[ie][iCell] = currentRhoe_[ie];
    }

    forAll(this->currentXs_, is) {
      XsCells[is][iCell]   = currentXs_[is];
      YsCells[is][iCell]   = currentYs_[is];
    }

    TCells[iCell]     = currentT_;
    pCells[iCell]   	= currentP_;
    psiCells[iCell]   = currentPsi_;
    muCells[iCell]    = currentMu_;
    nuCells[iCell]    = currentNu_;
    gammaCells[iCell] = currentGamma_;
    kappaCells[iCell] = currentKappa_;
    alphaCells[iCell] = currentAlpha_;
    sigmaCells[iCell] = currentSigma_;

    forAll(this->omegas_, is) {
      omegasCells[is][iCell] = currentOmegas_[is];

      tmpRhos[is] = rhosCells[is][iCell];
      tmpGradXs[is] = gradXsCells[is][iCell];
    }

    getDiffVelocities(gradPCells[iCell], gradTCells[iCell], tmpGradXs, tmpDiffVs);

    forAll(tmpDiffVs, is) {
      diffVsCells[is][iCell] = tmpDiffVs[is];
    }

    getDiffHeatFlux(tmpRhos, tmpDiffVs, qDiffCells[iCell]);
  }

  this->updateBoundaries();
}

// ====================================================== //

void Foam::FiniteRateMutation::updateBoundaries()
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

  scalarList tmpRhos(this->nSpecies());
  List<vector> tmpDiffVs(this->nSpecies());
  List<vector> tmpGradXs(this->nSpecies());

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

    if (patchT.fixesValue()) {
      forAll(patchT, iFace) {
        currentT_ = patchT[iFace];
        currentP_ = patchP[iFace];

        forAll(this->currentYe_,ie) {
          currentYe_[ie] = this->Ye_[ie].boundaryFieldRef()[iPatch][iFace];
        }

        setState(currentP_, currentT_, currentYs_, 1);

        this->updateCurrent();

        forAll(this->currentXe_, ie) {
          this->Xe_[ie].boundaryFieldRef()[iPatch][iFace]   = currentXe_[ie];
          this->Ye_[ie].boundaryFieldRef()[iPatch][iFace]   = currentYe_[ie];
          this->rhoe_[ie].boundaryFieldRef()[iPatch][iFace] = currentRhoe_[ie];
        }

        forAll(this->currentXs_, is) {
          this->Xs_[is].boundaryFieldRef()[iPatch][iFace]   = currentXs_[is];
          this->Ys_[is].boundaryFieldRef()[iPatch][iFace]   = currentYs_[is];
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

        forAll(this->currentOmegas_, is) {
          this->omegas_[is].boundaryFieldRef()[iPatch][iFace] = currentOmegas_[is];

          tmpRhos[is] = this->rhos_[is].boundaryFieldRef()[iPatch][iFace];
          tmpGradXs[is] = (fvc::grad(this->Xs_[is]))->boundaryField()[iPatch][iFace];
        }

        getDiffVelocities(tmpGradP, tmpGradT, tmpGradXs, tmpDiffVs);
        forAll(tmpDiffVs, is) {
          this->diffVs_[is].boundaryFieldRef()[iPatch][iFace] = tmpDiffVs[is];
        }

        getDiffHeatFlux(tmpRhos, tmpDiffVs, this->qDiff_->boundaryFieldRef()[iPatch][iFace]);
      }
    } else {
      forAll(patchT, iFace) {
        forAll(this->currentYe_,ie) {
          currentYe_[ie] = this->Ye_[ie].boundaryFieldRef()[iPatch][iFace];
        }

        currentE_   = patchE[iFace];
        currentRho_ = patchRho[iFace];

        setState(currentRho_, currentE_, currentYs_, 0);

        this->updateCurrent();

        forAll(this->currentXe_, ie) {
          this->Xe_[ie].boundaryFieldRef()[iPatch][iFace]   = currentXe_[ie];
          this->Ye_[ie].boundaryFieldRef()[iPatch][iFace]   = currentYe_[ie];
          this->rhoe_[ie].boundaryFieldRef()[iPatch][iFace] = currentRhoe_[ie];
        }

        forAll(this->currentXs_, is) {
          this->Xs_[is].boundaryFieldRef()[iPatch][iFace]   = currentXs_[is];
          this->Ys_[is].boundaryFieldRef()[iPatch][iFace]   = currentYs_[is];
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
        List<vector> tmpGradXs(this->nSpecies());

        forAll(this->currentOmegas_, is) {
          this->omegas_[is].boundaryFieldRef()[iPatch][iFace] = currentOmegas_[is];

          tmpRhos[is] = this->rhos_[is].boundaryFieldRef()[iPatch][iFace];
          tmpGradXs[is] = (fvc::grad(this->Xs_[is]))->boundaryField()[iPatch][iFace];
        }

        getDiffVelocities(tmpGradP, tmpGradT, tmpGradXs, tmpDiffVs);
        forAll(tmpDiffVs, is) {
          this->diffVs_[is].boundaryFieldRef()[iPatch][iFace] = tmpDiffVs[is];
        }

        getDiffHeatFlux(tmpRhos, tmpDiffVs, this->qDiff_->boundaryFieldRef()[iPatch][iFace]);
      }
    }
  }
}

// ====================================================== //

void Foam::FiniteRateMutation::updateCurrent()
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
  getOmegasMass(currentOmegas_);
}

// ====================================================== //

void Foam::FiniteRateMutation::setState(const scalar& var1, const scalar& var2, const scalarList& Ys, const int vars = 0)
{
  scalar* pSpecies = new scalar[this->nSpecies()];

  forAll(Ys, is) {
    pSpecies[is] = Ys[is];
  }

  switch (vars) {
  case 0: {
    //- Set the state using rhos and e

    scalar* rhoE = new scalar[1];

    rhoE[0] = var1 * var2;

    forAll(Ys, is) {
      pSpecies[is] *= var1;
    }

    this->pMutationMix_->setState(pSpecies, rhoE, 0);

    delete[] rhoE;
    break;
  }

  case 1: {
    //- Set the state using p and T and Ys

    scalar* pVars    = new scalar[2];

    pVars[0] = var1;
    pVars[1] = var2;

    this->pMutationMix_->setState(pSpecies, pVars, 2);

    delete[] pVars;
    break;
  }

  default:
    FatalErrorInFunction
        << "Unknown combination of parameters for setState"
        << nl
        << exit(FatalError);
    // ATM: fatal error TO-DO
  }

  delete[] pSpecies;
}

// ====================================================== //

void Foam::FiniteRateMutation::getGamma(scalar& gamma)
{
  gamma = this->pMutationMix_->mixtureFrozenGamma();
}

// ====================================================== //

void Foam::FiniteRateMutation::getKappa(scalar& kappa)
{
  kappa = this->pMutationMix_->frozenThermalConductivity();
}

// ====================================================== //

void Foam::FiniteRateMutation::getAlpha(scalar& alpha)
{
  scalar kappa;
  this->getKappa(kappa);

  scalar Cv = this->pMutationMix_->mixtureFrozenCvMass();

  alpha = kappa / Cv;
}

// ====================================================== //

void Foam::FiniteRateMutation::getOmegasMass(scalarList& omegas)
{
  scalar* pOmegas = new scalar[this->nSpecies()];

  this->pMutationMix_->netProductionRates(pOmegas);

  forAll(omegas, is) {
    omegas[is] = pOmegas[is];
  }

  delete[] pOmegas;
}

// ====================================================== //

void Foam::FiniteRateMutation::getDiffDrivingForces(vector& gradP, vector& gradT, List<vector>& gradXs, List<vector>& dFs)
{
  scalar p(0.);
  scalar T(0.);
  scalarList Xs(this->nSpecies());
  scalarList Ys(this->nSpecies());
  scalarList rhos(this->nSpecies());
  scalarList thermalDiffusionRatios(this->nSpecies());

  this->getP(p);
  this->getT(T);
  this->getSpeciesComposition(Xs, Ys, rhos);
  this->getThermalDiffRatios(thermalDiffusionRatios);

  forAll(dFs, is) {
    dFs[is] = gradXs[is] + (Xs[is] - Ys[is]) * (gradP / p) + (thermalDiffusionRatios[is]) * (gradT / T);
  }
}

// ====================================================== //

void Foam::FiniteRateMutation::getDiffVelocities(vector& gradP, vector& gradT, List<vector>& gradXs, List<vector>& diffVs)
{
  scalar* pDrivingForces = new scalar[this->nSpecies()];
  scalar* pDiffVs = new scalar[this->nSpecies()];

  List<vector> drivingForces(this->nSpecies());

  getDiffDrivingForces(gradP, gradT, gradXs, drivingForces);

  vector fieldE(0., 0., 0.);

  for(int iComp = 0; iComp < 3; iComp++) {
    forAll(drivingForces, is) {
      pDrivingForces[is] = drivingForces[is].component(iComp);
    }

    this->pMutationMix_->stefanMaxwell(pDrivingForces, pDiffVs, fieldE.component(iComp));

    forAll(diffVs, is) {
      diffVs[is].component(iComp) = pDiffVs[is];
    }
  }

  delete[] pDrivingForces;
  delete[] pDiffVs;
}

// ====================================================== //

void Foam::FiniteRateMutation::getDiffHeatFlux(const scalarList& rhos, const List<vector>& diffVs, vector& qDiff)
{
  scalar T(0.);
  scalar rho(0.);
  scalar Rgas(0.);
  scalarList Hs(this->nSpecies());
  scalarList thermalDiffusionRatios(this->nSpecies());

  this->getT(T);
  this->getRho(rho);
  this->getMixtureR(Rgas);
  this->getThermalDiffRatios(thermalDiffusionRatios);

  scalar* pHs = new scalar[this->nSpecies()];

  this->pMutationMix_->getEnthalpiesMass(pHs);

  qDiff = vector(0.,0.,0.);

  forAll(rhos, is) {
    Hs[is] = pHs[is];
    qDiff += rhos[is]*Hs[is]*diffVs[is] + rho * Rgas * T * thermalDiffusionRatios[is] * diffVs[is];
  }

  delete[] pHs;
}

// ************************************************************************* //

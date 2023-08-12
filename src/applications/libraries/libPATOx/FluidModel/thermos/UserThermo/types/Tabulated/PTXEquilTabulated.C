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

#include "PTXEquilTabulated.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PTXEquilTabulated::PTXEquilTabulated
(
    const fvMesh& mesh,
    const word& dictName
)
  :
basicTabulated(mesh, dictName)
{
  stateModelName_ = "Equil";

  configureMixture();

  this->currentFluxesP_.resize(nElements_+1);
  this->currentFluxesT_.resize(nElements_+1);
  this->currentFluxesX_.resize((nElements_+1) * nElements_);

  this->elementsMassDiffFluxesGradP_.resize(nElements_+1);
  this->elementsMassDiffFluxesGradT_.resize(nElements_+1);
  this->elementsMassDiffFluxesGradX_.resize(nElements_ * nElements_);
  this->elementsEnergyDiffFluxesGradX_.resize(nElements_);

  createCompositionFields(mesh);
  createElementFluxFields(mesh);

  indexAr_=-1;
  forAll(elementsList_, eI) {
    if (elementsList_[eI]=="Ar") {
      indexAr_=eI;
    }
  }

  if(indexAr_<0) {
    FatalErrorInFunction << "Ar element not found in the mixture." << exit(FatalError);
  }

  if(!isFile(basicTabulated::tableFileName_)) {
    generateTable(basicTabulated::tableName_);
  } else {
    readTable(basicTabulated::tableName_);
  }
  update();

  //- ATM: the following is needed because Mutation++ has both models in 'EquilStateModel'
  //- ATM: over-writing number of mass conservation equations for elements
  this->nMassEqns_ = this->nElements_;

  isInConstructor_ = 0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PTXEquilTabulated::~PTXEquilTabulated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PTXEquilTabulated::update()
{
  scalarField& TCells = this->T_.primitiveFieldRef();
  scalarField& pCells = this->p_.primitiveFieldRef();
  scalarField& eCells = this->e_.primitiveFieldRef();
  scalarField& rhoCells = this->rho_.primitiveFieldRef();

  scalarField& psiCells = this->psi_.primitiveFieldRef();
  scalarField& sigmaCells = this->sigma_.primitiveFieldRef();
  scalarField& alphaCells = this->alpha_.primitiveFieldRef();
  scalarField& gammaCells = this->gamma_.primitiveFieldRef();
  scalarField& muCells = this->mu_.primitiveFieldRef();
  scalarField& nuCells = this->nu_.primitiveFieldRef();

  PtrList<volVectorField>& diffMassFluxesCells = this->diffMassFluxes_;
  vectorField& qDiffCells = this->qDiff_();

  PtrList<volScalarField>& YeCells = this->Ye_;
  PtrList<volScalarField>& rhoeCells = this->rhoe_;

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

  forAll(TCells, iCell) {
    scalarList Y(nElements_);
    forAll(Y, eI) {
      Y[eI]=YeCells[eI][iCell];
    }

    if (isInConstructor_) {
      updateCurrentPTY(pCells[iCell], TCells[iCell],Y);
      eCells[iCell] = currentE_;
      rhoCells[iCell] = currentRho_;

      forAll(currentRhoe_, ie) {
        rhoeCells[ie]=currentRhoe_[ie];
      }
    } else {
      updateCurrentRhoEY(rhoCells[iCell], eCells[iCell],Y);
      TCells[iCell] = currentT_;
      pCells[iCell] = currentP_;
    }

    psiCells[iCell] = currentPsi_;
    sigmaCells[iCell] = currentSigma_;
    alphaCells[iCell] = currentAlpha_;
    gammaCells[iCell] = currentGamma_;
    muCells[iCell] = currentMu_;
    nuCells[iCell] = currentNu_;

    forAll(currentRhoe_, ie) {
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

    getDiffMassFluxes(tmpMassFluxP, tmpMassFluxT, tmpMassFluxXe, gradPCells[iCell], gradTCells[iCell], tmpGradXe, tmpDiffMassFluxes);

    forAll(tmpDiffMassFluxes, ie) {
      diffMassFluxesCells[ie][iCell] = tmpDiffMassFluxes[ie];
    }

    getDiffHeatFlux(energyFluxPCells[iCell], energyFluxTCells[iCell], tmpEnergyFluxXe, gradPCells[iCell], gradTCells[iCell], tmpGradXe, qDiffCells[iCell]);
  }

  updateBoundaries();
}

// ====================================================== //

void Foam::PTXEquilTabulated::updateCurrentPTY(scalar& p, scalar& T, scalarList& Y)
{
  currentP_ = p;
  currentT_ = T;
  currentYe_ = Y;
  currentYAr_ = currentYe_[indexAr_];
  currentPsi_ = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, *psiTable, currentP_, currentT_, currentYAr_);
  currentRho_ = currentPsi_*currentP_;
  currentRhoe_ = currentYe_*currentRho_;
  currentE_ = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, *eTable, currentP_, currentT_, currentYAr_);
  currentCv_ = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, *cvTable, currentP_, currentT_, currentYAr_);
  currentCp_ = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, *cpTable, currentP_, currentT_, currentYAr_);
  currentGamma_ = currentCp_/currentCv_;
  currentAlpha_ = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, *alphaTable, currentP_, currentT_, currentYAr_);
  currentSigma_ = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, *sigmaTable, currentP_, currentT_, currentYAr_);
  currentMu_ = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, *muTable, currentP_, currentT_, currentYAr_);
  currentNu_ = currentMu_/currentRho_;
  forAll(this->currentFluxesP_, i) {
    List<List<List<scalarList>>>& fp = *equilDiffFluxPTable;
    List<List<List<scalarList>>>& ft = *equilDiffFluxTTable;
    currentFluxesP_[i] = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, fp[i], currentP_, currentT_, currentYAr_);
    currentFluxesT_[i] = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, ft[i], currentP_, currentT_, currentYAr_);
  }
  forAll(this->currentFluxesP_, i) {
    List<List<List<scalarList>>>& fx = *equilDiffFluxXTable;
    currentFluxesX_[i] = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, fx[i], currentP_, currentT_, currentYAr_);
  }
}

// ====================================================== //

void Foam::PTXEquilTabulated::updateCurrentRhoEY(scalar& rho, scalar& Es, scalarList& Y)
{
  currentRho_ = rho;
  currentE_ = Es;
  currentYe_ = Y;
  currentYAr_ = currentYe_[indexAr_];
  getPTYfromRhoEY(rho, Es, currentP_, currentT_, currentYe_);
  currentPsi_ = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, *psiTable, currentP_, currentT_, currentYAr_);
  currentCv_ = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, *cvTable, currentP_, currentT_, currentYAr_);
  currentCp_ = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, *cpTable, currentP_, currentT_, currentYAr_);
  currentGamma_ = currentCp_/currentCv_;
  currentAlpha_ = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, *alphaTable, currentP_, currentT_, currentYAr_);
  currentSigma_ = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, *sigmaTable, currentP_, currentT_, currentYAr_);
  currentMu_ = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, *muTable, currentP_, currentT_, currentYAr_);
  currentNu_ = currentMu_/currentRho_;
  forAll(this->currentFluxesP_, i) {
    List<List<List<scalarList>>>& fp = *equilDiffFluxPTable;
    List<List<List<scalarList>>>& ft = *equilDiffFluxTTable;
    currentFluxesP_[i] = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, fp[i], currentP_, currentT_, currentYAr_);
    currentFluxesT_[i] = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, ft[i], currentP_, currentT_, currentYAr_);
  }
  forAll(this->currentFluxesP_, i) {
    List<List<List<scalarList>>>& fx = *equilDiffFluxXTable;
    currentFluxesX_[i] = trilinearInterpolation(*pressureTable, *temperatureTable, *argonMassFractionTable, fx[i], currentP_, currentT_, currentYAr_);
  }
  currentP_ = currentRho_/currentPsi_; // update p_{n+1}

}

// ====================================================== //

void Foam::PTXEquilTabulated::updateBoundaries()
{
  volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
  volScalarField::Boundary& pBf = this->p_.boundaryFieldRef();
  volScalarField::Boundary& eBf = this->e_.boundaryFieldRef();
  volScalarField::Boundary& rhoBf = this->rho_.boundaryFieldRef();
  volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
  volScalarField::Boundary& sigmaBf = this->sigma_.boundaryFieldRef();
  volScalarField::Boundary& alphaBf = this->alpha_.boundaryFieldRef();
  volScalarField::Boundary& gammaBf = this->gamma_.boundaryFieldRef();
  volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
  volScalarField::Boundary& nuBf = this->nu_.boundaryFieldRef();

  forAll(TBf, iPatch) {
    fvPatchScalarField& T_patch = TBf[iPatch];
    fvPatchScalarField& p_patch = pBf[iPatch];
    fvPatchScalarField& e_patch = eBf[iPatch];
    fvPatchScalarField& rho_patch = rhoBf[iPatch];
    fvPatchScalarField& psi_patch = psiBf[iPatch];
    fvPatchScalarField& sigma_patch = sigmaBf[iPatch];
    fvPatchScalarField& alpha_patch = alphaBf[iPatch];
    fvPatchScalarField& gamma_patch = gammaBf[iPatch];
    fvPatchScalarField& mu_patch = muBf[iPatch];
    fvPatchScalarField& nu_patch = nuBf[iPatch];

    forAll(T_patch, iFace) {
      if (T_patch.fixesValue()||isInConstructor_) {
        scalarList Y(nElements_);
        forAll(Y, eI) {
          Y[eI]=Ye_[eI].boundaryField()[iPatch][iFace];
        }
        updateCurrentPTY(p_patch[iFace], T_patch[iFace], Y);
        e_patch[iFace] = currentE_;
        rho_patch[iFace] = currentRho_;
        psi_patch[iFace] = currentPsi_;
        sigma_patch[iFace] = currentSigma_;
        alpha_patch[iFace] = currentAlpha_;
        gamma_patch[iFace] = currentGamma_;
        mu_patch[iFace] = currentMu_;
        nu_patch[iFace] = currentNu_;
      } else if (isA<zeroGradientFvPatchScalarField>(T_patch)||isA<fixedGradientFvPatchScalarField>(T_patch)) {
        scalarList Y(nElements_);
        forAll(Y, eI) {
          Y[eI]=Ye_[eI].boundaryField()[iPatch][iFace];
        }
        fixedGradientFvPatchScalarField& grade_patch=refCast<fixedGradientFvPatchScalarField>(e_patch);
        T_patch.evaluate();
        updateCurrentPTY(p_patch[iFace], T_patch[iFace], Y);
        rho_patch[iFace] = currentRho_;
        psi_patch[iFace] = currentPsi_;
        sigma_patch[iFace] = currentSigma_;
        alpha_patch[iFace] = currentAlpha_;
        gamma_patch[iFace] = currentGamma_;
        mu_patch[iFace] = currentMu_;
        nu_patch[iFace] = currentNu_;
        scalar cv = currentCv_;
        scalar e = currentE_;
        label faceCellIndex = mesh_.boundaryMesh()[iPatch].faceCells()[iFace];
        updateCurrentPTY(p_[faceCellIndex], T_[faceCellIndex], Y);
        scalar efc = currentE_;
        grade_patch.gradient()[iFace]= cv * fvc::snGrad(T_)->boundaryField()[iPatch][iFace]
                                       + mesh_.deltaCoeffs().boundaryField()[iPatch][iFace] * (e-efc);
        grade_patch.updateCoeffs();
      } else {
        FatalErrorInFunction << "This BC is not implemented yet." << exit(FatalError);
      }
    }
  }
}

// ====================================================== //

Foam::scalar Foam::PTXEquilTabulated::trilinearInterpolation(scalarList x, scalarList y, scalarList z, List<List<scalarList> > f_xyz, scalar x0, scalar y0, scalar z0)
{
  if(x.size()<2) {
    FatalErrorInFunction << "trilinearInterpolation problem: x.size() must have a size greater than 1." << exit(FatalError);
  }
  if(y.size()<2) {
    FatalErrorInFunction << "trilinearInterpolation problem: y.size() must have a size greater than 1." << exit(FatalError);
  }
  if(z.size()<2) {
    FatalErrorInFunction << "trilinearInterpolation problem: z.size() must have a size greater than 1." << exit(FatalError);
  }
  // order the x, y, z and f_xyz tables
  labelList visitOrder_x;
  sortedOrder(x, visitOrder_x);
  x = scalarList(x, visitOrder_x);
  f_xyz = List<List<scalarList>>(f_xyz, visitOrder_x);
  labelList visitOrder_y;
  sortedOrder(y, visitOrder_y);
  y = scalarList(y, visitOrder_y);
  labelList visitOrder_z;
  sortedOrder(z, visitOrder_z);
  z = scalarList(z, visitOrder_z);
  forAll(f_xyz, yI) {
    f_xyz[yI] = List<scalarList>(f_xyz[yI], visitOrder_y);
  }
  forAll(f_xyz, yI) {
    forAll(f_xyz[yI], zI) {
      f_xyz[yI][zI] = scalarList(f_xyz[yI][zI], visitOrder_z);
    }
  }
  // lower and upper values
  if      (x0<=x[0]) x0 = x[0]+1e-10;
  else if (x0>=x[x.size()-1])  x0 = x[x.size()-1]-1e-10;

  if      (y0<=y[0]) y0 = y[0]+1e-10;
  else if (y0>=y[y.size()-1])  y0 = y[y.size()-1]-1e-10;

  if      (z0<=z[0]) z0 = z[0]+1e-10;
  else if (z0>=z[z.size()-1])  z0 = z[z.size()-1]-1e-10;

  // x,y,z indexes
  int idz = -1;
  const int zsize = z.size();
  forAll(z, i) {
    if (i < (zsize-1)) {
      if (z0 >= z[i] && z0 < z[i+1]) {
        idz = i;
      }
    }
  }
  assert(idz>=0);
  assert(idz<zsize-1);
  assert(z[idz]<=z0);
  assert(z[idz+1]>z0);

  int idy = -1;
  const int ysize = y.size();
  forAll(y, i) {
    if (i < (ysize-1)) {
      if (y0 >= y[i] && y0 < y[i+1]) {
        idy = i;
      }
    }
  }

  assert(idy>=0);
  assert(idy<ysize-1);
  assert(y[idy]<=y0);
  assert(y[idy+1]>y0);

  int idx = -1;
  const int xsize = x.size();
  forAll(x, i) {
    if (i < (xsize-1)) {
      if (x0 >= x[i] && x0 < x[i+1]) {
        idx = i;
      }
    }
  }
  assert(idx>=0);
  assert(idx<xsize-1);
  assert(x[idx]<=x0);
  assert(x[idx+1]>x0);

  // linear interpolation
  scalar percent_x = (x0 - x[idx]) / (x[idx+1] - x[idx]);
  scalar percent_y = (y0 - y[idy]) / (y[idy+1] - y[idy]);
  scalar percent_z = (z0 - z[idz]) / (z[idz+1] - z[idz]);

  scalar yLzL =  f_xyz[idx][idy][idz] + percent_x * ( f_xyz[idx+1][idy][idz] - f_xyz[idx][idy][idz] ) ;
  scalar yHzL =  f_xyz[idx][idy+1][idz] + percent_x * ( f_xyz[idx+1][idy+1][idz] - f_xyz[idx][idy+1][idz] ) ;

  scalar yLzH =  f_xyz[idx][idy][idz+1] + percent_x * ( f_xyz[idx+1][idy][idz+1] - f_xyz[idx][idy][idz+1] ) ;
  scalar yHzH =  f_xyz[idx][idy+1][idz+1] + percent_x * ( f_xyz[idx+1][idy+1][idz+1] - f_xyz[idx][idy+1][idz+1] ) ;

  scalar zL = yLzL + percent_y * (yHzL - yLzL);
  scalar zH = yLzH + percent_y * (yHzH - yLzH);

  scalar f_x0y0z0 =  zL + percent_z * (zH - zL) ;
  return f_x0y0z0;


}

// ====================================================== //

void Foam::PTXEquilTabulated::generateTable(const word tableName)
{
  if (isFile(basicTabulated::tableFileName_)) {
    FatalErrorInFunction << "Problem to generate the thermo table: the file " << basicTabulated::tableFileName_ << " exists already." << exit(FatalError);
  }

  dictionary generateTableDict_(basicTabulated::thermoDict_.subDict("generateThermoTable"));
  scalarList pressureList(generateTableDict_.lookup("pressureList"));
  scalarList temperatureList(generateTableDict_.lookup("temperatureList"));
  List<scalarList> elementsMassFractionList(generateTableDict_.lookup("elementsMassFractionList"));
  List<List<scalarList>> eTable_(pressureList.size());
  List<List<scalarList>> psiTable_(pressureList.size());
  List<List<scalarList>> alphaTable_(pressureList.size());
  List<List<scalarList>> cpTable_(pressureList.size());
  List<List<scalarList>> cvTable_(pressureList.size());
  List<List<scalarList>> sigmaTable_(pressureList.size());
  List<List<scalarList>> muTable_(pressureList.size());
  List<List<List<scalarList>>> equilDiffFluxPTable_(nElements_+1);
  List<List<List<scalarList>>> equilDiffFluxTTable_(nElements_+1);
  List<List<List<scalarList>>> equilDiffFluxXTable_(nElements_*(nElements_+1));
  forAll(equilDiffFluxPTable_, i) {
    equilDiffFluxPTable_[i].resize(pressureList.size());
    equilDiffFluxTTable_[i].resize(pressureList.size());
  }
  forAll(equilDiffFluxXTable_, i) {
    equilDiffFluxXTable_[i].resize(pressureList.size());
  }

  scalar* pVars = new scalar[2];
  scalar* pXe   = new scalar[this->nElements()];
  scalar* pYe   = new scalar[this->nElements()];

  forAll(pressureList, pI) {
    eTable_[pI].resize(temperatureList.size());
    psiTable_[pI].resize(temperatureList.size());
    alphaTable_[pI].resize(temperatureList.size());
    cpTable_[pI].resize(temperatureList.size());
    cvTable_[pI].resize(temperatureList.size());
    sigmaTable_[pI].resize(temperatureList.size());
    muTable_[pI].resize(temperatureList.size());
    forAll(equilDiffFluxPTable_, i) {
      equilDiffFluxPTable_[i][pI].resize(temperatureList.size());
      equilDiffFluxTTable_[i][pI].resize(temperatureList.size());
    }
    forAll(equilDiffFluxXTable_, i) {
      equilDiffFluxXTable_[i][pI].resize(temperatureList.size());
    }

    forAll(temperatureList, tI) {
      eTable_[pI][tI].resize(elementsMassFractionList.size());
      psiTable_[pI][tI].resize(elementsMassFractionList.size());
      alphaTable_[pI][tI].resize(elementsMassFractionList.size());
      cpTable_[pI][tI].resize(elementsMassFractionList.size());
      cvTable_[pI][tI].resize(elementsMassFractionList.size());
      sigmaTable_[pI][tI].resize(elementsMassFractionList.size());
      muTable_[pI][tI].resize(elementsMassFractionList.size());
      forAll(equilDiffFluxPTable_, i) {
        equilDiffFluxPTable_[i][pI][tI].resize(elementsMassFractionList.size());
        equilDiffFluxTTable_[i][pI][tI].resize(elementsMassFractionList.size());
      }
      forAll(equilDiffFluxXTable_, i) {
        equilDiffFluxXTable_[i][pI][tI].resize(elementsMassFractionList.size());
      }

      forAll(elementsMassFractionList, eI) {

        pVars[0] = pressureList[pI];
        pVars[1] = temperatureList[tI];
        forAll(elementsMassFractionList[eI], ie) {
          pYe[ie] = elementsMassFractionList[eI][ie];
        }
        this->pMutationMix_->convert<Mutation::Thermodynamics::ConversionType::YE_TO_XE>(pYe, pXe);
        this->pMutationMix_->setState(pXe, pVars, 2);
        eTable_[pI][tI][eI] = this->pMutationMix_->mixtureEnergyMass();
        psiTable_[pI][tI][eI] = this->pMutationMix_->density()/pressureList[pI];
        alphaTable_[pI][tI][eI] = this->pMutationMix_->equilibriumThermalConductivity() / this->pMutationMix_->mixtureEquilibriumCvMass();
        cpTable_[pI][tI][eI] = this->pMutationMix_->mixtureEquilibriumCpMass();
        cvTable_[pI][tI][eI] = this->pMutationMix_->mixtureEquilibriumCvMass();
        sigmaTable_[pI][tI][eI] = pMutationMix_->electricConductivity();
        muTable_[pI][tI][eI] =  pMutationMix_->viscosity();
        scalarList fp(nElements_+1);
        scalarList ft(nElements_+1);
        scalarList fx(nElements_*(nElements_+1));
        getElementFluxes(fp,ft,fx);
        forAll(equilDiffFluxPTable_, i) {
          equilDiffFluxPTable_[i][pI][tI][eI] = fp[i];
          equilDiffFluxTTable_[i][pI][tI][eI] = ft[i];
        }
        forAll(equilDiffFluxXTable_, i) {
          equilDiffFluxXTable_[i][pI][tI][eI] = fx[i];
        }
      }
    }
  }

  pressureTable = new scalarList(pressureList);
  temperatureTable=new scalarList(temperatureList);
  elementsMassFractionTable=new List<scalarList>(elementsMassFractionList);
  argonMassFractionTable=new scalarList();
  List<scalarList>& Y = *elementsMassFractionTable;
  forAll(Y, tableI) {
    argonMassFractionTable->append(Y[tableI][indexAr_]);
  }
  eTable=new List<List<scalarList>>(eTable_);
  psiTable=new List<List<scalarList>>(psiTable_);
  alphaTable=new List<List<scalarList>>(alphaTable_);
  cpTable=new List<List<scalarList>>(cpTable_);
  cvTable=new List<List<scalarList>>(cvTable_);
  sigmaTable=new List<List<scalarList>>(sigmaTable_);
  muTable=new List<List<scalarList>>(muTable_);
  equilDiffFluxPTable=new List<List<List<scalarList>>>(equilDiffFluxPTable_);
  equilDiffFluxTTable=new List<List<List<scalarList>>>(equilDiffFluxTTable_);
  equilDiffFluxXTable=new List<List<List<scalarList>>>(equilDiffFluxXTable_);

  // Writing the thermo table
  IOdictionary thermoTableDict_
  (
      IOobject
      (
          tableName,
          basicUserThermo::mesh_.time().constant(),
          basicUserThermo::mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE,
          false
      )
  );
  thermoTableDict_.set("p", pressureList);
  thermoTableDict_.set("T", temperatureList);
  thermoTableDict_.set("Elements", getElementsList());
  thermoTableDict_.set("Y", elementsMassFractionList);
  thermoTableDict_.set("e", eTable_);
  thermoTableDict_.set("psi", psiTable_);
  thermoTableDict_.set("alpha", alphaTable_);
  thermoTableDict_.set("cp", cpTable_);
  thermoTableDict_.set("cv", cvTable_);
  thermoTableDict_.set("sigma", sigmaTable_);
  thermoTableDict_.set("mu", muTable_);
  thermoTableDict_.set("equilDiffFluxP", equilDiffFluxPTable_);
  thermoTableDict_.set("equilDiffFluxT", equilDiffFluxTTable_);
  thermoTableDict_.set("equilDiffFluxX", equilDiffFluxXTable_);
  thermoTableDict_.regIOobject::write();

  Info << "Writing the thermo table in " << basicTabulated::tableFileName_ << endl;
}

// ====================================================== //

void Foam::PTXEquilTabulated::readTable(const word tableName)
{
  Info << "Reading the thermo table in " << basicTabulated::tableFileName_ << endl;
  IOdictionary thermoTableDict_
  (
      IOobject
      (
          tableName,
          basicUserThermo::mesh_.time().constant(),
          basicUserThermo::mesh_,
          IOobject::MUST_READ,
          IOobject::NO_WRITE,
          false
      )
  );
  wordList elementsTable(thermoTableDict_.lookup("Elements"));
  wordList elementsMutation(getElementsList());
  if (elementsMutation.size() != elementsTable.size()) {
    FatalErrorInFunction << "Problem in " << tableFileName_ << endl;
    FatalErrorInFunction << "Number of elements in mixture (" <<elementsMutation.size()  << ") different than the number of elements in table (" << elementsTable.size() << ")" << exit(FatalError);
  }
  forAll(elementsMutation, eI) {
    if (elementsMutation[eI]!=elementsTable[eI]) {
      FatalErrorInFunction << "Problem in " << tableFileName_ << endl;
      FatalErrorInFunction << "The elements are different between the mixture:" << elementsMutation << endl;
      FatalErrorInFunction << "and the table:" << elementsTable << exit(FatalError);
    }
  }
  pressureTable = new scalarList(thermoTableDict_.lookup("p"));
  temperatureTable=new scalarList(thermoTableDict_.lookup("T"));
  elementsMassFractionTable=new List<scalarList>(thermoTableDict_.lookup("Y"));
  argonMassFractionTable=new scalarList();
  List<scalarList>& Y = *elementsMassFractionTable;
  forAll(Y, tableI) {
    argonMassFractionTable->append(Y[tableI][indexAr_]);
  }
  eTable=new List<List<scalarList>>(thermoTableDict_.lookup("e"));
  psiTable=new List<List<scalarList>>(thermoTableDict_.lookup("psi"));
  alphaTable=new List<List<scalarList>>(thermoTableDict_.lookup("alpha"));
  cpTable=new List<List<scalarList>>(thermoTableDict_.lookup("cp"));
  cvTable=new List<List<scalarList>>(thermoTableDict_.lookup("cv"));
  sigmaTable=new List<List<scalarList>>(thermoTableDict_.lookup("sigma"));
  muTable=new List<List<scalarList>>(thermoTableDict_.lookup("mu"));
  equilDiffFluxPTable=new List<List<List<scalarList>>>(thermoTableDict_.lookup("equilDiffFluxP"));
  equilDiffFluxTTable=new List<List<List<scalarList>>>(thermoTableDict_.lookup("equilDiffFluxT"));
  equilDiffFluxXTable=new List<List<List<scalarList>>>(thermoTableDict_.lookup("equilDiffFluxX"));
}

// ====================================================== //

void Foam::PTXEquilTabulated::getPTYfromRhoEY(scalar rho, scalar Es, scalar p, scalar T, scalarList Y, scalar tol, scalar maxIter)
{
  scalar Test = T;
  scalar Tnew = T;
  scalar Ttol = T*tol;
  int    iter = 0;
  do {
    Test = Tnew;
    updateCurrentPTY(p, Test, Y);
    Tnew = Test - (currentE_ - Es)/currentCv_;

    if (iter++ > maxIter) {
      FatalErrorInFunction
          << "Maximum number of iterations exceeded"
          << abort(FatalError);
    }

  } while (mag(Tnew - Test) > Ttol);
//  int    iter = 0;
//  int    iterT = 0;
//  int    iterP = 0;

//	updateCurrentPTY(p, T, Y);

//	scalar Rgas = currentCp_ - currentCv_;
//	scalar rhoE = rho * Es;

//	scalar deltaRhoE = rho * (Es - currentE_);
//	scalar deltaRho  = rho - currentRho_;
//	scalar dRhoEdT = rho * currentCv_; 										//- ATM: this is not accurate, but probably works
//	scalar dRhodT  = - rho * currentCv_ / Es;							//- ATM: this is not accurate, but probably works
//	scalar dRhoEdP = 1 / (currentGamma_ - 1);							//- ATM: this is not accurate, but probably works
//	scalar dRhodP  = 1 / (Rgas * (T - Es / currentCv_));	//- ATM: this is not accurate, but probably works
//	scalar deltaT = 0.;
//	scalar deltaP = 0.;

//	while (std::max(std::abs(deltaRhoE / rhoE), std::abs(deltaRho / rho)) > tol)
//	{
//		deltaT = deltaRhoE / dRhoEdT + deltaRho / dRhodT;
//		deltaP = deltaRhoE / dRhoEdP + deltaRho / dRhodP;

//        iterT = 0;
//        iterP = 0;
//        while (T + deltaT < 0)	{ deltaT *= 0.5;
//            if (iterT++ > maxIter)
//                {
//              FatalErrorInFunction
//                << "Maximum number of iterations exceeded in TE function (computing p and T for given rho and e)."
//                << abort(FatalError);
//                }}
//        while (p + deltaP < 0)	{ deltaP *= 0.5;
//            if (iterP++ > maxIter)
//                {
//              FatalErrorInFunction
//                << "Maximum number of iterations exceeded in TE function (computing p and T for given rho and e)."
//                << abort(FatalError);
//                }}

//		T += deltaT;
//		p += deltaP;

//		updateCurrentPTY(p, T, Y);

//		Rgas = currentCp_ - currentCv_;

//		deltaRhoE = rho * (Es - currentE_);
//		deltaRho  = rho - currentRho_;     //rho - (p / (Rgas * T));
//		dRhoEdT = rho * currentCv_;
//		dRhodT  = - rho * currentCv_ / Es;
//		dRhoEdP = 1 / (currentGamma_ - 1);
//		dRhodP  = 1 / (Rgas * T - Rgas * Es / currentCv_);

//    if (iter++ > maxIter)
//		{
//      FatalErrorInFunction
//        << "Maximum number of iterations exceeded in TE function (computing p and T for given rho and e)."
//        << abort(FatalError);
//		}
//	}

//	this->currentP_ = p;
//	this->currentT_ = T;
}

// ====================================================== //

void Foam::PTXEquilTabulated::getElementFluxes(scalarList& fluxP, scalarList& fluxT, scalarList& fluxX)
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

void Foam::PTXEquilTabulated::getDiffMassFluxes(const scalarList& fluxP, const scalarList& fluxT, const scalarList& fluxX, const vector& gradP, const vector& gradT, const List<vector>& gradXe, List<vector>& diffMF)
{
  forAll(fluxP, ie) {
    diffMF[ie] = fluxP[ie] * gradP + fluxT[ie] * gradT;
    forAll(fluxP, ie2) {
      diffMF[ie] += fluxX[ie*this->nElements() + ie2] * gradXe[ie2]; 			//- ATM: check the order of the terms TO-DO
    }
  }
}

// ====================================================== //

void Foam::PTXEquilTabulated::getDiffHeatFlux(const scalar& eFluxP, const scalar& eFluxT, const scalarList& eFluxXe, vector& gradP, vector& gradT, List<vector>& gradXe, vector& qDiff)
{
  qDiff = eFluxP * gradP + eFluxT * gradT;
  forAll(eFluxXe, ie) {
    qDiff += eFluxXe[ie] * gradXe[ie];
  }
}
// ************************************************************************* //

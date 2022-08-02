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

#include "PTEquilTabulated.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PTEquilTabulated::PTEquilTabulated
(
    const fvMesh& mesh,
    const word& dictName
)
  :
basicTabulated(mesh, dictName)
{
  stateModelName_ = "Equil";

  configureMixture();
  createCompositionFields(mesh);
  if(!isFile(basicTabulated::tableFileName_)) {
    generateTable(basicTabulated::tableName_);
  } else {
    readTable(basicTabulated::tableName_);
  }
  update();
  isInConstructor_ = 0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PTEquilTabulated::~PTEquilTabulated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PTEquilTabulated::update()
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

  forAll(TCells, iCell) {
    if (isInConstructor_) {
      updateCurrentPT(pCells[iCell], TCells[iCell]);
      eCells[iCell] = currentE_;
      rhoCells[iCell] = currentRho_;
      psiCells[iCell] = currentPsi_;
      sigmaCells[iCell] = currentSigma_;
      alphaCells[iCell] = currentAlpha_;
      gammaCells[iCell] = currentGamma_;
      muCells[iCell] = currentMu_;
      nuCells[iCell] = currentNu_;
    } else {
      updateCurrentRhoE(rhoCells[iCell], eCells[iCell]);
      TCells[iCell] = currentT_;
      pCells[iCell] = currentP_;
      psiCells[iCell] = currentPsi_;
      sigmaCells[iCell] = currentSigma_;
      alphaCells[iCell] = currentAlpha_;
      gammaCells[iCell] = currentGamma_;
      muCells[iCell] = currentMu_;
      nuCells[iCell] = currentNu_;
    }
  }
  updateBoundaries();
}

// ====================================================== //

void Foam::PTEquilTabulated::updateCurrentPT(scalar& p, scalar& T)
{
  currentP_ = p;
  currentT_ = T;
  currentPsi_ = bilinearInterpolation(*pressureTable, *temperatureTable, *psiTable, currentP_, currentT_);
  currentRho_ = currentPsi_*currentP_;
  currentE_ = bilinearInterpolation(*pressureTable, *temperatureTable, *eTable, currentP_, currentT_);
  currentCv_ = bilinearInterpolation(*pressureTable, *temperatureTable, *cvTable, currentP_, currentT_);
  currentCp_ = bilinearInterpolation(*pressureTable, *temperatureTable, *cpTable, currentP_, currentT_);
  currentGamma_ = currentCp_/currentCv_;
  currentAlpha_ = bilinearInterpolation(*pressureTable, *temperatureTable, *alphaTable, currentP_, currentT_);
  currentSigma_ = bilinearInterpolation(*pressureTable, *temperatureTable, *sigmaTable, currentP_, currentT_);
  currentMu_ = bilinearInterpolation(*pressureTable, *temperatureTable, *muTable, currentP_, currentT_);
  currentNu_ = currentMu_/currentRho_;
}

// ====================================================== //

void Foam::PTEquilTabulated::updateCurrentRhoE(scalar& rho, scalar& Es)
{
  currentRho_ = rho;
  currentE_ = Es;
  getPTfromRhoE(rho, Es, currentP_, currentT_);
  currentPsi_ = bilinearInterpolation(*pressureTable, *temperatureTable, *psiTable, currentP_, currentT_);
  currentCv_ = bilinearInterpolation(*pressureTable, *temperatureTable, *cvTable, currentP_, currentT_);
  currentCp_ = bilinearInterpolation(*pressureTable, *temperatureTable, *cpTable, currentP_, currentT_);
  currentGamma_ = currentCp_/currentCv_;
  currentAlpha_ = bilinearInterpolation(*pressureTable, *temperatureTable, *alphaTable, currentP_, currentT_);
  currentSigma_ = bilinearInterpolation(*pressureTable, *temperatureTable, *sigmaTable, currentP_, currentT_);
  currentMu_ = bilinearInterpolation(*pressureTable, *temperatureTable, *muTable, currentP_, currentT_);
  currentNu_ = currentMu_/currentRho_;
}


// ====================================================== //

void Foam::PTEquilTabulated::updateBoundaries()
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
      if (T_patch.fixesValue()) {
        updateCurrentPT(p_patch[iFace], T_patch[iFace]);
        e_patch[iFace] = currentE_;
        rho_patch[iFace] = currentRho_;
        psi_patch[iFace] = currentPsi_;
        sigma_patch[iFace] = currentSigma_;
        alpha_patch[iFace] = currentAlpha_;
        gamma_patch[iFace] = currentGamma_;
        mu_patch[iFace] = currentMu_;
        nu_patch[iFace] = currentNu_;
      } else if (isA<zeroGradientFvPatchScalarField>(T_patch)||isA<fixedGradientFvPatchScalarField>(T_patch)) {
        fixedGradientFvPatchScalarField& grade_patch=refCast<fixedGradientFvPatchScalarField>(e_patch);
        T_patch.evaluate();
        updateCurrentPT(p_patch[iFace], T_patch[iFace]);
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
        updateCurrentPT(p_[faceCellIndex], T_[faceCellIndex]);
        scalar efc = currentE_;
        grade_patch.gradient()[iFace]= cv * fvc::snGrad(T_)->boundaryField()[iPatch][iFace]
                                       + mesh_.deltaCoeffs().boundaryField()[iPatch][iFace] * (e-efc);
        grade_patch.updateCoeffs();
      } else {
        FatalErrorInFunction << "This BC is not implemented yet." << exit(FatalError);
//                e_.correctBoundaryConditions();
//                rho_.correctBoundaryConditions();
//                updateCurrentRhoE(rho_patch[iFace], e_patch[iFace]); // BC energy to verify
//                T_patch[iFace] = currentT_;
//                p_patch[iFace] = currentP_;
//                psi_patch[iFace] = currentPsi_;
//                sigma_patch[iFace] = currentSigma_;
//                alpha_patch[iFace] = currentAlpha_;
//                gamma_patch[iFace] = currentGamma_;
//                mu_patch[iFace] = currentMu_;
//                nu_patch[iFace] = currentNu_;
      }
    }
  }
}

// ====================================================== //

Foam::scalar Foam::PTEquilTabulated::bilinearInterpolation(scalarList x, scalarList y, List<scalarList> f_xy, scalar x0, scalar y0)
{
  if(x.size()<2) {
    FatalErrorInFunction << "bilinearInterpolation problem: x.size() must have a size greater than 1." << exit(FatalError);
  }
  if(y.size()<2) {
    FatalErrorInFunction << "bilinearInterpolation problem: y.size() must have a size greater than 1." << exit(FatalError);
  }

  // order the x, y and f_xy tables
  labelList visitOrder;
  sortedOrder(x, visitOrder);
  x = scalarList(x, visitOrder);
  f_xy = List<scalarList>(f_xy, visitOrder);
  sortedOrder(y, visitOrder);
  y = scalarList(y, visitOrder);
  forAll(f_xy, yI) {
    f_xy[yI] = scalarList(f_xy[yI], visitOrder);
  }

  // lower and upper values
  if      (x0<=x[0]) x0 = x[0]+1e-10;
  else if (x0>=x[x.size()-1])  x0 = x[x.size()-1]-1e-10;

  if      (y0<=y[0]) y0 = y[0]+1e-10;
  else if (y0>=y[y.size()-1])  y0 = y[y.size()-1]-1e-10;

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

  scalar percent_x = (x0 - x[idx]) / (x[idx+1] - x[idx]);
  scalar percent_y = (y0 - y[idy]) / (y[idy+1] - y[idy]);
  scalar yL =  f_xy[idx][idy] + percent_x * ( f_xy[idx+1][idy] - f_xy[idx][idy] ) ;
  scalar yH =  f_xy[idx][idy+1] + percent_x * ( f_xy[idx+1][idy+1] - f_xy[idx][idy+1] ) ;
  scalar f_x0y0 =  yL + percent_y * (yH - yL) ;
  return f_x0y0;


}

// ====================================================== //

void Foam::PTEquilTabulated::generateTable(const word tableName)
{
  if (isFile(basicTabulated::tableFileName_)) {
    FatalErrorInFunction << "Problem to generate the thermo table: the file " << basicTabulated::tableFileName_ << " exists already." << exit(FatalError);
  }

  dictionary generateTableDict_(basicTabulated::thermoDict_.subDict("generateThermoTable"));
  scalarList pressureList(generateTableDict_.lookup("pressureList"));
  scalarList temperatureList(generateTableDict_.lookup("temperatureList"));
  List<scalarList> eTable_(pressureList.size());
  List<scalarList> psiTable_(pressureList.size());
  List<scalarList> alphaTable_(pressureList.size());
  List<scalarList> cpTable_(pressureList.size());
  List<scalarList> cvTable_(pressureList.size());
  List<scalarList> sigmaTable_(pressureList.size());
  List<scalarList> muTable_(pressureList.size());
  forAll(pressureList, pI) {
    eTable_[pI].resize(temperatureList.size());
    psiTable_[pI].resize(temperatureList.size());
    alphaTable_[pI].resize(temperatureList.size());
    cpTable_[pI].resize(temperatureList.size());
    cvTable_[pI].resize(temperatureList.size());
    sigmaTable_[pI].resize(temperatureList.size());
    muTable_[pI].resize(temperatureList.size());
    forAll(temperatureList, tI) {
      this->pMutationMix_->setState(&pressureList[pI], &temperatureList[tI], 1);
      eTable_[pI][tI] = this->pMutationMix_->mixtureEnergyMass();
      psiTable_[pI][tI] = this->pMutationMix_->density()/pressureList[pI];
      alphaTable_[pI][tI] = this->pMutationMix_->equilibriumThermalConductivity() / this->pMutationMix_->mixtureEquilibriumCvMass();
      cpTable_[pI][tI] = this->pMutationMix_->mixtureEquilibriumCpMass();
      cvTable_[pI][tI] = this->pMutationMix_->mixtureEquilibriumCvMass();
      sigmaTable_[pI][tI] = pMutationMix_->electricConductivity();
      muTable_[pI][tI] =  pMutationMix_->viscosity();
    }
  }

  pressureTable = new scalarList(pressureList);
  temperatureTable=new scalarList(temperatureList);
  eTable=new List<scalarList>(eTable_);
  psiTable=new List<scalarList>(psiTable_);
  alphaTable=new List<scalarList>(alphaTable_);
  cpTable=new List<scalarList>(cpTable_);
  cvTable=new List<scalarList>(cvTable_);
  sigmaTable=new List<scalarList>(sigmaTable_);
  muTable=new List<scalarList>(muTable_);

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
  thermoTableDict_.set("e", eTable_);
  thermoTableDict_.set("psi", psiTable_);
  thermoTableDict_.set("alpha", alphaTable_);
  thermoTableDict_.set("cp", cpTable_);
  thermoTableDict_.set("cv", cvTable_);
  thermoTableDict_.set("sigma", sigmaTable_);
  thermoTableDict_.set("mu", muTable_);
  thermoTableDict_.regIOobject::write();

  Info << "Writing the thermo table in " << basicTabulated::tableFileName_ << endl;
}

// ====================================================== //

void Foam::PTEquilTabulated::readTable(const word tableName)
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
  pressureTable = new scalarList(thermoTableDict_.lookup("p"));
  temperatureTable=new scalarList(thermoTableDict_.lookup("T"));
  eTable=new List<scalarList>(thermoTableDict_.lookup("e"));
  psiTable=new List<scalarList>(thermoTableDict_.lookup("psi"));
  alphaTable=new List<scalarList>(thermoTableDict_.lookup("alpha"));
  cpTable=new List<scalarList>(thermoTableDict_.lookup("cp"));
  cvTable=new List<scalarList>(thermoTableDict_.lookup("cv"));
  sigmaTable=new List<scalarList>(thermoTableDict_.lookup("sigma"));
  muTable=new List<scalarList>(thermoTableDict_.lookup("mu"));
}

// ====================================================== //

void Foam::PTEquilTabulated::getPTfromRhoE(scalar& rho, scalar& Es, scalar& p, scalar& T, scalar tol, scalar maxIter)
{
  scalar Test = T;
  scalar Tnew = T;
  scalar Ttol = T*tol;
  int    iter = 0;
  do {
    Test = Tnew;
    updateCurrentPT(p, Test);
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

//	updateCurrentPT(p, T);

//	scalar Rgas = currentCp_ - currentCv_;
//	scalar rhoE = rho * Es;

//	scalar deltaRhoE = rho * (Es - currentE_);
//	scalar deltaRho  = rho - currentRho_;
//	scalar dRhoEdT = rho * currentCv_;                                  //- ATM: this is not accurate, but probably works
//	scalar dRhodT  = - rho * currentCv_ / Es;							//- ATM: this is not accurate, but probably works
//	scalar dRhoEdP = 1 / (currentGamma_ - 1);							//- ATM: this is not accurate, but probably works
//	scalar dRhodP  = 1 / (Rgas * (T - Es / currentCv_));                //- ATM: this is not accurate, but probably works
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

//		updateCurrentPT(p, T);

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
}

// ************************************************************************* //

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
    along with OpenFOAM.  If LinearArrheniust, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LinearArrheniusPyrolysisModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LinearArrheniusPyrolysisModel::LinearArrheniusPyrolysisModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simplePyrolysisModel(mesh, regionName),
nSolidPhases_(createScalarProp("nSolidPhases")),
tau(createVolField<scalar>("tau",dimensionedScalar("1", dimless, 1.0),wordList(nPatches_,"zeroGradient"))),
piTotal_(createVolField<scalar>("piTotal",dimensionedScalar("0",dimMass/dimVolume/dimTime,0))),
rho_v_(createVolField<scalar>("rho_v",dimensionedScalar("0",dimMass/dimVolume,0))),
rho_c_(createVolField<scalar>("rho_c",dimensionedScalar("0",dimMass/dimVolume,0))),
normFp_sum(createVolField<scalar>("normFp_sum", dimensionedScalar("0", dimMass/dimVolume, 0.0))),
materialChemistryModel(refModel<simpleMaterialChemistryModel>()),
speciesNames_(materialChemistryModel.speciesNames()),
elementNames_(materialChemistryModel.elementNames()),
energyModel(refModel<simpleEnergyModel>()),
Ta_(energyModel.refVolField<scalar>("Ta")),
rho_s_(energyModel.refVolField<scalar>("rho_s")),
R_(::constant::physicoChemical::R),
sizeXsi(0),
hp(simplePyrolysisModel::hp_),
MaterialChemistryType_(simplePyrolysisModel::materialDict_.isDict("MaterialChemistry")?word(simplePyrolysisModel::materialDict_.subDict("MaterialChemistry").lookup("MaterialChemistryType")):""),
oneVolScalarField(createVolField<scalar>("one",dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0))),
minPiPyro_(simplePyrolysisModel::materialDict_.subDict("Pyrolysis").template lookupOrDefault("minPiPyro",1e-9))
{

  if (this->debug_) {
    Info << getTabLevel() << "debug: start --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }

  if (MaterialChemistryType_!="") {
    if ( !(MaterialChemistryType_ == "ConstantEquilibrium" || MaterialChemistryType_ == "SpeciesConservation"
           || MaterialChemistryType_ == "EquilibriumElement" || MaterialChemistryType_ == "ConstantFiniteRate" )) {

      FatalErrorInFunction << "Unknown MaterialChemistryTransport type \n\n"
                           <<  "Valid MaterialChemistryTransport types using Pyrolysis are:\n\n"
                           <<  "   ConstantEquilibrium \n"
                           <<  "       (MaterialChemistry : equilibrium, transport : averaged-momentum) \n"
                           <<  "   EquilibriumElement \n"
                           <<  "       (MaterialChemistry : equilibrium, transport : averaged-momentum + element conservation) \n"
                           <<  "   ConstantFiniteRate \n"
                           <<  "      (MaterialChemistry : finite rate, transport : averaged-momentum) \n"
                           <<  "   SpeciesConservation \n"
                           <<  "      (MaterialChemistry : finite rate, transport : averaged-momentum + species conservation) "
                           << exit(FatalError);
    }
  }


  if (nSolidPhases_ < 1) {
    FatalErrorInFunction
        << "The number of solid phases 'nSolidPhases' must be > 0. "
        "If not, just use a CFD solver!" << nl
        << exit(FatalError);
  }

  if (this->debug_) {
    Info << getTabLevel() << "debug: nPyroReac --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  nPyroReac.resize(nSolidPhases_);
  forAll(nPyroReac, solidI) {
    const word namePyro = "nPyroReac[" + Foam::name(solidI + 1) + "]";
    nPyroReac.set(solidI, createScalarProp(namePyro,"yes")); // default value = 0
    sizeXsi += nPyroReac[solidI];
  }
  hashXsi.setSize(sizeXsi, 2);
  int nindex = 0;

  if (this->debug_) {
    Info << getTabLevel() << "debug: hashXsi --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  for (int solidI = 0; solidI < nSolidPhases_; solidI++) {
    for (int reacI = 0; reacI < nPyroReac[solidI]; reacI++) {
      hashXsi[nindex + reacI][0] = solidI + 1;
      hashXsi[nindex + reacI][1] = reacI + 1;
    }
    nindex += nPyroReac[solidI];
  }

  /* Pyrolysis advancements */

  wordList zeroGradBC(nPatches_,"zeroGradient");

  if (this->debug_) {
    Info << getTabLevel() << "debug: Xsi --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }

  Xsi.resize(sizeXsi);
  forAll(Xsi, reacI) {
    word nameXsi = "Xsi[" + std::to_string(hashXsi[reacI][0]) + "][" +
                   std::to_string(hashXsi[reacI][1]) + "]";
    Xsi.set(reacI,createVolField<scalar>(nameXsi,dimensionedScalar("zero", dimless, 0.0),zeroGradBC));
  }
  /* END:  Pyrolysis advancements  */

  if (this->debug_) {
    Info << getTabLevel() << "debug: F,A,E,m,n,T,h --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  /* Pyrolysis reactions parameters */

  forAll(Xsi, i) {
    word phase_reac="[" + std::to_string(hashXsi[i][0]) + "][" + std::to_string(hashXsi[i][1]) + "]";
    Ap.append(createDimScalarProp("A"+phase_reac,"yes"));
    Fp.append(createDimScalarProp("F"+phase_reac,"yes"));
    Ep.append(createDimScalarProp("E"+phase_reac,"yes"));
    mp.append(createDimScalarProp("m"+phase_reac,"yes"));
    np.append(createDimScalarProp("n"+phase_reac,"yes"));
    Tp.append(createDimScalarProp("T"+phase_reac,"yes"));
    hp.append(createDimScalarProp("h"+phase_reac,"yes"));
    piPyroReac_.append(createVolField<scalar>("piPyroReac[" + Foam::name(i) + "]",piTotal_));
    AField.append(createVolField<scalar>("AField[" + Foam::name(i) + "]",Ap[i]));
  }

  /* END: Pyrolysis reactions parameters */

  if (debug_) {
    Info << getTabLevel() << "debug: solidEps and solidRho  --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }

  // Initialize solidEps_, solidEpsI_, solidRho_, solidRhoI_
  // The solid phases are defined from value 1 (0 is attributed to the gas).
  for(int i = 0; i < nSolidPhases_; i++) {
    word name_epsI = "epsI_s["+Foam::name(i+1)+"]";
    word name_eps = "eps_s["+Foam::name(i+1)+"]";
    dimensionedScalar solidEpsI_value = createDimScalarProp("epsI[" + std::to_string(i+1) + "]",true,dimensionedScalar("0",dimless,0.0));
    solidEpsI_.append(createVolField<scalar>(name_epsI.c_str(),solidEpsI_value));
    solidEps_.append(createVolField<scalar>(name_eps.c_str(),solidEpsI_value));
    word name_rhoI = "rhoI_s["+Foam::name(i+1)+"]";
    word name_rho = "rho_s["+Foam::name(i+1)+"]";
    dimensionedScalar solidRhoI_value = createDimScalarProp("rhoI[" + std::to_string(i+1) + "]",true,dimensionedScalar("0",dimMass/dimVolume,0.0));
    solidRhoI_.append(createVolField<scalar>(name_rhoI.c_str(),solidRhoI_value));
    solidRho_.append(createVolField<scalar>(name_rho.c_str(),solidRhoI_value));
  }

  if (this->debug_) {
    Info << getTabLevel() << "debug: Xsi --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
    Info << getTabLevel() << "debug: solidEpsI_.size() = "<< solidEpsI_.size() <<" --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
    Info << getTabLevel() << "debug: solidRhoI_.size() = "<< solidRhoI_.size() <<" --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }


  // Initialize normFp
  normFp.resize(sizeXsi);
  forAll(normFp, i) {
    word name_normFp = "normFp["+Foam::name(i)+"]";
    normFp.set(i, createVolField<scalar>(name_normFp.c_str(),dimensionedScalar("0", dimless, scalar(0.0))));
  }

  // Initialize normFp, rho_c and rho_v
  initialize();

  if ( MaterialChemistryType_ == "ConstantFiniteRate" || MaterialChemistryType_ == "SpeciesConservation") {
    if (this->debug_) {
      Info << getTabLevel() << "debug: finiteRateSpeciesConservation --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
    }
    // Read pyrolysis-production SPECIES mass fractions - nb: only species present in the MaterialChemistry mechanism are looked up.
    const int speciesNames_size = materialChemistryModel.speciesNames().size();

    for (int i = 0; i < sizeXsi; i++) {
      Yp.append(new PList<scalar>());
      for (int j = 0; j < speciesNames_size; j++) {
        word nameYp = "gamma[" + std::to_string(hashXsi[i][0]) + "][" +
                      std::to_string(hashXsi[i][1]) + "][" +
                      speciesNames_[j] + "]";
        Yp[i].append(createScalarProp(nameYp,"yes"));
      }
    }
    pi_.resize(speciesNames_.size());
    forAll(pi_, specieI) {
      pi_.set(specieI,createVolField<scalar>("pi[" + speciesNames_[specieI] + "]",piTotal_));
    }
  }

  if (  MaterialChemistryType_ == "ConstantEquilibrium"
        || MaterialChemistryType_ == "EquilibriumElement" ) {
    //     Read pyrolysis-production ELEMENT mass fractions - nb: only elements present
    //     in the equilibrium file (Mutation) are looked up.
    if (this->debug_) {
      Info << getTabLevel() << "debug: ConstantEquilibrium || EquilibriumElement --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
    }
    const int ne_mix = elementNames_.size();
    for (int i = 0; i < sizeXsi; i++) {
      Zp.append(new PList<scalar>());
      for (int j = 0; j < ne_mix; j++) {
        word nameZp =  "zeta[" + std::to_string(hashXsi[i][0]) + "][" +
                       std::to_string(hashXsi[i][1]) + "][" +
                       elementNames_[j] + "]";
        Zp[i].append(createScalarProp(nameZp,"yes"));
      }
    }
    pi_.resize(elementNames_.size());
    forAll(pi_, elemI) {
      pi_.set(elemI,createVolField<scalar>("pi[" + elementNames_[elemI] + "]",piTotal_));
    }
  }

  if (this->debug) {
    Info << getTabLevel() << "debug: end --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  modelInitialized();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LinearArrheniusPyrolysisModel::~LinearArrheniusPyrolysisModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LinearArrheniusPyrolysisModel::update()
{
  //  The advancement of the pyrolysis reactions is computed
  // Update pyrolysis zone (nb: if Ta_<Tp[i], then there is no pyrolysis at all because the temperature is below the threshold)
  forAll(Xsi, i) {
    forAll(Ta_, cellI) {
      if (Ta_[cellI] < Tp[i].value()) {
        AField[i][cellI] = 0.0;
      } else {
        AField[i][cellI] = Ap[i].value();
      }
    }
  }
  if(this->debug_) {
    Info << getTabLevel() << "debug: solving Xsi --- Foam::linearArrheniusPyrolysisModel::solveXsi()" << endl;
  }
  for (int i = 0; i < Xsi.size(); i++) {
    volScalarField& Xsii = Xsi[i];
    if(!this->dynamicMesh_) {
      solve
      (
          fvm::ddt(Xsii)
          - pow(mag(1 - Xsii), mp[i]) * AField[i] * pow(Ta_, np[i]) * exp(-Ep[i] / (R_ * Ta_)),
          "Xsii"
      );
      Xsii.min(1.0);
      // Compute pyrolysis rate
      piTotal_ = piTotal_ * 0.;

      forAll(piPyroReac_, i) {
        piPyroReac_[i] =
            solidEpsI_[hashXsi[i][0]-1] *
            solidRhoI_[hashXsi[i][0]-1] *
            Fp[i] *
            (fvc::ddt(Xsi[i]));
        piPyroReac_[i].max(0);
        forAll(piPyroReac_[i], iCell) {
          if(piPyroReac_[i][iCell] < minPiPyro_) {
            piPyroReac_[i][iCell] = 0;
          }
        }
        forAll(piPyroReac_[i].boundaryField(), iPatch) {
          forAll(piPyroReac_[i].boundaryField()[iPatch], iFace) {
            if(piPyroReac_[i].boundaryField()[iPatch][iFace] < minPiPyro_) {
              piPyroReac_[i].boundaryFieldRef()[iPatch][iFace] = 0;
            }
          }
        }
        piTotal_ += piPyroReac_[i];
      }

    } else {
      solve
      (
          fvm::ddt(Xsii)
          - fvm::div(mesh_.phi(), Xsii)
          - pow(mag(1 - Xsii), mp[i]) * AField[i] * pow(Ta_, np[i]) * exp(-Ep[i] / (R_ * Ta_)),
          "Xsii"
      );
      Xsii.min(1.0);
      // Compute pyrolysis rate
      piTotal_ = piTotal_ * 0.;

      forAll(piPyroReac_, i) {
        piPyroReac_[i] =
            solidEpsI_[hashXsi[i][0]-1] *
            solidRhoI_[hashXsi[i][0]-1] *
            Fp[i] *
            (fvc::ddt(Xsi[i]) - fvc::div(mesh_.phi(), Xsi[i]));
        piPyroReac_[i].max(0);
        forAll(piPyroReac_[i], iCell) {
          if(piPyroReac_[i][iCell] < minPiPyro_) {
            piPyroReac_[i][iCell] = 0;
          }
        }
        forAll(piPyroReac_[i].boundaryField(), iPatch) {
          forAll(piPyroReac_[i].boundaryField()[iPatch], iFace) {
            if(piPyroReac_[i].boundaryField()[iPatch][iFace] < minPiPyro_) {
              piPyroReac_[i].boundaryFieldRef()[iPatch][iFace] = 0;
            }
          }
        }
        piTotal_ += piPyroReac_[i];
      }

    }
  }

  if (  MaterialChemistryType_ == "ConstantEquilibrium"
        || MaterialChemistryType_ == "EquilibriumElement" ) {

    forAll(pi_, elemI) {
      pi_[elemI]=  pi_[elemI]*0;
      forAll(piPyroReac_, i) {
        pi_[elemI] += Zp[i][elemI] * piPyroReac_[i] ;
      }
    }
  }
  if ( MaterialChemistryType_ == "ConstantFiniteRate"  || MaterialChemistryType_ == "SpeciesConservation" ) {
    forAll(pi_, specieI) {
      pi_[specieI]=  pi_[specieI]*0;
      forAll(piPyroReac_, i) {
        pi_[specieI] += Yp[i][specieI] * piPyroReac_[i] ;
      }
    }
  }

  // Solid densities update
  forAll(solidRho_, phaseI) {
    solidRho_[phaseI] = solidRhoI_[phaseI] * oneVolScalarField ;
    forAll(piPyroReac_, i) {
      if (hashXsi[i][0]-1 == phaseI) {
        solidRho_[phaseI] -= solidRhoI_[phaseI] * Fp[i] * Xsi[i];
      }
    }
  }
  // rho_s is updated in the material properties model
  if (Xsi.size()>0) {
    tau *= 0;
    forAll(Xsi, i) {
      tau += normFp[i] * (1. - Xsi[i]);       // average advancement of the pyrolysis reactions
    }
    tau.min(1.0);
    tau.max(0.0);
    tau.correctBoundaryConditions();
  }

  if (this->debug_) {
    forAll(tau, cellI) {
      if (tau[cellI]>1 || tau[cellI]<0) {
        FatalErrorInFunction << "Problem in pyrolysis model: tau[" << cellI << "] = " <<  tau[cellI] << exit(FatalError);
      }
    }
    forAll(tau.boundaryField(), patchI) {
      forAll(tau.boundaryField()[patchI], faceI) {
        if (tau.boundaryField()[patchI][faceI] >1 || tau.boundaryField()[patchI][faceI] <0) {
          FatalErrorInFunction << "Problem in pyrolysis model: tau.boundaryField()[" << patchI <<"][" << faceI << "] = " << tau.boundaryField()[patchI][faceI]  << exit(FatalError);
        }
      }
    }
  }
}

void Foam::LinearArrheniusPyrolysisModel::initialize()
{

  if (this->debug_) {
    Info << getTabLevel() << "debug: initialize normFp, rho_v, rho_c --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }

  normFp_sum*=0;
  forAll(Xsi, i) {
    normFp_sum +=
        Fp[i].value() *
        solidEpsI_[hashXsi[i][0]-1] *
        solidRhoI_[hashXsi[i][0]-1];
  }

  forAll(Xsi, i) {
    forAll(normFp_sum, cellI) {
      if (normFp_sum[cellI] > 0) {
        normFp[i][cellI] =
            Fp[i].value() *
            solidEpsI_[hashXsi[i][0]-1][cellI] *
            solidRhoI_[hashXsi[i][0]-1][cellI] /
            normFp_sum[cellI];
      } else {
        normFp[i][cellI] = 1;
      }
    }
    forAll(normFp_sum.boundaryField(), patchI) {
      forAll(normFp_sum.boundaryField()[patchI], faceI) {
        if (normFp_sum.boundaryField()[patchI][faceI] > 0) {
          normFp[i].boundaryFieldRef()[patchI][faceI]  =
              Fp[i].value() *
              solidEpsI_[hashXsi[i][0]-1].boundaryField()[patchI][faceI] *
              solidRhoI_[hashXsi[i][0]-1].boundaryField()[patchI][faceI] /
              normFp_sum.boundaryField()[patchI][faceI] ;
        } else {
          normFp[i].boundaryFieldRef()[patchI][faceI]  = 1;
        }
      }
    }
  }

  if (this->debug_) {
    Info << getTabLevel() << "debug: tau --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  if (Xsi.size()>0) {
    tau *= 0;
    forAll(Xsi, i) {
      tau += normFp[i] * (1 - Xsi[i]);
    }
  }

  if (this->debug_) {
    Info << getTabLevel() << "debug: rho_v, rho_c --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  rho_v_ *= 0;

  for (int i = 0; i < nSolidPhases_; i++) {
    rho_v_ += solidRhoI_[i] * solidEpsI_[i];
  }

  rho_c_ *= 0;
  for (int i = 0; i < nSolidPhases_; i++) {
    rho_c_ += solidRhoI_[i] * solidEpsI_[i];
  }
  forAll(Xsi, i) {
    rho_c_ -=
        Fp[i] *
        solidEpsI_[hashXsi[i][0]-1] *
        solidRhoI_[hashXsi[i][0]-1];
  }
}

// ************************************************************************* //

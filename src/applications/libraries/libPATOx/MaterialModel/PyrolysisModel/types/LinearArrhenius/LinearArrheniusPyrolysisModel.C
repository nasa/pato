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
    const word& dictName
)
  :
simplePyrolysisModel(mesh, dictName),
mesh_(mesh),
dictName_(dictName),
materialPropertiesModel(meshLookupOrConstructModel<simpleMaterialPropertiesModel>(mesh_,dictName_,"MaterialProperties")),
materialChemistryModel(meshLookupOrConstructModel<simpleMaterialChemistryModel>(mesh_,dictName_,"MaterialChemistry")),
nSolidPhases_(materialPropertiesModel.nSolidPhases()),
sizeXsi(0),
nPyroReac(nSolidPhases_),
hp(simplePyrolysisModel::hp_),
tau(simplePyrolysisModel::tau_),
solidEpsI(materialPropertiesModel.solidEpsI()),
solidRhoI(materialPropertiesModel.solidRhoI()),
Ta_(simplePyrolysisModel::Ta_),
R_(::constant::physicoChemical::R),
piTotal_(simplePyrolysisModel::piTotal_),
pi_(simplePyrolysisModel::pi_),
piPyroReac_(simplePyrolysisModel::piPyroReac_),
materialPropertiesDirectory(fileName(simplePyrolysisModel::materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand()),
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
MaterialChemistryType_(simplePyrolysisModel::materialDict_.isDict("MaterialChemistry")?word(simplePyrolysisModel::materialDict_.subDict("MaterialChemistry").lookup("MaterialChemistryType")):""),
solidRho_(materialPropertiesModel.solidRho()),
solidEps_(materialPropertiesModel.solidEps()),
rho_s_(simplePyrolysisModel::rho_s_),
rho_v_(simplePyrolysisModel::rho_v_),
rho_c_(simplePyrolysisModel::rho_c_),
oneVolScalarField
(
    IOobject
    (
        "one",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0)
),
minPiPyro_(simplePyrolysisModel::materialDict_.subDict("Pyrolysis").template lookupOrDefault("minPiPyro",1e-9))
{

  if (this->debug_) {
    Info << "--- start --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
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
    Info << "--- nPyroReac --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  forAll(nPyroReac, solidI) {
    nPyroReac.set
    (
        solidI,
        new scalar
        (
            constantPropertiesDictionary.lookupOrDefault<scalar>
            (
                "nPyroReac[" + std::to_string(solidI + 1) + "]",
                0.0
            )
        )
    );
    sizeXsi += nPyroReac[solidI];
  }
  hashXsi.setSize(sizeXsi, 2);
  int nindex = 0;

  if (this->debug_) {
    Info << "--- hashXsi --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  for (int solidI = 0; solidI < nSolidPhases_; solidI++) {
    for (int reacI = 0; reacI < nPyroReac[solidI]; reacI++) {
      hashXsi[nindex + reacI][0] = solidI + 1;
      hashXsi[nindex + reacI][1] = reacI + 1;
    }
    nindex += nPyroReac[solidI];
  }

  /* Pyrolysis advancements */
  Xsi.resize(sizeXsi);
  wordList zeroGradBC(mesh.boundaryMesh().size());
  forAll(zeroGradBC, reacI) {
    zeroGradBC = "zeroGradient";
  }

  if (this->debug_) {
    Info << "--- Xsi --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }

  forAll(Xsi, reacI) {
    Xsi.set
    (
        reacI,
        new volScalarField
        (
            IOobject
            (
                "Xsi[" + std::to_string(hashXsi[reacI][0]) + "][" +
                std::to_string(hashXsi[reacI][1]) + "]",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradBC
        )
    );
  }
  /* END:  Pyrolysis advancements  */

  if (this->debug_) {
    Info << "--- F,A,E,m,n,T,h --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  /* Pyrolysis reactions parameters */
  Fp.resize(sizeXsi);
  Ap.resize(sizeXsi);
  Ep.resize(sizeXsi);
  mp.resize(sizeXsi);
  np.resize(sizeXsi);
  Tp.resize(sizeXsi);
  hp.resize(sizeXsi);
  AField.resize(sizeXsi);
  piPyroReac_.resize(sizeXsi);

  forAll(Xsi, i) {
    Fp.set(i, new dimensionedScalar("1", dimless, 0.0));
    Ep.set(i, new dimensionedScalar("1", dimensionSet(1, 2, -2, 0, -1, 0, 0), 1.0));
    mp.set(i, new dimensionedScalar("1", dimless, 1.0));
    np.set(i, new dimensionedScalar("1", dimless, 1.0));
    Tp.set(i, new dimensionedScalar("1", dimensionSet(0, 0, 0, 1, 0, 0, 0), 10000.0));
    hp.set(i, new dimensionedScalar("1", dimensionSet(0, 2, -2, 0, 0, 0, 0), 1.0));
    piPyroReac_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "piPyroReac[" + std::to_string(i) + "]",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            piTotal_
        )
    );

  }

  forAll(Xsi, i) {
    Fp[i] = constantPropertiesDictionary.lookup("F[" + std::to_string(hashXsi[i][0]) + "][" + std::to_string(hashXsi[i][1]) + "]");
    Ep[i] = constantPropertiesDictionary.lookup("E[" + std::to_string(hashXsi[i][0]) + "][" + std::to_string(hashXsi[i][1]) + "]");
    mp[i] = constantPropertiesDictionary.lookup("m[" + std::to_string(hashXsi[i][0]) + "][" + std::to_string(hashXsi[i][1]) + "]");
    np[i] = constantPropertiesDictionary.lookup("n[" + std::to_string(hashXsi[i][0]) + "][" + std::to_string(hashXsi[i][1]) + "]");
    Tp[i] = constantPropertiesDictionary.lookup("T[" + std::to_string(hashXsi[i][0]) + "][" + std::to_string(hashXsi[i][1]) + "]");
    hp[i] = constantPropertiesDictionary.lookup("h[" + std::to_string(hashXsi[i][0]) + "][" + std::to_string(hashXsi[i][1]) + "]");
  }

  forAll(Xsi, i) {
    dimensionedScalar Ap_value = constantPropertiesDictionary.lookup("A[" + std::to_string(hashXsi[i][0]) + "][" + std::to_string(hashXsi[i][1]) + "]");
    Ap.set(i, new dimensionedScalar(Ap_value));
    AField.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "AField[" + std::to_string(i) + "]",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(Ap_value)
        )
    );
  }

  forAll(Xsi, i) {
    Info << Fp[i] <<  nl;
    Info << Ap[i] << nl;
    Info << Ep[i] << nl;
    Info << mp[i] << nl;
    Info << np[i] << nl;
    Info << Tp[i] << nl;
    Info << hp[i] << nl;
  }

  /* END: Pyrolysis reactions parameters */

  if (this->debug_) {
    Info << "--- normFp --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }

  normFp.resize(sizeXsi);
  forAll(normFp, i) {
    word name_normFp = "normFp["+Foam::name(i)+"]";
    normFp.set(i, new volScalarField
               (
                   IOobject
                   (
                       name_normFp.c_str(),
                       mesh_.time().timeName(),
                       mesh_,
                       IOobject::NO_READ,
                       IOobject::NO_WRITE
                   ),
                   mesh_,
                   dimensionedScalar("0", dimless, scalar(0.0))
               )
              );
  }
  volScalarField& normFp_sum(meshLookupOrConstructScalar(mesh,"normFp_sum", dimensionedScalar("0", dimMass/dimVolume, 0.0)));


  if (this->debug_) {
    Info << "--- Xsi --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
    Info << "--- solidEpsI.size() = "<< solidEpsI.size() <<" --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
    Info << "--- solidRhoI.size() = "<< solidRhoI.size() <<" --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;

  }

  forAll(Xsi, i) {
    normFp_sum +=
        Fp[i].value() *
        solidEpsI[hashXsi[i][0]-1] *
        solidRhoI[hashXsi[i][0]-1];
  }

  forAll(Xsi, i) {
    forAll(normFp_sum, cellI) {
      if (normFp_sum[cellI] > 0) {
        normFp[i][cellI] =
            Fp[i].value() *
            solidEpsI[hashXsi[i][0]-1][cellI] *
            solidRhoI[hashXsi[i][0]-1][cellI] /
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
              solidEpsI[hashXsi[i][0]-1].boundaryField()[patchI][faceI] *
              solidRhoI[hashXsi[i][0]-1].boundaryField()[patchI][faceI] /
              normFp_sum.boundaryField()[patchI][faceI] ;
        } else {
          normFp[i].boundaryFieldRef()[patchI][faceI]  = 1;
        }
      }
    }
  }

  if (this->debug_) {
    Info << "--- tau --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  if (Xsi.size()>0) {
    tau *= 0;
    forAll(Xsi, i) {
      tau += normFp[i] * (1 - Xsi[i]);
    }
  }


  if (this->debug_) {
    Info << "--- rho_v, rho_c --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
  rho_v_ *= 0;

  for (int i = 0; i < nSolidPhases_; i++) {
    rho_v_ += solidRhoI[i] * solidEpsI[i];
  }

  rho_c_ *= 0;
  for (int i = 0; i < nSolidPhases_; i++) {
    rho_c_ += solidRhoI[i] * solidEpsI[i];
  }
  forAll(Xsi, i) {
    rho_c_ -=
        Fp[i] *
        solidEpsI[hashXsi[i][0]-1] *
        solidRhoI[hashXsi[i][0]-1];
  }

  if ( MaterialChemistryType_ == "ConstantFiniteRate" || MaterialChemistryType_ == "SpeciesConservation") {
    if (this->debug_) {
      Info << "--- finiteRateSpeciesConservation --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
    }
    // Read pyrolysis-production SPECIES mass fractions - nb: only species present in the MaterialChemistry mechanism are looked up.
    const wordList& speciesNames_ = materialChemistryModel.speciesNames();
    Yp.setSize(sizeXsi, speciesNames_.size());

    for (int i = 0; i < sizeXsi; i++) {
      for (int j = 0; j < speciesNames_.size(); j++) {
        Yp[i][j] =
            constantPropertiesDictionary.lookupOrDefault<scalar>
            (
                "gamma[" + std::to_string(hashXsi[i][0]) + "][" +
                std::to_string(hashXsi[i][1]) + "][" +
                speciesNames_[j] + "]",
                0.0
            );
        if(Yp[i][j] != 0) {
          Info <<   "gamma[" << std::to_string(hashXsi[i][0]) << "][" <<
               std::to_string(hashXsi[i][1]) << "][" <<
               speciesNames_[j] << "]" << " = " << Yp[i][j] << endl;
        }
      }
    }
    pi_.resize(speciesNames_.size());
    forAll(pi_, specieI) {
      pi_.set
      (
          specieI,
          new volScalarField
          (
              IOobject
              (
                  "pi[" + speciesNames_[specieI] + "]",
                  mesh.time().timeName(),
                  mesh,
                  IOobject::READ_IF_PRESENT,
                  IOobject::NO_WRITE
              ),
              piTotal_
          )
      );

    }
  }

  if (  MaterialChemistryType_ == "ConstantEquilibrium"
        || MaterialChemistryType_ == "EquilibriumElement" ) {
    //     Read pyrolysis-production ELEMENT mass fractions - nb: only elements present
    //     in the equilibrium file (Mutation) are looked up.
    if (this->debug_) {
      Info << "--- ConstantEquilibrium || EquilibriumElement --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
    }
    const wordList& elementNames_= materialChemistryModel.elementNames();
    const int ne_mix = elementNames_.size();
    Zp.setSize(sizeXsi, ne_mix);
    for (int i = 0; i < sizeXsi; i++) {
      for (int j = 0; j < ne_mix; j++) {
        word nameZp =  "zeta[" + std::to_string(hashXsi[i][0]) + "][" +
                       std::to_string(hashXsi[i][1]) + "][" +
                       elementNames_[j] + "]";
        Zp[i][j] =
            constantPropertiesDictionary.lookupOrDefault<scalar>
            (
                nameZp,
                0.0
            );
        if (Zp[i][j] != 0) {
          Info << "zeta[" << std::to_string(hashXsi[i][0]) << "][" <<
               std::to_string(hashXsi[i][1]) << "][" <<
               elementNames_[j] << "]" << " = " << Zp[i][j] << endl;
        }
      }
    }
    pi_.resize(elementNames_.size());
    forAll(pi_, elemI) {
      pi_.set
      (
          elemI,
          new volScalarField
          (
              IOobject
              (
                  "pi[" + elementNames_[elemI] + "]",
                  mesh.time().timeName(),
                  mesh,
                  IOobject::READ_IF_PRESENT,
                  IOobject::NO_WRITE
              ),
              piTotal_
          )
      );

    }
  }

  if (this->debug) {
    Info << "--- end --- Foam::linearArrheniusPyrolysisModel::linearArrheniusPyrolysisSolver" << endl;
  }
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
    Info << "--- solving Xsi --- Foam::linearArrheniusPyrolysisModel::solveXsi()" << endl;
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
            solidEpsI[hashXsi[i][0]-1] *
            solidRhoI[hashXsi[i][0]-1] *
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
            solidEpsI[hashXsi[i][0]-1] *
            solidRhoI[hashXsi[i][0]-1] *
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
    solidRho_[phaseI] = solidRhoI[phaseI] * oneVolScalarField ;
    forAll(piPyroReac_, i) {
      if (hashXsi[i][0]-1 == phaseI) {
        solidRho_[phaseI] -= solidRhoI[phaseI] * Fp[i] * Xsi[i];
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


// ************************************************************************* //

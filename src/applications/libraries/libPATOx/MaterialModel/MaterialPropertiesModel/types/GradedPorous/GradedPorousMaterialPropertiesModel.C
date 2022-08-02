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
    along with OpenFOAM.  If Porous_t, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "GradedPorousMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GradedPorousMaterialPropertiesModel::GradedPorousMaterialPropertiesModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
simpleMaterialPropertiesModel(mesh, dictName),
materialPropertiesDirectory(fileName(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand()),
constantPropertiesDictionary
(
    IOobject
    (
        materialPropertiesDirectory+"/constantProperties",
        mesh.time().db().parent(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
nSolidPhases_(simpleMaterialPropertiesModel::nSolidPhases_=readScalar(constantPropertiesDictionary.lookup("nSolidPhases"))),
solidEpsI_(simpleMaterialPropertiesModel::solidEpsI_),
solidRhoI_(simpleMaterialPropertiesModel::solidRhoI_),
solidRho_(simpleMaterialPropertiesModel::solidRho_),
solidEps_(simpleMaterialPropertiesModel::solidEps_),
virginOrChar(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").template lookupOrDefault<word>("virginOrChar","char")),
multiRegions(true),
rho_s_(simpleMaterialPropertiesModel::rho_s_),
K_(simpleMaterialPropertiesModel::K_),
absorptivity_(simpleMaterialPropertiesModel::absorptivity_),
emissivity_(simpleMaterialPropertiesModel::emissivity_),
T_(meshLookupOrConstructScalar(mesh_, "Ta")),
p_(meshLookupOrConstructScalar(mesh_, "p")),
h_c_(simpleMaterialPropertiesModel::h_c_),
h_bar_(simpleMaterialPropertiesModel::h_bar_),
cp_(simpleMaterialPropertiesModel::cp_),
k_(simpleMaterialPropertiesModel::k_),
kijk_
(
    IOobject
    (
        "kijk",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("kI",dimensionSet(1, 1, -3, -1, 0, 0, 0),tensor(1, 0, 0, 0, 1, 0, 0, 0, 1))
),
k_abl_sym_
(
    IOobject
    (
        "k_abl_sym",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    symm(kijk_)

),
kiCoef_(constantPropertiesDictionary.lookupOrDefault<dimensionedScalar>("kiCoef", dimensionedScalar("0",dimensionSet(0,0,0,0,0),0))),
kjCoef_(constantPropertiesDictionary.lookupOrDefault<dimensionedScalar>("kjCoef", dimensionedScalar("0",dimensionSet(0,0,0,0,0),0))),
kkCoef_(constantPropertiesDictionary.lookupOrDefault<dimensionedScalar>("kkCoef", dimensionedScalar("0",dimensionSet(0,0,0,0,0),0))),
tP_(constantPropertiesDictionary.lookupOrDefault<dimensionedTensor>("tP",tensor(1,0,0,0,1,0,0,0,1))),
Kc_(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").found("K_c")?simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").lookup("K_c"):constantPropertiesDictionary.lookup("K_c")),
Kv_(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").found("K_v")?simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").lookup("K_v"):constantPropertiesDictionary.lookup("K_v")),
pyrolysisFlux_(simpleMaterialPropertiesModel::pyrolysisFlux_),
detailedSolidEnthalpies_(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").template lookupOrDefault<word>("detailedSolidEnthalpies","no")),
PyrolysisModel_(meshLookupOrConstructModel<simplePyrolysisModel>(mesh,dictName,"Pyrolysis")),
piPyroReac_(PyrolysisModel_.piPyroReac()),
piTotal_(PyrolysisModel_.piTotal()),
tau_(PyrolysisModel_.tau()),
hp_(PyrolysisModel_.hp()),
rho_v_(PyrolysisModel_.rho_v()),
rho_c_(PyrolysisModel_.rho_c()),
init_(init())
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GradedPorousMaterialPropertiesModel::~GradedPorousMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::GradedPorousMaterialPropertiesModel::update()
{
  if(simpleMaterialPropertiesModel::debug_) {
    Info << "--- update K --- Foam::porousMaterialPropertiesModel::update()" << endl;
  }
  // update permeability (used in MassModel)
  K_ = Kc_ + (Kv_ - Kc_) * tau_;

  if(simpleMaterialPropertiesModel::debug_) {
    Info << "--- update cp --- Foam::porousMaterialPropertiesModel::update()" << endl;
  }
  // update solid specific heat (used in EnergyModel)
  updateCp();

  if(simpleMaterialPropertiesModel::debug_) {
    Info << "--- update k --- Foam::porousMaterialPropertiesModel::update()" << endl;
  }
  // update solid thermal conductivity (used in EnergyModel)
  updatek();

  if(simpleMaterialPropertiesModel::debug_) {
    Info << "--- update h_bar --- Foam::porousMaterialPropertiesModel::update()" << endl;
  }
  // update averaged solid enthalpy (used in PyrolysisModel)
  updateHbar();

  if(simpleMaterialPropertiesModel::debug_) {
    Info << "--- update emissivity absorptivity --- Foam::porousMaterialPropertiesModel::update()" << endl;
  }
  // update emissivity and absorptivity
  updateEmissivityAbsorptivity();

  if(simpleMaterialPropertiesModel::debug_) {
    Info << "--- update h_c --- Foam::porousMaterialPropertiesModel::update()" << endl;
  }
  // update charred enthalpy
  updateHc();

  // update solid density
  rho_s_ *= 0;
  forAll(solidRho_, phaseI) {
    rho_s_ += solidRho_[phaseI] * solidEps_[phaseI];
  }
}

Switch Foam::GradedPorousMaterialPropertiesModel::init()
{
  if(this->debug_) {
    Info << "--- begin --- Foam::porousMaterialPropertiesModel::init()" << endl;
  }

  if(this->debug_) {
    Info << "--- nSolidPhases_=" << nSolidPhases_ << " --- Foam::porousMaterialPropertiesModel::init()" << endl;
  }

  if (!isFile(constantPropertiesDictionary.path()+"/constantProperties")) {
    FatalErrorInFunction << "Unknown materialPropertiesDirectory " << constantPropertiesDictionary.path()+"/constantProperties"
                         << exit(FatalError);

  }

  if (virginOrChar != "char" && virginOrChar != "virgin") {
    FatalErrorInFunction
        <<  "Problem in " << dictName_ << "MaterialProperties :" << nl << "\tvirginOrChar " << virginOrChar << ";" << nl << "\t"
        << "virginOrChar has to be \"virgin\" or \"char\" "
        << exit(FatalError);

  }

  if(this->debug_) {
    Info << "--- readEps() readRhoI() --- Foam::porousMaterialPropertiesModel::init()" << endl;
  }

  readEpsI();
  readRhoI();

  virginObject.reset(new LinearInterpolationMaterialPropertiesObject(readCharredVirginDict("virgin")));
  charObject.reset(new LinearInterpolationMaterialPropertiesObject(readCharredVirginDict("char")));

  if(this->debug_) {
    Info << "--- rho_s --- Foam::porousMaterialPropertiesModel::init()" << endl;
  }
  rho_s_*=0;
  forAll(solidRho_, phaseI) {
    rho_s_+= solidEps_[phaseI]*solidRho_[phaseI];
  }

  if(this->debug_) {
    Info << "--- noPyrolysisModel --- Foam::porousMaterialPropertiesModel::init()" << endl;
  }
  word typeNamePyrolysisModel = "no";
  if (this->materialDict_.isDict("Pyrolysis")) {
    typeNamePyrolysisModel = this->materialDict_.subDict("Pyrolysis").template lookupOrDefault<word>("PyrolysisType","no");
  }

  if(this->debug_) {
    Info << "--- solidRho_ --- Foam::porousMaterialPropertiesModel::init()" << endl;
  }

  if (typeNamePyrolysisModel=="no") {
    rho_v_ *= 0;
    forAll(solidRho_, phaseI) {
      rho_v_ += solidEpsI_[phaseI]*solidRhoI_[phaseI];
    }
    rho_c_=rho_v_;

    if (virginOrChar == "char") {
      forAll(tau_, cellI) {
        tau_[cellI]=0;
      }
      forAll(mesh_.boundaryMesh(), patchI) {
        forAll(tau_.boundaryFieldRef()[patchI], faceI) {
          tau_.boundaryFieldRef()[patchI][faceI] = 0;
        }
      }
    }
  }

  if(this->debug_) {
    Info << "--- end --- Foam::porousMaterialPropertiesModel::init()" << endl;
  }

  return true;
}

void Foam::GradedPorousMaterialPropertiesModel::updatek()
{
  LinearInterpolationMaterialPropertiesObject& virgin = virginObject();
  LinearInterpolationMaterialPropertiesObject& charred = charObject();
  forAll(kijk_, i) {
    kijk_[i].xx() =
        kiCoef_.value() *
        (
            virgin.ki(p_[i], T_[i]) *
            tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i])
            +
            charred.ki(p_[i], T_[i]) *
            (1. - tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i]))
        );
    kijk_[i].yy() =
        kjCoef_.value() *
        (
            virgin.kj(p_[i], T_[i]) *
            tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i])
            +
            charred.kj(p_[i], T_[i]) *
            (1. - tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i]))
        );
    kijk_[i].zz() =
        kkCoef_.value() *
        (
            virgin.kk(p_[i], T_[i]) *
            tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i])
            +
            charred.kk(p_[i], T_[i]) *
            (1. - tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i]))
        );
  }

  forAll(mesh_.boundaryMesh(), patchI) {
    forAll(kijk_.boundaryFieldRef()[patchI], faceI) {
      kijk_.boundaryFieldRef()[patchI][faceI].xx() =
          kiCoef_.value() *
          (
              virgin.ki(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
              tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI])
              +
              charred.ki(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
              (1. - tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI]))
          );
      kijk_.boundaryFieldRef()[patchI][faceI].yy() =
          kjCoef_.value() *
          (
              virgin.kj(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
              tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI])
              +
              charred.kj(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
              (1. - tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI]))
          );
      kijk_.boundaryFieldRef()[patchI][faceI].zz() =
          kkCoef_.value() *
          (
              virgin.kk(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
              tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI])
              +
              charred.kk(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
              (1. - tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI]))
          );
    }
  }


  k_ = tP_.value() & kijk_ & tP_.value().T(); // passage from basis ijk (main directions of the conductivity tensor) to basis xyz (mesh)
  if (multiRegions) {
    k_abl_sym_ = symm(k_);
  }
}

void Foam::GradedPorousMaterialPropertiesModel::updateCp()
{
  LinearInterpolationMaterialPropertiesObject& virgin = virginObject();
  LinearInterpolationMaterialPropertiesObject& charred = charObject();
  forAll(cp_, i) {
    cp_[i] = virgin.cp(p_[i], T_[i]) *
             tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i]) +
             charred.cp(p_[i], T_[i]) *
             (1. - tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i]));
  }

  forAll(cp_.boundaryFieldRef(), patchI) {
    forAll(cp_.boundaryFieldRef()[patchI], faceI) {

      cp_.boundaryFieldRef()[patchI][faceI]=
          virgin.cp(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
          tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI])
          +
          charred.cp(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
          (1. - tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI]));
    }

  }

}

void Foam::GradedPorousMaterialPropertiesModel::updateHbar()
{
  LinearInterpolationMaterialPropertiesObject& virgin = virginObject();
  LinearInterpolationMaterialPropertiesObject& charred = charObject();

  if(!detailedSolidEnthalpies_) {
    forAll (h_bar_, i) {
      if (rho_v_[i] == rho_c_[i]) {
        h_bar_[i] = virgin.h(p_[i], T_[i]);
      } else {
        h_bar_[i] =
            (rho_v_[i] * virgin.h(p_[i], T_[i]) - rho_c_[i] * charred.h(p_[i], T_[i])) /
            (rho_v_[i] - rho_c_[i]);
      }
    }
    forAll(mesh_.boundaryMesh(), patchI) {
      forAll (h_bar_.boundaryFieldRef()[patchI], faceI) {
        if (rho_v_.boundaryField()[patchI][faceI] == rho_c_.boundaryField()[patchI][faceI]) {
          h_bar_.boundaryFieldRef()[patchI][faceI] = virgin.h(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]);
        } else {
          h_bar_.boundaryFieldRef()[patchI][faceI] =
              (
                  rho_v_.boundaryField()[patchI][faceI] * virgin.h(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) -
                  rho_c_.boundaryField()[patchI][faceI] * charred.h(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI])
              ) /
              (rho_v_.boundaryField()[patchI][faceI] - rho_c_.boundaryField()[patchI][faceI]);
        }
      }
    }
    pyrolysisFlux_ = - h_bar_ * piTotal_;    // NASA Apollo approach (Kendall 1960)
  } else {

    if (piPyroReac_.size()!=0 && hp_.size()!=0) {

      pyrolysisFlux_ = pyrolysisFlux_ * 0.; // Multiphase approach (Lachaud 2017)
      tmp<volScalarField> hs_tmp = h_bar_ * 0.;
      volScalarField& hs_ = const_cast<volScalarField&>(hs_tmp());

      forAll (h_bar_, i) {
        if (rho_v_[i] == rho_c_[i]) {
          hs_[i] = virgin.hs(p_[i], T_[i]) - virgin.hs(p_[i], 298);
        } else {
          hs_[i] =
              (
                  rho_v_[i] * (virgin.hs(p_[i], T_[i]) - virgin.hs(p_[i], 298)) -
                  rho_c_[i] * (charred.hs(p_[i], T_[i]) - charred.hs(p_[i], 298))
              ) /
              (rho_v_[i] - rho_c_[i]);
        }
      }

      forAll(mesh_.boundaryMesh(), patchI) {
        forAll (h_bar_.boundaryFieldRef()[patchI], faceI) {
          if (rho_v_.boundaryField()[patchI][faceI] == rho_c_.boundaryField()[patchI][faceI]) {
            hs_.boundaryFieldRef()[patchI][faceI] = virgin.hs(p_.boundaryFieldRef()[patchI][faceI], T_.boundaryFieldRef()[patchI][faceI]) - virgin.hs(p_.boundaryFieldRef()[patchI][faceI], 298);
          } else {
            hs_.boundaryFieldRef()[patchI][faceI] =
                (
                    rho_v_.boundaryField()[patchI][faceI] * (virgin.hs(p_.boundaryFieldRef()[patchI][faceI], T_.boundaryFieldRef()[patchI][faceI]) - virgin.hs(p_.boundaryFieldRef()[patchI][faceI], 298)) -
                    rho_c_.boundaryField()[patchI][faceI] * (charred.hs(p_.boundaryFieldRef()[patchI][faceI], T_.boundaryFieldRef()[patchI][faceI]) - charred.hs(p_.boundaryFieldRef()[patchI][faceI], 298))
                ) /
                (rho_v_.boundaryField()[patchI][faceI] - rho_c_.boundaryField()[patchI][faceI]);
          }
        }
      }

      forAll(piPyroReac_, i) {
        pyrolysisFlux_ += - piPyroReac_[i] * (hp_[i] + hs_);
      }
    }
  }
}

void Foam::GradedPorousMaterialPropertiesModel::updateHc()
{
  LinearInterpolationMaterialPropertiesObject& charred = charObject();
  forAll(h_c_, cellI) {
    h_c_[cellI] = charred.h(p_[cellI], T_[cellI]);
  }


  forAll(mesh_.boundaryMesh(), patchI) {
    fvPatchScalarField& Tp = T_.boundaryFieldRef()[patchI];
    fvPatchScalarField& pp = p_.boundaryFieldRef()[patchI];
    fvPatchScalarField& rho_sp = rho_s_.boundaryFieldRef()[patchI];
    fvPatchScalarField& h_cp = h_c_.boundaryFieldRef()[patchI];

    forAll(h_cp, faceI) {
      h_cp[faceI] = charred.h(pp[faceI], Tp[faceI]);
    }
  }
}

void Foam::GradedPorousMaterialPropertiesModel::updateEmissivityAbsorptivity()
{
  LinearInterpolationMaterialPropertiesObject& virgin = virginObject();
  LinearInterpolationMaterialPropertiesObject& charred = charObject();

  forAll(mesh_.boundaryMesh(), patchI) {
    fvPatchScalarField& Tp = T_.boundaryFieldRef()[patchI];
    fvPatchScalarField& emissivityp = emissivity_.boundaryFieldRef()[patchI];
    fvPatchScalarField& absorptivityp = absorptivity_.boundaryFieldRef()[patchI];
    fvPatchScalarField& pp = p_.boundaryFieldRef()[patchI];
    fvPatchScalarField& rho_sp = rho_s_.boundaryFieldRef()[patchI];
    fvPatchScalarField& taup = tau_.boundaryFieldRef()[patchI];

    forAll(emissivityp, faceI) {
      emissivityp[faceI] = virgin.eps(pp[faceI], Tp[faceI]) * taup[faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_sp[faceI], rho_c_.boundaryField()[patchI][faceI])
                           + charred.eps(pp[faceI], Tp[faceI]) * (1 - taup[faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_sp[faceI], rho_c_.boundaryField()[patchI][faceI]));
    }

    forAll(absorptivityp, faceI) {
      absorptivityp[faceI] = virgin.alpha(pp[faceI], Tp[faceI]) * taup[faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_sp[faceI], rho_c_.boundaryField()[patchI][faceI])
                             + charred.alpha(pp[faceI], Tp[faceI]) * (1 - taup[faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_sp[faceI], rho_c_.boundaryField()[patchI][faceI]));
    }
  }
}

void Foam::GradedPorousMaterialPropertiesModel::readEpsI()
{
  fileName constPropFile = materialPropertiesDirectory+"/constantProperties";
  if (isFile(constPropFile)) {
    Info << "Reading " << constPropFile << endl;
  }

  if (nSolidPhases_ < 1) {
    FatalErrorInFunction
        << "The number of solid phases 'nSolidPhases' must be > 0. "
        "If not, just use a CFD solver!" << nl
        << exit(FatalError);
  }

  solidEpsI_.resize(nSolidPhases_);
  solidEps_.resize(nSolidPhases_);
  forAll(solidEpsI_, i) {
    word name_epsI = "epsI_s_phase"+Foam::name(i+1);
    word name_eps = "eps_s_phase"+Foam::name(i+1);
    solidEpsI_.set(
        i,
        new volScalarField
        (
            IOobject
            (
                name_epsI.c_str(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("epsI", dimless, scalar(0.0))
        )
    );
    solidEps_.set(
        i,
        new volScalarField
        (
            IOobject
            (
                name_eps.c_str(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("eps", dimless, scalar(0.0))
        )
    );
  }

  gradingFiles_.resize(nSolidPhases_);
  gradingPatchNames_.resize(nSolidPhases_);
  gradingPatchID_.resize(nSolidPhases_);

  for (int i = 0; i < nSolidPhases_; i++) { // the solid phases are defined from value 1 (0 is attributed to the gas).
    gradingFiles_.set(i,new fileName(constantPropertiesDictionary.lookupOrDefault<fileName>("gradingFileEpsI["+std::to_string(i+1)+"]","")));
    if(gradingFiles_[i]=="") {
      dimensionedScalar solidEpsI_value(constantPropertiesDictionary.lookup("epsI[" + std::to_string(i+1) + "]"));
      solidEpsI_[i] = solidEpsI_value;
      solidEps_[i] = solidEpsI_value;
      Info << "epsI[" << i+1 << "]=" << solidEpsI_value.value() <<  nl;
    } else {
      // Name of the patches for the grading distance calculation
      gradingPatchNames_.set(i, new word(constantPropertiesDictionary.lookup("gradingPatchNameEpsI["+std::to_string(i+1)+"]")));
      label patchID=mesh_.boundaryMesh().findPatchID(gradingPatchNames_[i]);
      if (patchID<0) {
        FatalErrorInFunction << "The grading patch " << gradingPatchNames_[i] << " is not found in mesh." << exit(FatalError);
      }
      gradingPatchID_.set(i, new label(patchID)); // ID of the grading patches

      List<scalarList> data_(readFileData(gradingFiles_[i])); // distance and phase volume fraction

      if(data_.size()!=2) {
        FatalErrorInFunction << gradingFiles_[i] << " must have two columns (distance and phase volume fraction)" << exit(FatalError);
      }
      if(data_[0].size()!=data_[1].size()) {
        FatalErrorInFunction << gradingFiles_[i] << " must have two columns (distance and phase volume fraction) of the same length" << exit(FatalError);
      }

      List<vector> global_Cf = globalFaceCenters(mesh_, gradingPatchNames_[i], this->dictName_); // global face centers

      forAll(solidEpsI_[i], cellI) {
        scalar distance_=-1;
        forAll(global_Cf, faceI) {
          scalar d = dist(mesh_.C()[cellI], global_Cf[faceI]);
          if (d<distance_ || distance_<0) {
            distance_=d;
          }
        }
        if ( (distance_>max(data_[0])) || (distance_ < (min(data_[0])))) {
          solidEpsI_[i][cellI] = 0;
        } else {
          solidEpsI_[i][cellI] = linearInterpolation(data_[0],data_[1],distance_);
        }
      }

      forAll(solidEpsI_[i].boundaryField(), patchI) {
        forAll(solidEpsI_[i].boundaryField()[patchI], faceI) {
          scalar distance_=-1;
          if (patchI==gradingPatchID_[i]) {
            distance_=0;
          } else {
            forAll(global_Cf, faceJ) {
              scalar d = dist(mesh_.Cf().boundaryField()[patchI][faceI], global_Cf[faceJ]);
              if (d<distance_ || distance_<0) {
                distance_=d;
              }
            }
          }
          if ( (distance_>max(data_[0])) || (distance_ < (min(data_[0])))) {
            solidEpsI_[i].boundaryFieldRef()[patchI][faceI] = 0;
          } else {
            solidEpsI_[i].boundaryFieldRef()[patchI][faceI] = linearInterpolation(data_[0],data_[1],distance_);
          }
        }
      }
      forAll(solidEps_[i], cellI) {
        solidEps_[i][cellI]=solidEpsI_[i][cellI];
      }
      forAll(solidEps_[i].boundaryField(), patchI) {
        forAll(solidEps_[i].boundaryField()[patchI], faceI) {
          solidEps_[i].boundaryFieldRef()[patchI][faceI]=solidEpsI_[i].boundaryField()[patchI][faceI];
        }
      }

      Info << "epsI[" << i+1 << "] graded using " << gradingFiles_[i] <<  nl;
    }
  }
}

void Foam::GradedPorousMaterialPropertiesModel::readRhoI()
{
  fileName constPropFile = materialPropertiesDirectory+"/constantProperties";
  if (isFile(constPropFile)) {
    Info << "Reading " << constPropFile << endl;
  }

  if (nSolidPhases_ < 1) {
    FatalErrorInFunction
        << "The number of solid phases 'nSolidPhases' must be > 0. "
        "If not, just use a CFD solver!" << nl
        << exit(FatalError);
  }

  if(this->debug_) {
    Info << "--- solidRhoI_ solidRho_ --- Foam::porousMaterialPropertiesModel::readRhoI()" << endl;
  }

  solidRhoI_.resize(nSolidPhases_);
  solidRho_.resize(nSolidPhases_);
  forAll(solidRhoI_, i) {
    word name_rhoI = "rhoI_s_phase"+Foam::name(i+1);
    word name_rho = "rho_s_phase"+Foam::name(i+1);
    solidRhoI_.set(
        i,
        new volScalarField
        (
            IOobject
            (
                name_rhoI.c_str(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("rhoI", dimMass/dimVolume, scalar(0.0))
        )
    );
    solidRho_.set(
        i,
        new volScalarField
        (
            IOobject
            (
                name_rho.c_str(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("rho", dimMass/dimVolume, scalar(0.0))
        )
    );
  }

  for (int i = 0; i < nSolidPhases_; i++) { // the solid phases are defined from value 1 (0 is attributed to the gas).
    dimensionedScalar solidRhoI_value = constantPropertiesDictionary.lookup("rhoI[" + std::to_string(i+1) + "]");
    solidRhoI_[i] =  solidRhoI_value;
    solidRho_[i] = solidRhoI_value;
    Info << "rhoI[" << i+1 << "]=" << solidRhoI_value.value() <<  nl;
  }
}

LinearInterpolationMaterialPropertiesObject Foam::GradedPorousMaterialPropertiesModel::readCharredVirginDict(word charredVirgin)
{
  // read and store the material property table into the RAM for faster access
  Info << "Reading " << charredVirgin << " material properties" << nl;

  fileName charredVirginFileName = materialPropertiesDirectory+"/"+charredVirgin;
  IFstream charredVirginThermo(charredVirginFileName);  // opens an input file
  IFstream charredVirginThermoTemp(charredVirginFileName);  // opens an input file

  if (charredVirginThermo.good()==false) { // checks the input file is opened
    FatalErrorInFunction
        << charredVirginFileName << " not found" << nl
        << exit(FatalError); // exits otherwise
  }

  //Info << "The " << charredVirgin << " material table must have 9 columns: p, T, cp, h, ki, kj, kk, eps, alpha." << nl;
  int columnTableV=9;
  int rawTableV=0;
  scalar tempV, rawTableFracV, rawTableIntV;
  int i_rawV=0;
  while(true) {
    charredVirginThermoTemp >> tempV;
    if(charredVirginThermoTemp.eof()==1) {
      break;
    }
    i_rawV++;
  }
  rawTableFracV=modf(static_cast<scalar>(i_rawV) / static_cast<scalar>(columnTableV), &rawTableIntV);
  if (rawTableFracV!=0) {
    FatalErrorInFunction <<
                         rawTableIntV << " lines and " << rawTableFracV * columnTableV << " 'scalar' have been read."
                         << "... and it does not seem to have 9 columns." << nl
                         << exit(FatalError); // exits
  } else {
    rawTableV = i_rawV / columnTableV;
    Info << "The " << charredVirgin << " material table has " << rawTableV << " lines."<< nl;
  }


  RectangularMatrix<scalar> charredVirginThermoTable(rawTableV, columnTableV+1);
  for(int x=0; x<rawTableV; x++) {
    for(int i=0; i<columnTableV; i++) {
      charredVirginThermo >> charredVirginThermoTable[x][i];
    }
  }

  // We precompute and add the sensible enthalpy as a last column hs = int(cp(T)dT)
  for(int x=0; x<rawTableV; x++) {
    charredVirginThermoTable[x][9]=0;
  }
  Info << charredVirginThermoTable[0][1] << " " << charredVirginThermoTable[0][2] << " " << charredVirginThermoTable[0][9] << endl;
  for(int x=1; x<rawTableV; x++) {

    if (charredVirginThermoTable[x][0] != charredVirginThermoTable[x-1][0]) {
      charredVirginThermoTable[x][9] = 0;
    } else {
      charredVirginThermoTable[x][9] += charredVirginThermoTable[x-1][9]
                                        + (charredVirginThermoTable[x][2] + charredVirginThermoTable[x-1][2])/2
                                        * (charredVirginThermoTable[x][1] - charredVirginThermoTable[x-1][1]);
    }
  }

  // charredVirgin is an object of the class materialProperties, which handles the material properties interpolations when called
  LinearInterpolationMaterialPropertiesObject objectCharredVirgin(charredVirginThermoTable);
  return objectCharredVirgin;
}

// ************************************************************************* //

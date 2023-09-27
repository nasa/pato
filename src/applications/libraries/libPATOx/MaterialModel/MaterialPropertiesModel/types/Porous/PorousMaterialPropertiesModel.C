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

#include "PorousMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PorousMaterialPropertiesModel::PorousMaterialPropertiesModel
(
    const fvMesh& mesh,
    const word& regionName
)
  :
simpleMaterialPropertiesModel(mesh, regionName),
CommonMaterialPropertiesModel(*this, modelName,"Porous"),
massModel_(refModel<simpleMassModel>()),
p_(massModel_.refVolField<scalar>("p")),
energyModel_(refModel<simpleEnergyModel>()),
T_(energyModel_.refVolField<scalar>("Ta")),
rho_s_(energyModel_.refVolField<scalar>("rho_s")),
pyrolysisModel_(refModel<simplePyrolysisModel>()),
tau_(pyrolysisModel_.refVolField<scalar>("tau")),
rho_v_(pyrolysisModel_.refVolField<scalar>("rho_v")),
rho_c_(pyrolysisModel_.refVolField<scalar>("rho_c"))
{
  listMapFields= {
    newMapField("K","volTensorField",initTensorLinear,updateTensorLinear),
    newMapField("Beta","volTensorField",initTensorLinear,updateTensorLinear),
    newMapField("nu","volScalarField",initScalarLinear,updateScalarLinear),
    newMapField("E","volScalarField",initScalarLinear,updateScalarLinear),
    newMapField("alpha","volScalarField",initScalarLinear,updateScalarLinear),
    newMapField("xi","volScalarField",initScalarConstant,updateScalarConstant),
    newMapField("cp","volScalarField",initVirginCharFiles,update_cp),
    newMapField("k","volTensorField",initVirginCharFiles,updatek),
    newMapField("pyrolysisFlux","volScalarField",initVirginCharFiles,updatePyrolysisFlux),
    newMapField("emissivity","volScalarField",initVirginCharFiles,updateEmissivityAbsorptivity),
    newMapField("absorptivity","volScalarField",initVirginCharFiles,updateEmissivityAbsorptivity),
    newMapField("h_c","volScalarField",initVirginCharFiles,updateHc),
    newMapField("rho_s","volScalarField",init_rho_s,update_rho_s),
  };
  modelInitialized();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PorousMaterialPropertiesModel::~PorousMaterialPropertiesModel()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PorousMaterialPropertiesModel::initialize()
{
  CommonMaterialPropertiesModel::initialize();
}

void Foam::PorousMaterialPropertiesModel::update()
{
  CommonMaterialPropertiesModel::update();
}

void Foam::PorousMaterialPropertiesModel::init_rho_s(word name)
{
  solidRho_ = pyrolysisModel_.refVolFieldPList<scalar>(name,1);
  solidEps_ = pyrolysisModel_.refVolFieldPList<scalar>("eps_s",1);
}

void Foam::PorousMaterialPropertiesModel::initVirginCharFiles(word name)
{
  // Columns of virgin/char files = {"p","Ta","cp","h","ki","kj","kk","emissivity","absorptivity"}
  const wordList availFieldNames( {"cp","h_c","pyrolysisFlux","k","emissivity","absorptivity"
  });
  if (!foundInList(name,availFieldNames)) {
    FatalError << "initVirginCharFiles(): " << name
               << " not available. Only the following fields are allowed: "
               << availFieldNames << exit(FatalError);
  }
  if (!initVirginCharFiles_) {
    // virgin/char objects
    virginObject.reset(new LinearInterpolationMaterialPropertiesObject(readCharredVirginDict("virgin")));
    charObject.reset(new LinearInterpolationMaterialPropertiesObject(readCharredVirginDict("char")));
    initVirginCharFiles_="yes";
  }
  if(name=="k") {
    dimensionedTensor kI = dimensionedTensor("kI",dimensionSet(1, 1, -3, -1, 0, 0, 0),I_);
    tensorFields_.insert("kijk",createVolField<tensor>("kijk",kI));
    symmTensorFields_.insert("k_abl_sym",createVolField<symmTensor>("k_abl_sym",symm(kI)));
    wordList dir= {"i","j","k"};
    forAll(dir,i) {
      word nameDir = "k"+dir[i]+"Coef";
      dimScalarMatProps_.insert(nameDir,createDimScalarProp(nameDir,"yes",dimensionedScalar("0",dimless,1)));
    }
    dimTensorMatProps_.insert("tP",createDimTensorProp("tP","yes",dimensionedTensor("tP",dimless,I_)));
  }
  if(name=="pyrolysisFlux") {
    scalarFields_.insert("h_bar",createVolField<scalar>("h_bar",dimensionedScalar("0",pow(dimLength,2)/pow(dimTime,2),scalar(0))));
    scalarFields_.insert("piTotal",pyrolysisModel_.refVolField<scalar>("piTotal"));
    switchMatProps_.insert("detailedSolidEnthalpies",createSwitchProp("detailedSolidEnthalpies","yes",Switch("no")));
    if (mesh_.objectRegistry::foundObject<volScalarField>("piPyroReac[0]")) {
      piPyroReac_=pyrolysisModel_.refVolFieldPList<scalar>("piPyroReac");
    }
  }
}

void Foam::PorousMaterialPropertiesModel::initScalarConstant(word name)
{
  dimScalarMatProps_.insert(name,createDimScalarProp(name));
}

void Foam::PorousMaterialPropertiesModel::initScalarLinear(word name)
{
  dimScalarMatProps_.insert(name+"_v",createDimScalarProp(name+"_v")); // virgin
  dimScalarMatProps_.insert(name+"_c",createDimScalarProp(name+"_c")); // char
}

void Foam::PorousMaterialPropertiesModel::initVectorLinear(word name)
{
  dimVectorMatProps_.insert(name+"_v",createDimVectorProp(name+"_v"));
  dimVectorMatProps_.insert(name+"_c",createDimVectorProp(name+"_c"));
}

void Foam::PorousMaterialPropertiesModel::initTensorLinear(word name)
{
  dimTensorMatProps_.insert(name+"_v",createDimTensorProp(name+"_v"));
  dimTensorMatProps_.insert(name+"_c",createDimTensorProp(name+"_c"));
}

void Foam::PorousMaterialPropertiesModel::updateScalarConstant(word name)
{
  // Update: field = field_prop;
  volScalarField& field = scalarFields_[name];
  field=dimScalarMatProps_[name];
}

void Foam::PorousMaterialPropertiesModel::updateScalarLinear(word name)
{
  // Update: field = field_c + (field_v - field_c) * tau_;
  volScalarField& field = scalarFields_[name];
  field=dimScalarMatProps_[name+"_c"]+(dimScalarMatProps_[name+"_v"]-dimScalarMatProps_[name+"_c"])*tau_;
}

void Foam::PorousMaterialPropertiesModel::updateVectorLinear(word name)
{
  // Update: field = field_c + (field_v - field_c) * tau_;
  volVectorField& field = vectorFields_[name];
  field=dimVectorMatProps_[name+"_c"]+(dimVectorMatProps_[name+"_v"]-dimVectorMatProps_[name+"_c"])*tau_;
}

void Foam::PorousMaterialPropertiesModel::updateTensorLinear(word name)
{
  // Update: field = field_c + (field_v - field_c) * tau_;
  volTensorField& field = tensorFields_[name];
  field=dimTensorMatProps_[name+"_c"]+(dimTensorMatProps_[name+"_v"]-dimTensorMatProps_[name+"_c"])*tau_;
}

void Foam::PorousMaterialPropertiesModel::update_cp(word name)
{
  LinearInterpolationMaterialPropertiesObject& virgin = virginObject();
  LinearInterpolationMaterialPropertiesObject& charred = charObject();

  volScalarField& cp_ = scalarFields_[name];

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

void Foam::PorousMaterialPropertiesModel::update_rho_s(word name)
{
  volScalarField& rho_s_ = scalarFields_[name];
  rho_s_ *= 0;
  forAll(solidRho_, phaseI) {
    rho_s_ += solidRho_[phaseI] * solidEps_[phaseI];
  }
}

void Foam::PorousMaterialPropertiesModel::updatek(word name)
{
  LinearInterpolationMaterialPropertiesObject& virgin = virginObject();
  LinearInterpolationMaterialPropertiesObject& charred = charObject();
  volTensorField& k_ = tensorFields_[name];
  volSymmTensorField& k_abl_sym_ = symmTensorFields_["k_abl_sym"];
  volTensorField& kijk_ = tensorFields_["kijk"];
  const dimensionedScalar& kiCoef_ = dimScalarMatProps_["kiCoef"];
  const dimensionedScalar& kjCoef_ = dimScalarMatProps_["kjCoef"];
  const dimensionedScalar& kkCoef_ = dimScalarMatProps_["kkCoef"];
  const dimensionedTensor& tP_ = dimTensorMatProps_["tP"];

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
              virgin.ki(p_.boundaryField()[patchI][faceI],
                        T_.boundaryField()[patchI][faceI])
              *tau_.boundaryField()[patchI][faceI]
              *rho_v_.boundaryField()[patchI][faceI]
              /max
              (
                  rho_s_.boundaryField()[patchI][faceI],
                  rho_c_.boundaryField()[patchI][faceI]
              )
              + charred.ki
              (
                  p_.boundaryField()[patchI][faceI],
                  T_.boundaryField()[patchI][faceI]
              )
              *(1. - tau_.boundaryField()[patchI][faceI]
                *rho_v_.boundaryField()[patchI][faceI]
                /max
                (
                    rho_s_.boundaryField()[patchI][faceI],
                    rho_c_.boundaryField()[patchI][faceI]
                ))
          );
      kijk_.boundaryFieldRef()[patchI][faceI].yy() =
          kjCoef_.value() *
          (
              virgin.kj
              (
                  p_.boundaryField()[patchI][faceI],
                  T_.boundaryField()[patchI][faceI]
              )
              *tau_.boundaryField()[patchI][faceI]
              *rho_v_.boundaryField()[patchI][faceI]
              /max
              (
                  rho_s_.boundaryField()[patchI][faceI],
                  rho_c_.boundaryField()[patchI][faceI]
              )
              + charred.kj
              (
                  p_.boundaryField()[patchI][faceI],
                  T_.boundaryField()[patchI][faceI]
              )
              *(1. - tau_.boundaryField()[patchI][faceI]
                *rho_v_.boundaryField()[patchI][faceI]
                /max
                (
                    rho_s_.boundaryField()[patchI][faceI],
                    rho_c_.boundaryField()[patchI][faceI])
               )
          );
      kijk_.boundaryFieldRef()[patchI][faceI].zz() =
          kkCoef_.value() *
          (
              virgin.kk(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI])
              *tau_.boundaryField()[patchI][faceI]
              *rho_v_.boundaryField()[patchI][faceI]
              /max
              (
                  rho_s_.boundaryField()[patchI][faceI],
                  rho_c_.boundaryField()[patchI][faceI]
              )
              + charred.kk
              (
                  p_.boundaryField()[patchI][faceI],
                  T_.boundaryField()[patchI][faceI]
              )
              *(1. - tau_.boundaryField()[patchI][faceI]
                *rho_v_.boundaryField()[patchI][faceI]
                /max
                (
                    rho_s_.boundaryField()[patchI][faceI],
                    rho_c_.boundaryField()[patchI][faceI]
                ))
          );
    }
  }

  // passage from basis ijk (main directions of the conductivity tensor) to basis xyz (mesh)
  k_.ref() = tP_.value() & kijk_() & tP_.value().T();
  forAll(mesh_.boundaryMesh(), patchI) {
    if (mesh_.boundaryMesh()[patchI].type()=="empty")
      continue;
    forAll(mesh_.boundaryMesh()[patchI], faceI) {
      k_.boundaryFieldRef()[patchI][faceI] = tP_.value() & kijk_.boundaryField()[patchI][faceI] & tP_.value().T();
    }
  }

  // symmetric tensor
  k_abl_sym_ = symm(k_);
}

void Foam::PorousMaterialPropertiesModel::updatePyrolysisFlux(word name)
{
  LinearInterpolationMaterialPropertiesObject& virgin = virginObject();
  LinearInterpolationMaterialPropertiesObject& charred = charObject();

  volScalarField& pyrolysisFlux_=scalarFields_[name];
  volScalarField& h_bar_=scalarFields_["h_bar"];
  volScalarField& piTotal_=scalarFields_["piTotal"];
  const PList<dimensionedScalar>& hp_=pyrolysisModel_.hp();
  const Switch& detailedSolidEnthalpies_=switchMatProps_["detailedSolidEnthalpies"];
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
          hs_[i] = virgin.hs(p_[i], T_[i]) - virgin.hs(p_[i], 298); // T_0 instead of 298?
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
            hs_.boundaryFieldRef()[patchI][faceI] = virgin.hs(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) - virgin.hs(p_.boundaryField()[patchI][faceI], 298);
          } else {
            hs_.boundaryFieldRef()[patchI][faceI] =
                (
                    rho_v_.boundaryField()[patchI][faceI] * (virgin.hs(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) - virgin.hs(p_.boundaryField()[patchI][faceI], 298)) -
                    rho_c_.boundaryField()[patchI][faceI] * (charred.hs(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) - charred.hs(p_.boundaryField()[patchI][faceI], 298))
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

void Foam::PorousMaterialPropertiesModel::updateHc(word name)
{
  LinearInterpolationMaterialPropertiesObject& charred = charObject();
  volScalarField& h_c_ = scalarFields_[name];
  forAll(h_c_, cellI) {
    h_c_[cellI] = charred.h(p_[cellI], T_[cellI]);
  }


  forAll(mesh_.boundaryMesh(), patchI) {
    const fvPatchScalarField& Tp = T_.boundaryField()[patchI];
    const fvPatchScalarField& pp = p_.boundaryField()[patchI];
    fvPatchScalarField& h_cp = h_c_.boundaryFieldRef()[patchI];

    forAll(h_cp, faceI) {
      h_cp[faceI] = charred.h(pp[faceI], Tp[faceI]);
    }
  }
}

void Foam::PorousMaterialPropertiesModel::updateEmissivityAbsorptivity(word name)
{
  LinearInterpolationMaterialPropertiesObject& virgin = virginObject();
  LinearInterpolationMaterialPropertiesObject& charred = charObject();

  volScalarField& ea_ = scalarFields_[name]; // emissivity or absorptivity
  forAll(mesh_.boundaryMesh(), patchI) {
    const fvPatchScalarField& Tp = T_.boundaryField()[patchI];
    const fvPatchScalarField& pp = p_.boundaryField()[patchI];
    const fvPatchScalarField& rho_sp = rho_s_.boundaryField()[patchI];
    const fvPatchScalarField& taup = tau_.boundaryField()[patchI];
    fvPatchScalarField& eap = ea_.boundaryFieldRef()[patchI];

    if (name != "emissivity" && name != "absorptivity") {
      FatalError << "void Foam::PorousMaterialPropertiesModel::updateEmissivityAbsorptivity(word name): name (" << name
                 << ") must be \"emissivity\" or \"absorptivity\"" << exit(FatalError);
    }
    forAll(eap, faceI) {
      if (name == "emissivity") {
        eap[faceI] = virgin.eps(pp[faceI], Tp[faceI]) * taup[faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_sp[faceI], rho_c_.boundaryField()[patchI][faceI])
                     + charred.eps(pp[faceI], Tp[faceI]) * (1 - taup[faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_sp[faceI], rho_c_.boundaryField()[patchI][faceI]));
      }
      if (name == "absorptivity") {
        eap[faceI] = virgin.alpha(pp[faceI], Tp[faceI]) * taup[faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_sp[faceI], rho_c_.boundaryField()[patchI][faceI])
                     + charred.alpha(pp[faceI], Tp[faceI]) * (1 - taup[faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_sp[faceI], rho_c_.boundaryField()[patchI][faceI]));
      }
    }
  }
}

LinearInterpolationMaterialPropertiesObject Foam::PorousMaterialPropertiesModel::readCharredVirginDict(word charredVirgin)
{
  // read and store the material property table into the RAM for faster access
  Info << getTabLevel() << "Reading " << charredVirgin << " material properties" << nl;

  fileName charredVirginFileName = materialDirectory_+"/"+charredVirgin;
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
    Info << getTabLevel() << "The " << charredVirgin << " material table has " << rawTableV << " lines."<< nl;
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
//  Info << charredVirginThermoTable[0][1] << " " << charredVirginThermoTable[0][2] << " " << charredVirginThermoTable[0][9] << endl;
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

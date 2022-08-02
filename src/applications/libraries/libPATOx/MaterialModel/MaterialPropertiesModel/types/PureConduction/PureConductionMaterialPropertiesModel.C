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
    along with OpenFOAM.  If PureConduction_t, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "PureConductionMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PureConductionMaterialPropertiesModel::PureConductionMaterialPropertiesModel
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
virginOrChar(simpleMaterialPropertiesModel::materialDict_.subDict("MaterialProperties").template lookupOrDefault<word>("virginOrChar","char")),
nSolidPhases_(simpleMaterialPropertiesModel::nSolidPhases_=readScalar(constantPropertiesDictionary.lookup("nSolidPhases"))),
solidEpsI_(simpleMaterialPropertiesModel::solidEpsI_),
solidRhoI_(simpleMaterialPropertiesModel::solidRhoI_),
kiCoef_(constantPropertiesDictionary.lookupOrDefault<dimensionedScalar>("kiCoef", dimensionedScalar("0",dimensionSet(0,0,0,0,0),0))),
kjCoef_(constantPropertiesDictionary.lookupOrDefault<dimensionedScalar>("kjCoef", dimensionedScalar("0",dimensionSet(0,0,0,0,0),0))),
kkCoef_(constantPropertiesDictionary.lookupOrDefault<dimensionedScalar>("kkCoef", dimensionedScalar("0",dimensionSet(0,0,0,0,0),0))),
tP_(constantPropertiesDictionary.lookupOrDefault<dimensionedTensor>("tP",tensor(1,0,0,0,1,0,0,0,1))),
solidRho_(simpleMaterialPropertiesModel::solidRho_),
solidEps_(simpleMaterialPropertiesModel::solidEps_),
PyrolysisModel_(meshLookupOrConstructModel<simplePyrolysisModel>(mesh,dictName,"Pyrolysis")),
tau_(PyrolysisModel_.tau()),
rho_s_(simpleMaterialPropertiesModel::rho_s_),
rho_v_(PyrolysisModel_.rho_v()),
rho_c_(PyrolysisModel_.rho_c()),
T_(meshLookupOrConstructScalar(mesh, "Ta")),
p_(meshLookupOrConstructScalar(mesh, "p")),
cp_(simpleMaterialPropertiesModel::cp_),
k_(simpleMaterialPropertiesModel::k_),
absorptivity_(simpleMaterialPropertiesModel::absorptivity_),
emissivity_(simpleMaterialPropertiesModel::emissivity_),
h_c_(simpleMaterialPropertiesModel::h_c_),
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

)

{
  if (!isFile(constantPropertiesDictionary.path()+"/constantProperties")) {
    FatalErrorInFunction << "Unknown materialPropertiesDirectory " << constantPropertiesDictionary.path()+"/constantProperties"
                         << exit(FatalError);

  }

  if (virginOrChar != "char" && virginOrChar != "virgin") {
    FatalErrorInFunction
        <<  "Problem in " << dictName << "MaterialProperties :" << nl << "\tvirginOrChar " << virginOrChar << ";" << nl << "\t"
        << "virginOrChar has to be \"virgin\" or \"char\" "
        << exit(FatalError);

  }

  readEpsI();
  readRhoI();

  virginObject.reset(new LinearInterpolationMaterialPropertiesObject(readCharredVirginDict("virgin")));
  charObject.reset(new LinearInterpolationMaterialPropertiesObject(readCharredVirginDict("char")));

  rho_s_*=0;
  forAll(solidRho_, phaseI) {
    rho_s_+= solidEps_[phaseI]*solidRho_[phaseI];
  }

  word typeNamePyrolysisModel = "noPyrolysisModel<simplePyrolysisModel>";
  if (this->materialDict_.isDict("Pyrolysis")) {
    typeNamePyrolysisModel = this->materialDict_.subDict("Pyrolysis").template lookupOrDefault<word>("PyrolysisType","noPyrolysisModel<simplePyrolysisModel>");
  }

  typeNamePyrolysisModel.replaceAll("PyrolysisModel<simplePyrolysisModel>","");

  if (typeNamePyrolysisModel=="no") {
    rho_v_ *= 0;
    forAll(solidRho_, phaseI) {
      rho_v_ += solidEpsI_[phaseI] * solidRhoI_[phaseI];
    }
    rho_c_=rho_v_;

    if (virginOrChar == "char") {
      forAll(tau_, cellI) {
        tau_[cellI]=0;
      }
      forAll(mesh.boundaryMesh(), patchI) {
        forAll(tau_.boundaryFieldRef()[patchI], faceI) {
          tau_.boundaryFieldRef()[patchI][faceI] = 0;
        }
      }
    }
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PureConductionMaterialPropertiesModel::~PureConductionMaterialPropertiesModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PureConductionMaterialPropertiesModel::update()
{
  // update solid specific heat
  updateCp();

  // update solid thermal conductivity
  updatek();

  // update emissivity and absorptivity
  updateEmissivityAbsorptivity();

  // update charred enthalpy
  updateHc();
}

void Foam::PureConductionMaterialPropertiesModel::updateHc()
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

void Foam::PureConductionMaterialPropertiesModel::updateEmissivityAbsorptivity()
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
      emissivityp[faceI] = virgin.eps(pp[faceI], Tp[faceI]) * taup[faceI] * rho_v_.boundaryField()[patchI][faceI] / max(rho_sp[faceI], rho_c_.boundaryField()[patchI][faceI] )
                           + charred.eps(pp[faceI], Tp[faceI]) * (1 - taup[faceI] * rho_v_.boundaryField()[patchI][faceI]  / max(rho_sp[faceI], rho_c_.boundaryField()[patchI][faceI] ));
    }

    forAll(absorptivityp, faceI) {
      absorptivityp[faceI] = virgin.alpha(pp[faceI], Tp[faceI]) * taup[faceI] * rho_v_.boundaryField()[patchI][faceI]  / max(rho_sp[faceI], rho_c_.boundaryField()[patchI][faceI] )
                             + charred.alpha(pp[faceI], Tp[faceI]) * (1 - taup[faceI] * rho_v_.boundaryField()[patchI][faceI]  / max(rho_sp[faceI], rho_c_.boundaryField()[patchI][faceI] ));
    }

  }
}

void Foam::PureConductionMaterialPropertiesModel::updatek()
{
  LinearInterpolationMaterialPropertiesObject& virgin = virginObject();
  LinearInterpolationMaterialPropertiesObject& charred = charObject();
  forAll(k_, i) {
    k_[i].xx() =
        kiCoef_.value() *
        (
            virgin.ki(p_[i], T_[i]) *
            tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i])
            +
            charred.ki(p_[i], T_[i]) *
            (1. - tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i]))
        );
    k_[i].yy() =
        kjCoef_.value() *
        (
            virgin.kj(p_[i], T_[i]) *
            tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i])
            +
            charred.kj(p_[i], T_[i]) *
            (1. - tau_[i] * rho_v_[i] / max(rho_s_[i], rho_c_[i]))
        );
    k_[i].zz() =
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
    forAll(k_.boundaryFieldRef()[patchI], faceI) {
      k_.boundaryFieldRef()[patchI][faceI].xx() =
          kiCoef_.value() *
          (
              virgin.ki(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
              tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI]  / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI] )
              +
              charred.ki(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
              (1. - tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI]  / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI] ))
          );
      k_.boundaryFieldRef()[patchI][faceI].yy() =
          kjCoef_.value() *
          (
              virgin.kj(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
              tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI]  / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI] )
              +
              charred.kj(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
              (1. - tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI]  / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI] ))
          );
      k_.boundaryFieldRef()[patchI][faceI].zz() =
          kkCoef_.value() *
          (
              virgin.kk(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
              tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI]  / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI] )
              +
              charred.kk(p_.boundaryField()[patchI][faceI], T_.boundaryField()[patchI][faceI]) *
              (1. - tau_.boundaryField()[patchI][faceI] * rho_v_.boundaryField()[patchI][faceI]  / max(rho_s_.boundaryField()[patchI][faceI], rho_c_.boundaryField()[patchI][faceI] ))
          );
    }
  }


  k_ = tP_.value() & k_ & tP_.value().T(); // passage from basis ijk (main directions of the conductivity tensor) to basis xyz (mesh)
  k_abl_sym_ = symm(k_);
}

void Foam::PureConductionMaterialPropertiesModel::updateCp()
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

void Foam::PureConductionMaterialPropertiesModel::readEpsI()
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

  for (int i = 0; i < nSolidPhases_; i++) { // the solid phases are defined from value 1 (0 is attributed to the gas).
    dimensionedScalar solidEpsI_value = constantPropertiesDictionary.lookup("epsI[" + std::to_string(i+1) + "]");
    solidEpsI_[i] = solidEpsI_value;
    solidEps_[i] = solidEpsI_value;
    Info << "epsI[" << i+1 << "]=" << solidEpsI_value.value() <<  nl;
  }
}

void Foam::PureConductionMaterialPropertiesModel::readRhoI()
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
    solidRhoI_[i] = solidRhoI_value;
    solidRho_[i] = solidRhoI_value;
    Info << "rhoI["<< i+1 << "]=" << solidRhoI_value.value() <<  nl;
  }
}

Foam::LinearInterpolationMaterialPropertiesObject Foam::PureConductionMaterialPropertiesModel::readCharredVirginDict(word charredVirgin)
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
  charredVirginThermoTable[0][9]=0;
  for(int x=1; x<rawTableV; x++) {
    charredVirginThermoTable[x][9] += charredVirginThermoTable[x-1][9]
                                      + (charredVirginThermoTable[x][2] + charredVirginThermoTable[x-1][2])/2
                                      * (charredVirginThermoTable[x][1] - charredVirginThermoTable[x-1][1]);
  }

//        // "virgin" is an object of the class materialProperties, which handles the material properties interpolations when called
  LinearInterpolationMaterialPropertiesObject objectCharredVirgin(charredVirginThermoTable);
  return objectCharredVirgin;
//
  //Info << "The " << charredVirgin << " material properties have been read." << nl;
  // END OF THE IO TEST

}


// ************************************************************************* //

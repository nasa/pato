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
    along with OpenFOAM.  If BprimeCoatingMixturet, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "BprimeCoatingMixtureBoundaryConditions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BprimeCoatingMixtureBoundaryConditions::BprimeCoatingMixtureBoundaryConditions
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& regionName,
    const label& currentPatchID,
    const dictionary dict
)
  :
mesh_(mesh),
phaseName_(phaseName),
regionName_(regionName),
energyModel_(meshLookupOrConstructModel<simpleEnergyModel>(mesh_,regionName_,simpleEnergyModel::modelName)),
initialized_(false),
currentPatchID_(currentPatchID),
dict_(dict),
materialDict_(energyModel_.materialDict()),
dynamicMesh_(isA<dynamicFvMesh>(mesh)),
mixtureName(dict_.lookup("mixtureName")),
environmentDirectory(fileName(dict_.lookup("environmentDirectory")).expand()),
environmentDictionary
(
    IOobject
    (
        "environmentComposition",
        environmentDirectory,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
minBc(dict_.lookupOrDefault<scalar>("minBc",1e-5)),
maxBc(dict_.lookupOrDefault<scalar>("maxBc",100)),
sigmaSB(::constant::physicoChemical::sigma),
initialPosition_(energyModel_.createVolField<vector>("initialPosition",dimensionedVector("0", dimLength, vector(0.0,0.0,0.0)))),
initialRecession_(energyModel_.createVolField<scalar>("initialRecession",dimensionedScalar("0", dimLength, scalar(0.0)))),
debug_(materialDict_.lookupOrDefault<Switch>("debug","no")),
neededFields_(neededFields()),
mode_surf(dict_.lookupOrDefault<int>("modeSurfComp",0)),
ne_mix(0),
ns_mix(0),
recessionFlag(dict_.lookupOrDefault<Switch>("recessionFlag","yes"))
{
  // Create new fields in Energy Model
  scalarFields_.insert("mDotCw",energyModel_.createVolField<scalar>("mDotCw",dimensionedScalar("0", dimMass/pow(dimLength,2)/dimTime, scalar(0.0))));
  scalarFields_.insert("BprimeCw",energyModel_.createVolField<scalar>("BprimeCw",dimensionedScalar("0", dimless, scalar(0.0))));
  scalarFields_.insert("BprimeGw",energyModel_.createVolField<scalar>("BprimeGw",dimensionedScalar("0", dimless, scalar(0.0))));
  scalarFields_.insert("qRadEmission",energyModel_.createVolField<scalar>("qRadEmission",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("qRadAbsorption",energyModel_.createVolField<scalar>("qRadAbsorption",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("qAdvPyro",energyModel_.createVolField<scalar>("qAdvPyro",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("qAdvChar",energyModel_.createVolField<scalar>("qAdvChar",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("qConv",energyModel_.createVolField<scalar>("qConv",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("qCond",energyModel_.createVolField<scalar>("qCond",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("recession",energyModel_.createVolField<scalar>("recession", dimensionedScalar("0", dimLength, scalar(0.0))));
  scalarFields_.insert("rhoeUeCH",energyModel_.createVolField<scalar>("rhoeUeCH", dimensionedScalar("0", dimMass/ pow3(dimLength)/ dimTime, scalar(0.0))));
  scalarFields_.insert("blowingCorrection",energyModel_.createVolField<scalar>("blowingCorrection",dimensionedScalar("1", dimless, scalar(1.0))));
  scalarFields_.insert("h_w",energyModel_.createVolField<scalar>("h_w",dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0.0))));
  scalarFields_.insert("h_r",energyModel_.createVolField<scalar>("h_r",dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0.0))));
  scalarFields_.insert("Tbackground",energyModel_.createVolField<scalar>("Tbackground", dimensionedScalar("0", dimTemperature, scalar(0.0))));
  scalarFields_.insert("lambda",energyModel_.createVolField<scalar>("lambda", dimensionedScalar("0", dimless, scalar(0.0))));
  scalarFields_.insert("qRad",energyModel_.createVolField<scalar>("qRad", dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("chemistryOn",energyModel_.createVolField<scalar>("chemistryOn", dimensionedScalar("0", dimless, scalar(0.0))));
  scalarFields_.insert("h_c",energyModel_.createVolField<scalar>("h_c",dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0)));
  scalarFields_.insert("emissivity",energyModel_.createVolField<scalar>("emissivity",dimensionedScalar("0", dimless, scalar(0))));
  scalarFields_.insert("absorptivity",energyModel_.createVolField<scalar>("absorptivity",dimensionedScalar("0", dimless, scalar(0))));

  // Initialize initialRecession_
  volScalarField& recession_=scalarFields_["recession"];
  forAll(initialRecession_, cellI) {
    initialRecession_[cellI]=recession_[cellI];
  }
  forAll(initialRecession_.boundaryField(), patchI) {
    forAll(initialRecession_.boundaryField()[patchI], faceI) {
      initialRecession_.boundaryFieldRef()[patchI][faceI]=recession_.boundaryField()[patchI][faceI];
    }
  }

  if(debug_) {
    Info << "--- environmentDirectory=" << environmentDirectory << " --- Foam::BprimeCoatingMixtureBoundaryConditions::BprimeCoatingMixtureBoundaryConditions" << endl;
    Info << "--- dict_=" << dict_ << " --- Foam::BprimeCoatingMixtureBoundaryConditions::BprimeCoatingMixtureBoundaryConditions" << endl;
  }

  initMixture();
  if(dynamicMesh_) {

    cellMotionU_ptr = &const_cast<volVectorField&>(mesh_.objectRegistry::lookupObject<volVectorField>("cellMotionU"));

    forAll(initialPosition_.boundaryField()[currentPatchID_], faceI) {
      initialPosition_.boundaryFieldRef()[currentPatchID_][faceI] =  mesh_.Cf().boundaryField()[currentPatchID_][faceI];
    }
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::BprimeCoatingMixtureBoundaryConditions::~BprimeCoatingMixtureBoundaryConditions()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::BprimeCoatingMixtureBoundaryConditions::initialize()
{
  word bold_on="\e[1m";
  word bold_off="\e[0m";
  Info << simpleModel::getTabLevel("| ") << bold_on << "Initialize Temperature BC (Bprime)" <<  bold_off << endl;
  simpleModel::tabLevel_++;
  // References to models
  const simpleMassModel& massModel=meshLookupOrConstructModel<simpleMassModel>(mesh_,regionName_,simpleMassModel::modelName);
  const simpleGasPropertiesModel& gasPropertiesModel=meshLookupOrConstructModel<simpleGasPropertiesModel>(mesh_,regionName_,simpleGasPropertiesModel::modelName);
  // References to fields from Mass Model
  scalarFields_.insert("p",massModel.refVolField<scalar>("p"));
  vectorFields_.insert("mDotG",massModel.refVolField<vector>("mDotG"));
  scalarFields_.insert("mDotGw",massModel.refVolField<scalar>("mDotGw"));
  // References to fields from Gas Properties Model
  scalarFields_.insert("h_g",gasPropertiesModel.refVolField<scalar>("h_g"));
  // References to fields from Energy Model
  scalarFields_.insert("Ta",energyModel_.refVolField<scalar>("Ta"));
  scalarFields_.insert("rho_s",energyModel_.refVolField<scalar>("rho_s"));
  tensorFields_.insert("k",energyModel_.refVolField<tensor>("k"));
  // Boundary Mapping
  boundaryMapping_=simpleBoundaryMappingModel::New(mesh_, neededFields_,dict_);
  boundaryMapping_ptr=&boundaryMapping_();
  // Initialized flag
  initialized_=true;
  simpleModel::tabLevel_--;
}

List<Tuple2<string,double>> Foam::BprimeCoatingMixtureBoundaryConditions::computeSiCSiO2SurfComp(const Mutation::Mixture& mix, const double& zC)
{
  Eigen::ArrayXd mm (mix.nSpecies());
  for (int i = 0; i < mix.nSpecies(); i++) {
    mm[i] = mix.speciesMw(i);
  }

  // Arranging Species
  const int ns_e = 3;
  const int eSi = 0;
  string s_Si = "Si";
  const int eO = 1;
  string s_O = "O";
  const int eC = 2;
  string s_C = "C";

  // Computing molar masses
  const double mmSi = mm(mix.speciesIndex(s_Si));
  const double mmO  = mm(mix.speciesIndex(s_O));
  const double mmC  = mm(mix.speciesIndex(s_C));
  const double mmSiC = mmSi + mmC;
  const double mmSiO2 = mmSi + 2.*mmO;

  // Mass Fraction of Surface Species
  double zSiC = zC + mmSi/mmC*zC;
  double zSiO2 = 1.-zSiC;

  if (zSiO2 < 0.) {
    FatalErrorInFunction << "Such a composition is not possible. yC = " << zC << exit(FatalError);
  }

  List<Tuple2<string,double>> zs_e(ns_e);
  zs_e[eSi].first() = s_Si;
  zs_e[eSi].second() = mmSi/mmSiC*zSiC + mmSi/mmSiO2*zSiO2;
  zs_e[eO].first() = s_O;
  zs_e[eO].second() = 2.*mmO/mmSiO2*zSiO2;
  zs_e[eC].first() = s_C;
  zs_e[eC].second() = zSiC*mmC/mmSiC;

  return zs_e;
}

void Foam::BprimeCoatingMixtureBoundaryConditions::initMixture()
{
  mixture_ptr= new Mutation::Mixture(mixtureName);
  Mutation::Mixture& mixture = *mixture_ptr ;

  ne_mix= mixture.nElements();
  ns_mix= mixture.nSpecies();
  Yke_ref = new double[ne_mix];
  Ykg_ref= new double[ne_mix];

  // Initialization of Mutation++ for surface equilibrium chemistry (B')
  // Default conditions (input in mole fractions -> converted in mass fractions)
  Info << "Reading boundary-layer edge elemental composition (mass fraction)" << nl;
  for (int i = 0; i < ne_mix; i++) {
    Yke_ref[i] =
        environmentDictionary.lookupOrDefault<scalar>
        (
            "Yke["+  mixture.elementName(i) + "]",
            0.0
        );

  }

  mixture.convert<Mutation::Thermodynamics::XE_TO_YE>(Yke_ref, Yke_ref);
  for (int i = 0; i < ne_mix; i++) {
    Info << "Yke[" <<  mixture.elementName(i) << "] = " << Yke_ref[i] << nl;
  }

  Info << "Reading initial elemental composition of the gas inside the material for B' solver (mass fraction)" << nl;

  dictionary Ykg_dict = dict_;
  if (materialDict_.isDict("MaterialProperties")) {
    //- material properties directory
    const fileName materialPropertiesDirectory
    (
        fileName(materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand()
    );

    //- constantProperties dictionary
    IOdictionary constantPropertiesDictionary
    (
        IOobject
        (
            "constantProperties",
            materialPropertiesDirectory,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    Ykg_dict=constantPropertiesDictionary;
  }
  bool atLeastOneElement = false;
  for (int i = 0; i < ne_mix; i++) {
    Ykg_ref[i] =
        Ykg_dict.lookupOrDefault<scalar>
        (
            "Zx[" +  mixture.elementName(i) + "]",
            0.0
        );
    if (Ykg_ref[i]!=0) {
      atLeastOneElement = true;
    }
  }

  if (!atLeastOneElement) {
    FatalErrorInFunction << "Zx[...] not found in " <<  Ykg_dict << exit(FatalError);
  }

  mixture.convert<Mutation::Thermodynamics::XE_TO_YE>(Ykg_ref, Ykg_ref);

  for (int i = 0; i < ne_mix; i++) {
    Info << "Ykg[" <<  mixture.elementName(i) << "] = " << Ykg_ref[i] << nl;
  }

  for (int elementI=0; elementI<ne_mix; elementI++) {
    word nameYkg_ = "Ykg["+mixture.elementName(elementI)+"]";
    Ykg_names.append(nameYkg_);
    volScalarField * Ykg_ = &energyModel_.createVolField<scalar>(nameYkg_.c_str(),dimensionedScalar("0", dimless,  Ykg_ref[elementI]),wordList(mesh_.boundaryMesh().size(),"fixedValue"));
    Ykg_->store(Ykg_);
  }
  for (int speciesI=0; speciesI<ns_mix; speciesI++) {
    word nameXw_ = "Xw["+mixture.speciesName(speciesI)+"]";
    Xw_names.append(nameXw_);
    volScalarField * Xw_ = &energyModel_.createVolField<scalar>(nameXw_.c_str(),dimensionedScalar("0", dimless, scalar(0)),wordList(mesh_.boundaryMesh().size(),"fixedValue"));
    Xw_->store(Xw_);
  }

  // Surface Composition based on different modes.
  if (mode_surf == 0) {
    List<Tuple2<string,double>> elemOnSurface_dict(dict_.lookup("elemOnSurface"));
    elemOnSurface.resize(elemOnSurface_dict.size());
    forAll(elemOnSurface_dict, i) {
      elemOnSurface[i] = elemOnSurface_dict[i]; // list of the elements on the surface + mass fractions
    }
  } else if (mode_surf == 1) {
    const double zC = readScalar(dict_.lookup("massCSurface"));
    elemOnSurface = computeSiCSiO2SurfComp(mixture, zC);
  } else {
    FatalErrorInFunction << "Unrecognised mode for imposing surface composition: " << mode_surf << exit(FatalError);
  }

  // Verify elemOnSurface
  for (int eI=0; eI < elemOnSurface.size(); eI++) {
    bool isElem = false;
    for (int i = 0; i < ne_mix; i++) {
      if (elemOnSurface[eI].first()==mixture.elementName(i)) {
        isElem = true;
        break;
      }
    }
    if (!isElem) {
      FatalErrorInFunction << elemOnSurface[eI].first() << " element not found in the mixture." << exit(FatalError);
    }
  }

  Info << "Surface elemental composition: " << endl;
  for (int i=0; i < elemOnSurface.size(); i++) {
    Info << elemOnSurface[i].first() << ": " << elemOnSurface[i].second() << endl;
  }
}

void Foam::BprimeCoatingMixtureBoundaryConditions::update()
{
  if (!initialized_) {
    initialize();
  }

  if (debug_) {
    Info << "--- update BoundaryMapping --- Foam::BprimeCoatingMixtureBoundaryConditions::update()" << endl;
  }
  updateBoundaryMapping();

  if (debug_) {
    Info << "--- update BprimeCoatingMixture --- Foam::BprimeCoatingMixtureBoundaryConditions::update()" << endl;
  }
  updateBprimeBC();
  if (debug_) {
    Info << "--- update temperature --- Foam::BprimeCoatingMixtureBoundaryConditions::update()" << endl;
  }
  updateTemperatureBC();

  if(dynamicMesh_) {
    if (debug_) {
      Info << "--- update motion --- Foam::BprimeCoatingMixtureBoundaryConditions::update()" << endl;
    }
    updateMotionBC();

    if(this->debug_) {
      Info << "--- end --- Foam::BprimeCoatingMixtureBoundaryConditions::update()" << endl;
    }
  }
}

void Foam::BprimeCoatingMixtureBoundaryConditions::updateBoundaryMapping()
{
  forAll(neededFields_, fieldI) {
    // Update all the needed fields with boundaryMappingFvPatchScalarField
    if (isA<boundaryMappingFvPatchScalarField>(const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(neededFields_[fieldI])).boundaryField()[currentPatchID_])) {
      const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(neededFields_[fieldI])).correctBoundaryConditions();
    } else {
      if (foundInList(neededFields_[fieldI], boundaryMapping_ptr->mappingFieldsName())) { // Update all the fields in mappingFields
        boundaryMapping_ptr->update(mesh_.time().value(),currentPatchID_,neededFields_[fieldI]);
      }
    }
  }
}

void Foam::BprimeCoatingMixtureBoundaryConditions::updateMotionBC()
{
  // Reference to fields
  volScalarField& rhoeUeCH=scalarFields_["rhoeUeCH"];
  volScalarField& blowingCorrection=scalarFields_["blowingCorrection"];
  volScalarField& BprimeCw_=scalarFields_["BprimeCw"];
  volScalarField& rho_s_=scalarFields_["rho_s"];
  volScalarField& recession_=scalarFields_["recession"];

  if(mesh_.objectRegistry::foundObject<volVectorField>("cellMotionU")) {
    volVectorField& cellMotionU_ = *cellMotionU_ptr;
    forAll(cellMotionU_.boundaryFieldRef()[currentPatchID_], faceI) {
      // not used yet
      scalar failureFraction=0;
      scalar erosionRate=0;

      // Already updated
      const scalar& rhoeUeCH_BF = rhoeUeCH.boundaryField()[currentPatchID_][faceI];
      const vector nf =
          - mesh_.Sf().boundaryField()[currentPatchID_][faceI]
          / mesh_.magSf().boundaryField()[currentPatchID_][faceI];
      const scalar& blowingCorrection_BF = blowingCorrection.boundaryField()[currentPatchID_][faceI];
      const scalar& BprimeCw_BF = BprimeCw_.boundaryField()[currentPatchID_][faceI];
      const scalar& rho_s_bf = rho_s_.boundaryField()[currentPatchID_][faceI];

      // fluxFactorThreshold (2D axi boundary mapping): Mesh motion direction is projected on the direction of pointMotionDirection when the value of the fluxFactor is higher than the set Threshold, otherwise motion is applied directly
      // along the face normal. This is useful to preserve the topology for large and strongly distorted deformations.
      vector normal_ = nf;
      if (isA<boundaryMappingFvPatchScalarField>(rhoeUeCH.boundaryField()[currentPatchID_])) {
        boundaryMappingFvPatchScalarField& rhoeUeCH_bm=refCast<boundaryMappingFvPatchScalarField>(rhoeUeCH.boundaryFieldRef()[currentPatchID_]);
        autoPtr<simpleBoundaryMappingModel>& boundaryMapping_rhoeUeCH = rhoeUeCH_bm.boundaryMapping_;
        if (boundaryMapping_rhoeUeCH->fluxFactor_ptr) {
          const scalar& fluxThreshold_ = boundaryMapping_rhoeUeCH->fluxFactor_ptr->fluxFactorThreshold().value();
          const scalar& fluxMap_ = boundaryMapping_rhoeUeCH->fluxFactor_ptr->fluxMap().boundaryField()[currentPatchID_][faceI];
          const vector& pointMotionDirection_ = boundaryMapping_rhoeUeCH->fluxFactor_ptr->pointMotionDirection();
          const Switch& moreThanTresholdFlag_ = boundaryMapping_rhoeUeCH->fluxFactor_ptr->moreThanThresholdFlag();
          if (moreThanTresholdFlag_) {
            if (fluxMap_ > fluxThreshold_) {
              normal_ = (1 / (nf & pointMotionDirection_)) * pointMotionDirection_;
            }
          } else {
            if (fluxMap_ < fluxThreshold_) {
              normal_ = (1 / (nf & pointMotionDirection_)) * pointMotionDirection_;
            }
          }
        }
      }

      // Will be updated
      vector&  cellMotionU_BF = cellMotionU_.boundaryFieldRef()[currentPatchID_][faceI];

      cellMotionU_BF=
          (BprimeCw_BF* rhoeUeCH_BF * blowingCorrection_BF * (1 + failureFraction)
           + erosionRate)
          * 1 / rho_s_bf * normal_;
      if (!recessionFlag) {
        cellMotionU_BF = vector(0,0,0);
      }
    }
  }

  forAll(recession_.boundaryField()[currentPatchID_], faceI) {
    recession_.boundaryFieldRef()[currentPatchID_][faceI] = initialRecession_.boundaryField()[currentPatchID_][faceI] + mag(initialPosition_.boundaryField()[currentPatchID_][faceI] -  mesh_.Cf().boundaryField()[currentPatchID_][faceI]);
  }
}

void Foam::BprimeCoatingMixtureBoundaryConditions::updateTemperatureBC()
{
  // Reference to fields
  volScalarField& T_=scalarFields_["Ta"];
  volScalarField& h_r=scalarFields_["h_r"];
  volScalarField& Tbackground=scalarFields_["Tbackground"];
  volScalarField& qRad=scalarFields_["qRad"];
  volScalarField& chemistryOn=scalarFields_["chemistryOn"];
  volScalarField& rhoeUeCH=scalarFields_["rhoeUeCH"];
  volScalarField& mDotGw_=scalarFields_["mDotGw"];
  volScalarField& h_g_=scalarFields_["h_g"];
  volScalarField& h_c_=scalarFields_["h_c"];
  volScalarField& emissivity_=scalarFields_["emissivity"];
  volScalarField& absorptivity_=scalarFields_["absorptivity"];
  volTensorField& k_=tensorFields_["k"];
  volScalarField& blowingCorrection=scalarFields_["blowingCorrection"];
  volScalarField& BprimeCw_=scalarFields_["BprimeCw"];
  volScalarField& h_w=scalarFields_["h_w"];
  volScalarField& qRadEmission_=scalarFields_["qRadEmission"];
  volScalarField& qRadAbsorption_=scalarFields_["qRadAbsorption"];
  volScalarField& qAdvPyro_=scalarFields_["qAdvPyro"];
  volScalarField& qAdvChar_=scalarFields_["qAdvChar"];
  volScalarField& qConv_=scalarFields_["qConv"];
  volScalarField& qCond_=scalarFields_["qCond"];

  forAll(T_.boundaryFieldRef()[currentPatchID_], faceI) {
    // not used yet
    scalar wallSpeciesDiffusion = 1;

    // Updated by BoundaryMapping
    scalar& h_r_BF = h_r.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& Tbackground_BF = Tbackground.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qRad_BF = qRad.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& chemistryOn_BF = chemistryOn.boundaryFieldRef()[currentPatchID_][faceI];

    // Already updated
    const vector nf =
        - mesh_.Sf().boundaryField()[currentPatchID_][faceI]
        / mesh_.magSf().boundaryField()[currentPatchID_][faceI];
    const scalar& rhoeUeCH_BF = rhoeUeCH.boundaryFieldRef()[currentPatchID_][faceI];
    const scalar invDx_BF = mesh_.deltaCoeffs().boundaryField()[currentPatchID_][faceI];
    const tmp<scalarField> Tint_tmp = T_.boundaryField()[currentPatchID_].patchInternalField();
    const scalarField& Tint_ = Tint_tmp();
    const scalar Tint_BF = Tint_[faceI];
    const scalar mDotGw_BF = mDotGw_.boundaryField()[currentPatchID_][faceI];
    const scalar h_g_BF = h_g_.boundaryField()[currentPatchID_][faceI];
    const scalar h_c_BF = h_c_.boundaryField()[currentPatchID_][faceI];
    const scalar emissivity_BF = emissivity_.boundaryField()[currentPatchID_][faceI];
    const scalar absorptivity_BF = absorptivity_.boundaryField()[currentPatchID_][faceI];
    const tensor k_BF = k_.boundaryField()[currentPatchID_][faceI];
    const scalar kProj_BF = nf & k_BF & nf; // Projection of k on the surface normal
    const scalar& blowingCorrection_BF = blowingCorrection.boundaryFieldRef()[currentPatchID_][faceI];
    const scalar& BprimeCw_BF = BprimeCw_.boundaryFieldRef()[currentPatchID_][faceI];
    const scalar& h_w_BF = h_w.boundaryFieldRef()[currentPatchID_][faceI];

    // Will be updated
    scalar& qRadEmission_BF = qRadEmission_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qRadAbsorption_BF = qRadAbsorption_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qAdvPyro_BF = qAdvPyro_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qAdvChar_BF = qAdvChar_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qConv_BF = qConv_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qCond_BF = qCond_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& T_BF = T_.boundaryFieldRef()[currentPatchID_][faceI];

    // Surface energy balance: heat fluxes and temperature
    qRadEmission_BF = - emissivity_BF * sigmaSB.value() * (pow4(T_BF) - pow4(Tbackground_BF));
    qRadAbsorption_BF = absorptivity_BF * qRad_BF;

    if(chemistryOn_BF == 1) {
      qAdvPyro_BF = mDotGw_BF * (h_g_BF - h_w_BF);
      qAdvChar_BF =  BprimeCw_BF * rhoeUeCH_BF * blowingCorrection_BF * (h_c_BF - h_w_BF);
      qConv_BF = rhoeUeCH_BF * blowingCorrection_BF * (h_r_BF - h_w_BF * wallSpeciesDiffusion) ;
      T_BF =
          Tint_BF
          + (
              1. / (kProj_BF * invDx_BF) *
              (
                  qConv_BF + qAdvPyro_BF + qAdvChar_BF + qRadEmission_BF + qRadAbsorption_BF
              )
          );
    } else {
      qAdvChar_BF = 0;
      qConv_BF = 0;
      qAdvPyro_BF = 0;
      T_BF =
          Tint_BF
          + (
              1. / (kProj_BF * invDx_BF) *
              (
                  qRadEmission_BF + qRadAbsorption_BF
              )
          );
    }


    qCond_BF = (kProj_BF * invDx_BF)* (T_BF - Tint_BF);
  }



}

void Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()
{
  // Reference to fields
  volScalarField& BprimeCw_=scalarFields_["BprimeCw"];
  volScalarField& mDotCw_=scalarFields_["mDotCw"];
  volScalarField& BprimeGw_=scalarFields_["BprimeGw"];
  volScalarField& mDotGw_=scalarFields_["mDotGw"];
  volScalarField& rhoeUeCH=scalarFields_["rhoeUeCH"];
  volScalarField& lambda=scalarFields_["lambda"];
  volScalarField& chemistryOn=scalarFields_["chemistryOn"];
  volScalarField& p_=scalarFields_["p"];
  volScalarField& T_=scalarFields_["Ta"];
  volScalarField& blowingCorrection=scalarFields_["blowingCorrection"];
  volScalarField& h_w=scalarFields_["h_w"];

  if (debug_) {
    Info << "--- currentPatchID_=" << currentPatchID_ << " --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()"  << endl;
    Info << "--- update Boundary Mapping for patch " << mesh_.boundaryMesh()[currentPatchID_].name() << " --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()"  << endl;
  }

  forAll(BprimeCw_.boundaryField()[currentPatchID_], faceI) {

    // Species mole fraction at the wall
    Eigen::VectorXd p_Xw_(Eigen::VectorXd::Zero(ns_mix));
    // Elements mass fraction inside the material
    Eigen::VectorXd p_Ykg_(Eigen::VectorXd::Zero(ne_mix));
    // Elements mass fraction in the environment
    Eigen::VectorXd p_Yke_(Eigen::VectorXd::Zero(ne_mix));

    if (debug_ && faceI==0) {
      Info << "--- updated by Boundary Mapping --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
    }
    // Updated by Boundary Mapping
    scalar& rhoeUeCH_BF = rhoeUeCH.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& lambda_BF = lambda.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& chemistryOn_BF = chemistryOn.boundaryFieldRef()[currentPatchID_][faceI];

    if (debug_ && faceI==0) {
      Info << "--- already updated --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
    }
    // Already updated
    const vector nf =
        - mesh_.Sf().boundaryField()[currentPatchID_][faceI]
        / mesh_.magSf().boundaryField()[currentPatchID_][faceI];
    const scalar mDotGw_BF = mDotGw_.boundaryField()[currentPatchID_][faceI];
    const scalar p_BF = p_.boundaryField()[currentPatchID_][faceI];
    const scalar T_BF = T_.boundaryField()[currentPatchID_][faceI];

    if (debug_ && faceI==0) {
      Info << "--- will be updated --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
    }
    // Will be updated
    scalar& mDotCw_BF = mDotCw_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& blowingCorrection_BF = blowingCorrection.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& BprimeGw_BF = BprimeGw_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& BprimeCw_BF = BprimeCw_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& h_w_BF = h_w.boundaryFieldRef()[currentPatchID_][faceI];

    if (debug_ && faceI==0) {
      Info << "--- update pyrolysis gas flow rate --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
    }

    if (debug_ && faceI==0) {
      Info << "--- blowing condition --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
      Info << "--- chemistryOn_BF=" << chemistryOn_BF << " --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()"<<endl;
      Info << "--- rhoeUeCH_BF" << rhoeUeCH_BF << endl;
    }
    if(chemistryOn_BF == 1) {
      if (rhoeUeCH_BF < 1e-6) {
        rhoeUeCH_BF=1e-6;
      }
      // Blowing correction
      bool condition_ = (rhoeUeCH_BF <= 1e-4);
      if (materialDict_.isDict("Pyrolysis")) {
        word typePyro = materialDict_.subDict("Pyrolysis").lookupOrDefault<word>("PyrolysisType", "noPyrolysisSolver<specifiedPyrolysisModel>");
        typePyro.replaceAll("PyrolysisSolver<specifiedPyrolysisModel>", "");
        if (typePyro != "no") {
          condition_ = (mDotGw_BF <= 1e-6) || (rhoeUeCH_BF <= 1e-4);
        }
      }
      if (debug_ && faceI==0) {
        Info << "--- update blowing correction --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- rhoeUeCH=" << rhoeUeCH_BF << " --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- BprimeCw_BF=" << BprimeCw_BF << " --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- mDotGw_BF=" << mDotGw_BF << " --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- condition_=" << condition_ << " --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
      }

      if (((mDotGw_BF / (rhoeUeCH_BF * blowingCorrection_BF) + BprimeCw_BF )==0)||lambda_BF==0) {
        condition_=1;
      }

      if (condition_) {
        blowingCorrection_BF= 1.0;
      } else {
        blowingCorrection_BF =
            ::log
            (
                1. + 2. * lambda_BF *
                (
                    mDotGw_BF / (rhoeUeCH_BF * blowingCorrection_BF)
                    + BprimeCw_BF
                )
            ) /
            (
                2. * lambda_BF *
                (
                    mDotGw_BF / (rhoeUeCH_BF * blowingCorrection_BF)
                    + BprimeCw_BF
                )
            );
      }

      if (debug_ && faceI==0) {
        Info << "--- update BprimeCoatingMixtureG --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
      }
      // BprimeCoatingMixtureG
      BprimeGw_BF = max(mDotGw_BF / (rhoeUeCH_BF * blowingCorrection_BF),0);
    }

    if(chemistryOn_BF == 0) {
      if (debug_ && faceI==0) {
        Info << "--- update (BprimeCoatingMixtureC, mDotC, hw) == 0 --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
      }
      BprimeCw_BF = 0;
      mDotCw_BF = 0;
      h_w_BF = 0;
      continue;
    }

    if (debug_ && faceI==0) {
      Info << "--- update BprimeCoatingMixtureC and hw --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
    }

    forAll(p_Yke_, elementI) {
      p_Yke_[elementI] = Yke_ref[elementI];
    }

    forAll(p_Ykg_, elementI) {
      volScalarField& Ykgw = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(Ykg_names[elementI]));
      p_Ykg_[elementI] = Ykgw.boundaryFieldRef()[currentPatchID_][faceI];
    }

    // Compute B'c and h_w
    computeSurfaceMassBalance(*mixture_ptr, elemOnSurface, p_Yke_, p_Ykg_, T_BF, p_BF, BprimeGw_BF, BprimeCw_BF, h_w_BF, p_Xw_);

    forAll(p_Xw_, speciesI) {
      volScalarField& Xw = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(Xw_names[speciesI]));
      Xw.boundaryFieldRef()[currentPatchID_][faceI] = p_Xw_[speciesI];
    }

    if (debug_ && faceI==0) {
      Info << "--- update char ablation rate --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
    }
    // char ablation rate
    mDotCw_BF = BprimeCw_BF* rhoeUeCH_BF * blowingCorrection_BF;

    if (!recessionFlag) {
      BprimeCw_BF = 0;
      mDotCw_BF = 0;
    }

    if(!dynamicMesh_) {
      if (debug_ && faceI==0) {
        Info << "--- update (BprimeCoatingMixtureC, mDotC) == 0 --- Foam::BprimeCoatingMixtureBoundaryConditions::updateBprimeBC()" << endl;
      }
      BprimeCw_BF = 0;
      mDotCw_BF = 0;
    }
  }
}

Foam::wordList Foam::BprimeCoatingMixtureBoundaryConditions::neededFields()
{
  wordList neededFields;
  neededFields.append("p");
  neededFields.append("h_r");
  neededFields.append("Tbackground");
  neededFields.append("lambda");
  neededFields.append("rhoeUeCH");
  neededFields.append("qRad");
  neededFields.append("chemistryOn");
  return neededFields;
}

Eigen::VectorXd Foam::BprimeCoatingMixtureBoundaryConditions::equilibrium(const Mutation::Mixture& mix, const double& p, const double& T, const Eigen::VectorXd& Z_wsg, Eigen::VectorXd& X_sg)
{
  int ne = mix.nElements(); // number of elements
  Eigen::VectorXd Zx_wsg(Eigen::VectorXd::Zero(ne)); // Gaseous mass fraction (solid + gas)
  Eigen::VectorXd Zx_g(Eigen::VectorXd::Zero(ne)); // Gaseous mole fraction (gas)
  Eigen::VectorXd Z_g(Eigen::VectorXd::Zero(ne)); // Gaseous mass fraction (gas)
  const int ng = mix.nGas(); // Number of gases
  mix.convert<Mutation::Thermodynamics::YE_TO_XE>(Z_wsg.data(), Zx_wsg.data()); // Convert to mass fraction
  mix.equilibriumComposition(T, p, Zx_wsg.data(), X_sg.data(), Mutation::Thermodynamics::IN_PHASE); //IN_PHASE); //GLOBAL); // Compute the equilibrium in Mutation++
  Zx_g = mix.elementMatrix().topRows(ng).transpose()*X_sg.head(ng); // Gaseous mass fraction from the solid, liquid and gaseous mole fraction
  double m_z = Zx_g.sum(); // sum
  Zx_g /= m_z; // sum = 1
  mix.convert<Mutation::Thermodynamics::XE_TO_YE>(Zx_g.data(), Z_g.data()); // Convert to mass fraction
  return Z_g;
}

void Foam::BprimeCoatingMixtureBoundaryConditions::computeSurfaceMassBalance(const Mutation::Mixture& mix, const List<Tuple2<string,double>>& elemOnSurface, const Eigen::VectorXd& z_e, const Eigen::VectorXd& z_pg, const double T,
    const double p, const double Bg, double &Bc, double &hw, Eigen::VectorXd& X_sg)
{
  const int ns = mix.nSpecies(); // number of species
  const int ne = mix.nElements(); // number of elements
  const int ng = mix.nGas(); // number of gaseous species

  Eigen::VectorXd z_wsg(Eigen::VectorXd::Zero(ne)); // Solid, liquid and gaseous mass fraction at the wall
  Eigen::VectorXd z_wg(Eigen::VectorXd::Zero(ne)); // Gaseous mass fraction at the wall

  // Initialize the wall element fractions to be the pyrolysis gas fractions
  double sum = 0.0;
  z_wsg=z_e+Bg*z_pg;

  // Use "large" amount to simulate infinite char
  double large = max(100.0*Bg, 200.0);
  for (int eI=0; eI<elemOnSurface.size(); eI++) {
    int i = mix.elementIndex(elemOnSurface[eI].first());
    double large_elem = (elemOnSurface[eI].second() * large);
    z_wsg[i] += large_elem;
  }
  sum=z_wsg.sum();
  z_wsg/=sum;
  z_wg = equilibrium(mix,  p, T, z_wsg, X_sg); // Gaseous mass fraction at the wall

  double z_w_surf = 0;
  double z_e_surf = 0;
  double z_pg_surf = 0;

  for (int eI=0; eI<elemOnSurface.size(); eI++) {
    int i = mix.elementIndex(elemOnSurface[eI].first());
    z_w_surf += z_wg[i];
    z_e_surf += z_e[i];
    z_pg_surf += z_pg[i];
  }

  if (z_w_surf == 1) {
    Bc = maxBc;
  } else {
    Bc = (z_e_surf + Bg*z_pg_surf - z_w_surf*(1.0 + Bg)) / (z_w_surf - 1.0); // Dimensionless char blowing rate
    Bc = max(Bc, 0.0);
  }

  if (Bc>maxBc) {
    Bc = maxBc;
  }
  if (Bc<minBc) {
    Bc = minBc;
  }

  double * hwi = new double[ns]; // Wall species enthalpies
  for(int i=0; i<ns; i++) {
    hwi[i]=0;
  }
  const double RU = Mutation::RU;

  // Compute the gas enthalpy
  mix.speciesHOverRT(T, hwi);
  double mwg = 0.0;
  for (int j = 0; j < ng; ++j) {
    mwg += mix.speciesMw(j) * X_sg[j];
  }
  hw = 0.0;
  for (int i = 0; i < ng; ++i) {
    hw += X_sg[i] * hwi[i];
  }
  hw *= RU * T / mwg;
}

void Foam::BprimeCoatingMixtureBoundaryConditions::write(Ostream& os) const
{
  if (boundaryMapping_ptr) {
    os << "\t// --- start --- Boundary Mapping Inputs" << endl;
    boundaryMapping_ptr->write(os);
    os << "\t// --- end --- Boundary Mapping Inputs" << endl;
  }

  os.writeKeyword("mixtureName") << mixtureName << token::END_STATEMENT << nl;
  os.writeKeyword("environmentDirectory") << environmentDirectory << token::END_STATEMENT << nl;
  if (minBc != 1e-5) {
    os.writeKeyword("minBc") << minBc << token::END_STATEMENT << nl;
  }
  if (maxBc != 100) {
    os.writeKeyword("maxBc") << maxBc << token::END_STATEMENT << nl;
  }
}

// ************************************************************************* //

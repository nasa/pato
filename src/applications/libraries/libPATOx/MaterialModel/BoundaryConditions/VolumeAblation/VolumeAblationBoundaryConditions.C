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
    along with OpenFOAM.  If Bprimet, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "VolumeAblationBoundaryConditions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Switch Foam::VolumeAblationBoundaryConditions::ptr_deleted="no";

Foam::VolumeAblationBoundaryConditions::VolumeAblationBoundaryConditions
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
currentPatchID_(currentPatchID),
dict_(dict),
boundaryName_(boundaryName()),
energyModel_(meshLookupOrConstructModel<simpleEnergyModel>(mesh_,regionName_,simpleEnergyModel::modelName)),
initialized_(false),
dynamicMesh_(isA<dynamicFvMesh>(mesh)),
p_pT(new double [2]),
physicsBasedErosionModel_(readScalar(dict_.lookup("physicsBasedErosionModel"))),
environmentDirectory(fileName(dict_.lookup("environmentDirectory")).expand()),
environmentDictionary_
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
materialDict_(
    IOobject
    (
        IOobject::groupName(regionName+"Properties", phaseName),
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    )
),
mixtureMutation(materialDict_.subDict("MaterialChemistry").lookup("mixture")),
materialPropertiesDirectory
(
    fileName(materialDict_.subDict("MaterialProperties").lookup("MaterialPropertiesDirectory")).expand()
),
constantPropertiesDictionary
(
    IOobject
    (
        "constantProperties",
        materialPropertiesDirectory,
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
),
rfFail(dimensionedScalar::lookupOrDefault("rfFail",constantPropertiesDictionary,dimensionSet(0, 1, 0, 0, 0, 0, 0),1e-6)),
sigmaSB(::constant::physicoChemical::sigma),
initialPosition_(energyModel_.createVolField<vector>("initialPosition",dimensionedVector("0", dimLength, vector(0.0,0.0,0.0)))),
debug_(materialDict_.lookupOrDefault<Switch>("debug","no")),
neededFields_(neededFields())
{
  // Create new fields in Energy Model
  scalarFields_.insert("qRadEmission",energyModel_.createVolField<scalar>("qRadEmission",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("qRadAbsorption",energyModel_.createVolField<scalar>("qRadAbsorption",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("qConv",energyModel_.createVolField<scalar>("qConv",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("qCond",energyModel_.createVolField<scalar>("qCond",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("recession",energyModel_.createVolField<scalar>("recession", dimensionedScalar("0", dimLength, scalar(0.0))));
  scalarFields_.insert("rhoeUeCH",energyModel_.createVolField<scalar>("rhoeUeCH", dimensionedScalar("0", dimMass/ pow3(dimLength)/ dimTime, scalar(0.0))));
  scalarFields_.insert("blowingCorrection",energyModel_.createVolField<scalar>("blowingCorrection",dimensionedScalar("1", dimless, scalar(1.0))));
  scalarFields_.insert("h_w",energyModel_.createVolField<scalar>("h_w",dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0.0))));
  scalarFields_.insert("h_ew",energyModel_.createVolField<scalar>("h_ew",dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0.0))));
  scalarFields_.insert("h_r",energyModel_.createVolField<scalar>("h_r",dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0.0))));
  scalarFields_.insert("Tbackground",energyModel_.createVolField<scalar>("Tbackground", dimensionedScalar("0", dimTemperature, scalar(0.0))));
  scalarFields_.insert("lambda",energyModel_.createVolField<scalar>("lambda", dimensionedScalar("0", dimless, scalar(0.0))));
  scalarFields_.insert("qRad",energyModel_.createVolField<scalar>("qRad", dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("chemistryOn",energyModel_.createVolField<scalar>("chemistryOn", dimensionedScalar("0", dimless, scalar(0.0))));
  scalarFields_.insert("emissivity",energyModel_.createVolField<scalar>("emissivity",dimensionedScalar("0", dimless, scalar(0))));
  scalarFields_.insert("absorptivity",energyModel_.createVolField<scalar>("absorptivity",dimensionedScalar("0", dimless, scalar(0))));

  optEq_ = new Mutation::MixtureOptions(mixtureMutation);
  optEq_->setStateModel("ChemNonEq1T");
  mix_ = new Mutation::Mixture(*optEq_);
  ns_mix = mix_->nSpecies();

  Yie_ref = new double [ns_mix];
  Info << energyModel_.getTabLevel() << "Reading boundary-layer edge species composition" << nl;
  for (int j = 0; j < ns_mix ; j++) {
    Yie_ref[j] =
        environmentDictionary_.lookupOrDefault<scalar>
        (
            "Yie[" + mix_->speciesName(j) + "]",
            0.0
        );

    if (Yie_ref[j]!=0.0) Info << energyModel_.getTabLevel() << "Yie[" << mix_->speciesName(j) << "] = " <<  Yie_ref[j] << nl;
  }


  if(dynamicMesh_) {

    cellMotionU_ptr = &const_cast<volVectorField&>(mesh_.objectRegistry::lookupObject<volVectorField>("cellMotionU"));

    forAll(initialPosition_.boundaryField()[currentPatchID_], faceI) {
      initialPosition_.boundaryFieldRef()[currentPatchID_][faceI] =  mesh_.Cf().boundaryField()[currentPatchID_][faceI];
    }
  }
  energyModel_.tabLevel_--;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::VolumeAblationBoundaryConditions::~VolumeAblationBoundaryConditions()
{
  if (!ptr_deleted) { // Avoid deleting empty pointers if the destructor is called twice.
    delete[] p_pT;
    delete[] Yie_ref;
    ptr_deleted="yes";
  }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

word Foam::VolumeAblationBoundaryConditions::boundaryName()
{
  word boundary_name = "\"Temperature Boundary Conditions\" type \"VolumeAblation\"";
  word bold_on="\e[1m";
  word bold_off="\e[0m";
  Info << simpleModel::getTabLevel("| ") << bold_on << "Create " << boundary_name << bold_off << endl;
  simpleModel::tabLevel_++;
  return boundary_name;
}

void Foam::VolumeAblationBoundaryConditions::initialize()
{
  word bold_on="\e[1m";
  word bold_off="\e[0m";
  Info << simpleModel::getTabLevel("| ") << bold_on << "Initialize Temperature BC (VolumeAblation)" <<  bold_off << endl;
  simpleModel::tabLevel_++;
  // References to models
  const simpleMassModel& massModel=meshLookupOrConstructModel<simpleMassModel>(mesh_,regionName_,simpleMassModel::modelName);
  const simpleGasPropertiesModel& gasPropertiesModel=meshLookupOrConstructModel<simpleGasPropertiesModel>(mesh_,regionName_,simpleGasPropertiesModel::modelName);
  const simpleMaterialChemistryModel& materialChemistryModel=meshLookupOrConstructModel<simpleMaterialChemistryModel>(mesh_,regionName_,simpleMaterialChemistryModel::modelName);
  const simpleVolumeAblationModel& volumeAblationModel=meshLookupOrConstructModel<simpleVolumeAblationModel>(mesh_,regionName_,simpleVolumeAblationModel::modelName);
  // References to fields from Mass Model
  scalarFields_.insert("p",massModel.refVolField<scalar>("p"));
  vectorFields_.insert("mDotG",massModel.refVolField<vector>("mDotG"));
  scalarFields_.insert("mDotGw",massModel.refVolField<scalar>("mDotGw"));
  // References to fields from Gas Properties Model
  scalarFields_.insert("h_g",gasPropertiesModel.refVolField<scalar>("h_g"));
  // References to fields from Material Chemistry Model
  scalarFields_.insert("h_c",materialChemistryModel.refVolField<scalar>("h_c"));
  scalarFields_.insert("omegaHeterogeneousRate",materialChemistryModel.refVolField<scalar>("omegaHeterogeneousRate"));
  // References to fields from Volume Ablation Model
  scalarFields_.insert("rT",volumeAblationModel.refVolField<scalar>("rT"));
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

void Foam::VolumeAblationBoundaryConditions::update()
{
  if(!initialized_) {
    initialize();
  }

  if (debug_) {
    Info << "--- update BoundaryMapping --- Foam::VolumeAblationBoundaryConditions::update()" << endl;
  }
  updateBoundaryMapping();

  if (debug_) {
    Info << "--- update temperature --- Foam::VolumeAblationBoundaryConditions::update()" << endl;
  }
  updateTemperatureBC();

  if(dynamicMesh_) {
    if (debug_) {
      Info << "--- update motion --- Foam::VolumeAblationBoundaryConditions::update()" << endl;
    }
    updateMotionBC();

    if(this->debug_) {
      Info << "--- end --- Foam::VolumeAblationBoundaryConditions::update()" << endl;
    }
  }
}

void Foam::VolumeAblationBoundaryConditions::updateBoundaryMapping()
{
  // Update all the other fields (different than Ta) needed for VolumeAblationBoundaryConditions (p, rhoeUeCH_, ...)
  forAll(boundaryMapping_ptr->mappingFieldsName(), fieldI) {
    if (!isA<boundaryMappingFvPatchScalarField>(const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(boundaryMapping_ptr->mappingFieldsName()[fieldI])).boundaryField()[currentPatchID_])) {
      boundaryMapping_ptr->update(mesh_.time().value(),currentPatchID_,boundaryMapping_ptr->mappingFieldsName()[fieldI]);
    } else { // correct the fields with a different mappingBoundary patch type
      const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(boundaryMapping_ptr->mappingFieldsName()[fieldI])).correctBoundaryConditions();
    }
  }
}

void Foam::VolumeAblationBoundaryConditions::updateMotionBC()
{
  // Reference to fields
  volScalarField& rho_s_=scalarFields_["rho_s"];
  volScalarField& rT=scalarFields_["rT"];
  volScalarField& omegaHeterogeneousRate_=scalarFields_["omegaHeterogeneousRate"];
  volScalarField& recession_=scalarFields_["recession"];

  if(mesh_.objectRegistry::foundObject<volVectorField>("cellMotionU")) {
    volVectorField& cellMotionU_ = *cellMotionU_ptr;
    forAll(cellMotionU_.boundaryFieldRef()[currentPatchID_], faceI) {
      // Already updated
      const vector nf =
          - mesh_.Sf().boundaryField()[currentPatchID_][faceI]
          / mesh_.magSf().boundaryField()[currentPatchID_][faceI];
      const vector normal_ = nf;
      const scalar& rho_s_bf = rho_s_.boundaryField()[currentPatchID_][faceI];
      const tmp<volVectorField> grad_rho_s_tmp = fvc::grad(rho_s_);
      const volVectorField& grad_rho_s = grad_rho_s_tmp();
      const tmp<vectorField> grad_rho_s_int_tmp = grad_rho_s.boundaryField()[currentPatchID_].patchInternalField();
      const vectorField& grad_rho_s_int = grad_rho_s_int_tmp();
      const scalarField& grad_rho_s_int_proj = grad_rho_s_int & normal_;
      const tmp<volVectorField> grad_rT_tmp = fvc::grad(rT);
      const volVectorField& grad_rT = grad_rT_tmp();
      const tmp<vectorField> grad_rT_int_tmp = grad_rT.boundaryField()[currentPatchID_].patchInternalField();
      const vectorField& grad_rT_int = grad_rT_int_tmp();
      const scalarField& grad_rT_int_proj = grad_rT_int & normal_;
      const tmp<scalarField> rTw_tmp = rT.boundaryField()[currentPatchID_].patchInternalField();
      const scalarField& rTw = rTw_tmp();
      const scalarField& invDxw = mesh_.deltaCoeffs().boundaryField()[currentPatchID_];
      const tmp<scalarField> rho_s_int_tmp = rho_s_.boundaryField()[currentPatchID_].patchInternalField();
      const scalarField& rho_s_int = rho_s_int_tmp();
      const tmp<scalarField> omegaHeterogeneousRate_int_tmp = omegaHeterogeneousRate_.boundaryField()[currentPatchID_].patchInternalField();
      const scalarField& omegaHeterogeneousRate_int = omegaHeterogeneousRate_int_tmp();

      // Will be updated
      vector&  cellMotionU_BF = cellMotionU_.boundaryFieldRef()[currentPatchID_][faceI];

      if (physicsBasedErosionModel_ <= 0 || physicsBasedErosionModel_ >= 7) {
        FatalErrorInFunction << "Physics based erosion models (\"physicsBasedErosionModel\") are implemented from 1 to 6." << exit(FatalError);
      }
      if (physicsBasedErosionModel_ == 1) { // Rate based on min rf and grad(rf)
        if ((rTw[faceI] <= rfFail.value()) & (grad_rT_int_proj[faceI] > 0)) {
          cellMotionU_BF = (rfFail.value() - rTw[faceI]) / (grad_rT_int_proj[faceI] *  mesh_.time().deltaT().value()) * normal_;
        }
      }

      if (physicsBasedErosionModel_ == 2) { // Rate based on density gradient
        if ((rho_s_bf <= 50.0) & (grad_rho_s_int_proj[faceI] > 0.0)) {
          cellMotionU_BF = (50 - rho_s_bf) / (grad_rho_s_int_proj[faceI] *  mesh_.time().deltaT().value()) * normal_;
        }
      }

      if (physicsBasedErosionModel_ == 3) { // remove cell-by-cell
        if (rTw[faceI]  <= rfFail.value()) {
          cellMotionU_BF=  cellMotionU_BF * 0.98 + 0.02 * 1 / (invDxw[faceI] * mesh_.time().deltaT().value()) * normal_;
        } else {
          cellMotionU_BF = cellMotionU_BF * 0.98;
        }
      }

      if (physicsBasedErosionModel_ == 4) { // Test smoothing
        if ((rTw[faceI]  <= rfFail.value()) ) {
          cellMotionU_BF
            = max(1.05 * cellMotionU_BF & normal_, 1e-7 / rho_s_bf) * normal_;
        } else {
          cellMotionU_BF = cellMotionU_BF * 0.95;
        }
      }

      if (physicsBasedErosionModel_ == 5) { // Test omega
        if ((rTw[faceI]  <= rfFail.value()) ) {
          cellMotionU_BF = 0.95 * cellMotionU_BF + 0.05 * omegaHeterogeneousRate_int[faceI] / rho_s_int[faceI] * normal_;
        } else {
          cellMotionU_BF *= 0.95 ;
        }
      }

      if (physicsBasedErosionModel_ == 6) { // basic

        if ((rTw[faceI]  <= rfFail.value()) ) {
          cellMotionU_BF
            = (1 / invDxw[faceI]) * (1 / mesh_.time().deltaT().value()) * normal_;
        } else {
          cellMotionU_BF = 0.0 * normal_;
        }
      }
    }
  }

  forAll(recession_.boundaryField()[currentPatchID_], faceI) {
    recession_.boundaryFieldRef()[currentPatchID_][faceI] = mag(initialPosition_.boundaryField()[currentPatchID_][faceI] -  mesh_.Cf().boundaryField()[currentPatchID_][faceI]);
  }
}

void Foam::VolumeAblationBoundaryConditions::updateTemperatureBC()
{
  // Reference to fields
  volScalarField& T_=scalarFields_["Ta"];
  volScalarField& h_r=scalarFields_["h_r"];
  volScalarField& Tbackground=scalarFields_["Tbackground"];
  volScalarField& qRad=scalarFields_["qRad"];
  volScalarField& chemistryOn=scalarFields_["chemistryOn"];
  volScalarField& lambda=scalarFields_["lambda"];
  volScalarField& rhoeUeCH_=scalarFields_["rhoeUeCH"];
  volScalarField& p_=scalarFields_["p"];
  volScalarField& emissivity_=scalarFields_["emissivity"];
  volScalarField& absorptivity_=scalarFields_["absorptivity"];
  volTensorField& k_=tensorFields_["k"];
  volScalarField& mDotGw_=scalarFields_["mDotGw"];
  volScalarField& h_ew=scalarFields_["h_ew"];
  volScalarField& qRadEmission_=scalarFields_["qRadEmission"];
  volScalarField& qRadAbsorption_=scalarFields_["qRadAbsorption"];
  volScalarField& qConv_=scalarFields_["qConv"];
  volScalarField& qCond_=scalarFields_["qCond"];
  volScalarField& blowingCorrection_=scalarFields_["blowingCorrection"];

  forAll(T_.boundaryFieldRef()[currentPatchID_], faceI) {

    // Updated by BoundaryMapping
    scalar& h_r_BF = h_r.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& Tbackground_BF = Tbackground.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qRad_BF = qRad.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& chemistryOn_BF = chemistryOn.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& lambda_BF = lambda.boundaryFieldRef()[currentPatchID_][faceI];

    // Already updated
    const scalar& rhoeUeCH_BF = rhoeUeCH_.boundaryFieldRef()[currentPatchID_][faceI];
    const vector nf =
        - mesh_.Sf().boundaryField()[currentPatchID_][faceI]
        / mesh_.magSf().boundaryField()[currentPatchID_][faceI];
    const scalar invDx_BF = mesh_.deltaCoeffs().boundaryField()[currentPatchID_][faceI];
    const scalar p_BF = p_.boundaryField()[currentPatchID_][faceI];
    const tmp<scalarField> Tint_tmp = T_.boundaryField()[currentPatchID_].patchInternalField();
    const scalarField& Tint_ = Tint_tmp();
    const scalar Tint_BF = Tint_[faceI];
    const scalar emissivity_BF = emissivity_.boundaryField()[currentPatchID_][faceI];
    const scalar absorptivity_BF = absorptivity_.boundaryField()[currentPatchID_][faceI];
    const tensor k_BF = k_.boundaryField()[currentPatchID_][faceI];
    const scalar kProj_BF = nf & k_BF & nf; // Projection of k on the surface normal
    const scalar mDotGw_BF = mDotGw_.boundaryField()[currentPatchID_][faceI];

    // Will be updated
    scalar& h_ew_BF = h_ew.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qRadEmission_BF = qRadEmission_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qRadAbsorption_BF = qRadAbsorption_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qConv_BF = qConv_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qCond_BF = qCond_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& T_BF = T_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& blowingCorrection_BF = blowingCorrection_.boundaryFieldRef()[currentPatchID_][faceI];

    // update h_ew_BF
    p_pT[0] = p_BF;
    p_pT[1] = T_BF;
    mix_->setState(Yie_ref, p_pT, 2); // set state using the composition of the environment
    h_ew_BF = mix_->mixtureHMass();

    if (debug_) {
      Info << "--- blowing condition --- Foam::VolumeAblationBoundaryConditions::updateBprimeBC()" << endl;
      Info << "--- chemistryOn_BF=" << chemistryOn_BF << " --- Foam::VolumeAblationBoundaryConditions::updateBprimeBC()"<<endl;
      Info << "--- rhoeUeCH_BF" << rhoeUeCH_BF << endl;
    }
    if(chemistryOn_BF == 1) {
      // Blowing correction
      bool condition_ = (rhoeUeCH_BF <= 1e-4);
      if (materialDict_.isDict("Pyrolysis")) {
        word typePyro = materialDict_.subDict("Pyrolysis").lookupOrDefault<word>("PyrolysisType", "noPyrolysisSolver<specifiedPyrolysisModel>");
        typePyro.replaceAll("PyrolysisSolver<specifiedPyrolysisModel>", "");
        if (typePyro != "no") {
          condition_ = (mDotGw_BF <= 1e-6) || (rhoeUeCH_BF <= 1e-4);
        }
      }
      if (debug_) {
        Info << "--- update blowing correction --- Foam::VolumeAblationBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- rhoeUeCH=" << rhoeUeCH_BF << " --- Foam::VolumeAblationBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- mDotGw_BF=" << mDotGw_BF << " --- Foam::VolumeAblationBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- condition_=" << condition_ << " --- Foam::VolumeAblationBoundaryConditions::updateBprimeBC()" << endl;
      }


      if (rhoeUeCH_BF == 0) {
        FatalErrorInFunction << "rhoeUeCH_BF == 0 for rhoeUeCH\n" << rhoeUeCH_ << exit(FatalError);
      }
      if ((mDotGw_BF / (rhoeUeCH_BF * blowingCorrection_BF))==0) {
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
                )
            ) /
            (
                2. * lambda_BF *
                (
                    mDotGw_BF / (rhoeUeCH_BF * blowingCorrection_BF)
                )
            );
      }
    } // if chemistryOn_BF == 1
    else {
      blowingCorrection_BF= 1.0;
    }

    // Surface energy balance: heat fluxes and temperature
    qRadEmission_BF = - emissivity_BF * sigmaSB.value() * (pow4(T_BF) - pow4(Tbackground_BF));
    qRadAbsorption_BF = absorptivity_BF * qRad_BF;

    if(chemistryOn_BF == 1) {

      qConv_BF = rhoeUeCH_BF * blowingCorrection_BF * (h_r_BF - h_ew_BF) ;

      T_BF =
          Tint_BF
          + (
              1. / (kProj_BF * invDx_BF) *
              (
                  qConv_BF  + qRadEmission_BF + qRadAbsorption_BF
              )
          );
    } else {
      qConv_BF = 0;
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

Foam::wordList Foam::VolumeAblationBoundaryConditions::neededFields()
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

void Foam::VolumeAblationBoundaryConditions::write(Ostream& os) const
{
  os << "\t// --- start --- Boundary Mapping Inputs" << endl;
  boundaryMapping_ptr->write(os);
  os << "\t// --- end --- Boundary Mapping Inputs" << endl;
  os.writeKeyword("environmentDirectory") << environmentDirectory << token::END_STATEMENT << nl;
  os.writeKeyword("physicsBasedErosionModel") << physicsBasedErosionModel_ << token::END_STATEMENT << nl;

}


// ************************************************************************* //

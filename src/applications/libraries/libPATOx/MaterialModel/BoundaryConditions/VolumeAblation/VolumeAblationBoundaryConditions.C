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

Foam::VolumeAblationBoundaryConditions::VolumeAblationBoundaryConditions
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName,
    const label& currentPatchID,
    const dictionary dict
)
  :
mesh_(mesh),
phaseName_(phaseName),
dictName_(dictName),
currentPatchID_(currentPatchID),
dict_(dict),
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
qRadEmission_(meshLookupOrConstructScalar(mesh, "qRadEmission",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
qRadAbsorption_(meshLookupOrConstructScalar(mesh, "qRadAbsorption",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
qConv_(meshLookupOrConstructScalar(mesh, "qConv",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
qCond_(meshLookupOrConstructScalar(mesh, "qCond",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
recession_(meshLookupOrConstructScalar(mesh, "recession", dimensionedScalar("0", dimLength, scalar(0.0)))),
rhoeUeCH_(meshLookupOrConstructScalar(mesh, "rhoeUeCH", dimensionedScalar("0", dimMass/ pow3(dimLength)/ dimTime, scalar(0.0)))),
blowingCorrection_(meshLookupOrConstructScalar(mesh, "blowingCorrection",dimensionedScalar("1", dimless, scalar(1.0)))),
h_w(meshLookupOrConstructScalar(mesh, "h_w",dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0.0)))),
h_ew(meshLookupOrConstructScalar(mesh, "h_ew",dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0.0)))),
h_r(meshLookupOrConstructScalar(mesh, "h_r",dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0.0)))),
Tbackground(meshLookupOrConstructScalar(mesh, "Tbackground", dimensionedScalar("0", dimTemperature, scalar(0.0)))),
lambda(meshLookupOrConstructScalar(mesh, "lambda", dimensionedScalar("0", dimless, scalar(0.0)))),
qRad(meshLookupOrConstructScalar(mesh, "qRad", dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
heatOn(meshLookupOrConstructScalar(mesh, "heatOn", dimensionedScalar("0", dimless, scalar(0.0)))),
//update in Mass Model
mDotG_(meshLookupOrConstructVector(mesh, "mDotG", dimensionedVector("zero", dimMass/(dimLength*dimLength*dimTime), vector(0, 0, 0)))),
mDotGFace_(meshLookupOrConstructSurfaceVector(mesh, "mDotGface", linearInterpolate(mDotG_))),
//update in Gas Properties Model
h_g_(meshLookupOrConstructScalar(mesh, "h_g", dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0))),
//update in Material Properties Model
h_c_(meshLookupOrConstructScalar(mesh, "h_c", dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0))),
emissivity_(meshLookupOrConstructScalar(mesh, "emissivity", dimensionedScalar("0", dimless, scalar(0)))),
absorptivity_(meshLookupOrConstructScalar(mesh, "absorptivity", dimensionedScalar("0", dimless, scalar(0)))),
k_(meshLookupOrConstructTensor(mesh, "k", dimensionedTensor("0", dimensionSet(1, 1, -3, -1, 0, 0, 0), tensor(1,0,0,0,1,0,0,0,1)))),
rT(meshLookupOrConstructScalar(mesh, "rT", dimensionedScalar("1", dimLength, 1.0))),
omegaHeterogeneousRate_(meshLookupOrConstructScalar(mesh, "omegaHeterogeneousRate", dimensionedScalar("0", dimMass/pow3(dimLength)/dimTime, 0.0))),
// MUST_READ fields
rho_s_(meshLookupOrConstructScalar(mesh, "rho_s")),
p_(meshLookupOrConstructScalar(mesh, "p")),
T_(meshLookupOrConstructScalar(mesh, "Ta")),
materialDict_(
    IOobject
    (
        IOobject::groupName(dictName+"Properties", phaseName),
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
initialPosition_
(
    meshLookupOrConstructVector(mesh,  "initialPosition",    dimensionedVector("0", dimLength, vector(0.0,0.0,0.0)))
),
debug_(materialDict_.lookupOrDefault<Switch>("debug","no")),
neededFields_(neededFields()),
boundaryMapping_
(
    simpleBoundaryMappingModel::New(
        mesh_,
        neededFields_,
        dict_
    )
),
boundaryMapping_ptr(&boundaryMapping_())
{
  optEq_ = new Mutation::MixtureOptions(mixtureMutation);
  optEq_->setStateModel("ChemNonEq1T");
  mix_ = new Mutation::Mixture(*optEq_);
  ns_mix = mix_->nSpecies();

  Yie_ref = new double [ns_mix];
  Info << "Reading boundary-layer edge species composition" << nl;
  for (int j = 0; j < ns_mix ; j++) {
    Yie_ref[j] =
        environmentDictionary_.lookupOrDefault<scalar>
        (
            "Yie[" + mix_->speciesName(j) + "]",
            0.0
        );

    Info << "Yie[" << mix_->speciesName(j) << "] = " <<  Yie_ref[j] << nl;
  }


  if(dynamicMesh_) {

    cellMotionU_ptr =
        new volVectorField
    (
        IOobject
        (
            "cellMotionU",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh
    );

    forAll(initialPosition_.boundaryField()[currentPatchID_], faceI) {
      initialPosition_.boundaryFieldRef()[currentPatchID_][faceI] =  mesh_.Cf().boundaryField()[currentPatchID_][faceI];
    }
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::VolumeAblationBoundaryConditions::~VolumeAblationBoundaryConditions()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::VolumeAblationBoundaryConditions::update()
{
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
  // OpenFOAM 4.x
  if(mesh_.objectRegistry::foundObject<volVectorField>("cellMotionU")) {
    volVectorField& cellMotionU_ = const_cast<volVectorField&>(meshLookupOrConstructVector(mesh_, "cellMotionU"));
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
  forAll(T_.boundaryFieldRef()[currentPatchID_], faceI) {

    // Updated by BoundaryMapping
    scalar& h_r_BF = h_r.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& Tbackground_BF = Tbackground.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qRad_BF = qRad.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& heatOn_BF = heatOn.boundaryFieldRef()[currentPatchID_][faceI];
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
    const scalar mDotGFace_BF = mDotGFace_.boundaryField()[currentPatchID_][faceI]&(-nf);

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
      Info << "--- heatOn_BF=" << heatOn_BF << " --- Foam::VolumeAblationBoundaryConditions::updateBprimeBC()"<<endl;
      Info << "--- rhoeUeCH_BF" << rhoeUeCH_BF << endl;
    }
    if(heatOn_BF == 1) {
      // Blowing correction
      bool condition_ = (rhoeUeCH_BF <= 1e-4);
      if (materialDict_.isDict("Pyrolysis")) {
        word typePyro = materialDict_.subDict("Pyrolysis").lookupOrDefault<word>("PyrolysisType", "noPyrolysisSolver<specifiedPyrolysisModel>");
        typePyro.replaceAll("PyrolysisSolver<specifiedPyrolysisModel>", "");
        if (typePyro != "no") {
          condition_ = (mDotGFace_BF <= 1e-6) || (rhoeUeCH_BF <= 1e-4);
        }
      }
      if (debug_) {
        Info << "--- update blowing correction --- Foam::VolumeAblationBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- rhoeUeCH=" << rhoeUeCH_BF << " --- Foam::VolumeAblationBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- mDotGFace_BF=" << mDotGFace_BF << " --- Foam::VolumeAblationBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- condition_=" << condition_ << " --- Foam::VolumeAblationBoundaryConditions::updateBprimeBC()" << endl;
      }


      if (rhoeUeCH_BF == 0) {
        FatalErrorInFunction << "rhoeUeCH_BF == 0 for rhoeUeCH\n" << rhoeUeCH_ << exit(FatalError);
      }
      if ((mDotGFace_BF / (rhoeUeCH_BF * blowingCorrection_BF))==0) {
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
                    mDotGFace_BF / (rhoeUeCH_BF * blowingCorrection_BF)
                )
            ) /
            (
                2. * lambda_BF *
                (
                    mDotGFace_BF / (rhoeUeCH_BF * blowingCorrection_BF)
                )
            );
      }
    } // if heatOn_BF == 1
    else {
      blowingCorrection_BF= 1.0;
    }

    // Surface energy balance: heat fluxes and temperature
    qRadEmission_BF = - emissivity_BF * sigmaSB.value() * (pow4(T_BF) - pow4(Tbackground_BF));
    qRadAbsorption_BF = absorptivity_BF * qRad_BF;

    if(heatOn_BF == 1) {

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
  neededFields.append("heatOn");
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

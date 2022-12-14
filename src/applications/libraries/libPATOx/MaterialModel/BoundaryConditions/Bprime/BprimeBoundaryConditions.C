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

#include "BprimeBoundaryConditions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BprimeBoundaryConditions::BprimeBoundaryConditions
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
mixtureMutationBprime(dict_.lookupOrDefault<word>("mixtureMutationBprime","none")),
environmentDirectory(fileName(dict_.lookupOrDefault<fileName>("environmentDirectory","none")).expand()),
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
moleFractionGasInMaterial_(dict_.lookupOrDefault<List<Tuple2<word,scalar> > >("moleFractionGasInMaterial",List<Tuple2<word,scalar> >(0))),
writeBprimeTable_(dict_.lookupOrDefault<Switch>("writeBprimeTable","no")),
minBc(dict_.lookupOrDefault<scalar>("BprimeTableBcmin",1e-5)),
maxBc(dict_.lookupOrDefault<scalar>("BprimeTableBcmax",100)),
readBprimeTable_("no"),
mDotGw_(meshLookupOrConstructScalar(mesh, "mDotGw",dimensionedScalar("0", dimMass/pow(dimLength,2)/dimTime, scalar(0.0)))),
mDotCw_(meshLookupOrConstructScalar(mesh, "mDotCw", dimensionedScalar("0", dimMass/pow(dimLength,2)/dimTime, scalar(0.0)))),
BprimeCw_(meshLookupOrConstructScalar(mesh, "BprimeCw",dimensionedScalar("0", dimless, scalar(0.0)))),
BprimeGw_(meshLookupOrConstructScalar(mesh, "BprimeGw",dimensionedScalar("0", dimless, scalar(0.0)))),
qRadEmission_(meshLookupOrConstructScalar(mesh, "qRadEmission",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
qRadAbsorption_(meshLookupOrConstructScalar(mesh, "qRadAbsorption",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
qAdvPyro_(meshLookupOrConstructScalar(mesh, "qAdvPyro",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
qAdvChar_(meshLookupOrConstructScalar(mesh, "qAdvChar",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
qConv_(meshLookupOrConstructScalar(mesh, "qConv",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
qCond_(meshLookupOrConstructScalar(mesh, "qCond",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
recession_(meshLookupOrConstructScalar(mesh, "recession", dimensionedScalar("0", dimLength, scalar(0.0)))),
rhoeUeCH(meshLookupOrConstructScalar(mesh, "rhoeUeCH", dimensionedScalar("0", dimMass/ pow3(dimLength)/ dimTime, scalar(0.0)))),
blowingCorrection(meshLookupOrConstructScalar(mesh, "blowingCorrection",dimensionedScalar("1", dimless, scalar(1.0)))),
h_w(meshLookupOrConstructScalar(mesh, "h_w",dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0.0)))),
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
boundaryMapping_ptr(&boundaryMapping_()),
optionsBC_(mesh_, dict_, currentPatchID_),
ne_Bprime(0),
ns_Bprime(0),
recessionFlag(dict_.lookupOrDefault<Switch>("recessionFlag","yes")),
nameBprimeFile_(dict_.lookupOrDefault<fileName>("BprimeFile",fileName(word::null).expand()))
{
  if(debug_) {
    Info << "--- environmentDirectory=" << environmentDirectory << " --- Foam::BprimeBoundaryConditions::BprimeBoundaryConditions" << endl;
    Info << "--- dict_=" << dict_ << " --- Foam::BprimeBoundaryConditions::BprimeBoundaryConditions" << endl;
  }
  if (nameBprimeFile_==word::null) {
    /*
    fileName dirName_(fileName(dict_.lookup("environmentDirectory")).expand());
    environmentDirectory = dirName_;
    if (!isFile(environmentDirectory+"/environmentComposition") ) {
      FatalErrorInFunction << (word) environmentDirectory<< "/environmentComposition not found."
                           << exit(FatalError);
    }

    IOdictionary dictEnvi_
    (
        IOobject
        (
            "environmentComposition",
            environmentDirectory,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );
    environmentDictionary=dictEnvi_;
    */
  } else {
    if (mixtureMutationBprime != "none") {
      FatalErrorInFunction << "Please remove \"mixtureMutationBprime\" or \"BprimeFile\"" << exit(FatalError);
    }
  }

  initBprimeMixture();
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


Foam::BprimeBoundaryConditions::~BprimeBoundaryConditions()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::BprimeBoundaryConditions::initBprimeMixture()
{
  if (writeBprimeTable_ && mixtureMutationBprime=="none") {
    FatalErrorInFunction << "If you want to write the Bprime tables, please specify \"mixtureMutationBprime\"" << exit(FatalError);
  }

  if (nameBprimeFile_!=word::null && !writeBprimeTable_) {
    readBprimeTable_="yes";
    readBprimeTable();
  } else {
    word mixName_(dict_.lookup("mixtureMutationBprime"));
    mixtureMutationBprime = mixName_;
    mixBprime_ptr= new Mutation::Mixture(mixtureMutationBprime);
    Mutation::Mixture& mixBprime_ = *mixBprime_ptr ;

    ne_Bprime= mixBprime_.nElements();
    ns_Bprime= mixBprime_.nSpecies();
    Yke_ref = new double[ne_Bprime];
    Ykg_ref= new double[ne_Bprime];

    // Initialization of Mutation++ for surface equilibrium chemistry (B')
    // Default conditions (input in mole fractions -> converted in mass fractions)
    Info << "Reading boundary-layer edge elemental composition (mass fraction)" << nl;
    for (int i = 0; i < ne_Bprime; i++) {
      Yke_ref[i] =
          environmentDictionary.lookupOrDefault<scalar>
          (
              "Yke["+  mixBprime_.elementName(i) + "]",
              0.0
          );

    }

    mixBprime_.convert<Mutation::Thermodynamics::XE_TO_YE>(Yke_ref, Yke_ref);
    for (int i = 0; i < ne_Bprime; i++) {
      Info << "Yke[" <<  mixBprime_.elementName(i) << "] = " << Yke_ref[i] << nl;
    }

    Info << "Reading initial elemental composition of the gas inside the material for B' solver (mass fraction)" << nl;


    if (moleFractionGasInMaterial_.size()>0) {
      Info << "Use \"Zx[...]\" from \"moleFractionGasInMaterial\" instead of \"" << (word) constantPropertiesDictionary.path() << "/constantProperties\"." << endl;

      scalar sumZx_ = 0;
      forAll(moleFractionGasInMaterial_, gasI) {
        word name_ = moleFractionGasInMaterial_[gasI].first();
        scalar value_ =  moleFractionGasInMaterial_[gasI].second();
        sumZx_ += value_;
        Info << "Zx[\"" << name_ << "\"] = " << value_ << endl;
        name_ = "Zx["+name_+"]";
        constantPropertiesDictionary.add(name_ ,value_);
      }
      if (sumZx_ != 1) {
        FatalErrorInFunction << "The sum of \"moleFractionGasInMaterial\" has to be 1." << exit(FatalError);
      }
    }
    bool atLeastOneElement = false;
    for (int i = 0; i < ne_Bprime; i++) {
      Ykg_ref[i] =
          constantPropertiesDictionary.lookupOrDefault<scalar>
          (
              "Zx[" +  mixBprime_.elementName(i) + "]",
              0.0
          );
      if (Ykg_ref[i]!=0) {
        atLeastOneElement = true;
      }
    }

    if (!atLeastOneElement) {
      FatalErrorInFunction << "Zx[...] not found in \"" <<   (word) materialPropertiesDirectory << "/constantProperties\"" << exit(FatalError);
    }

    mixBprime_.convert<Mutation::Thermodynamics::XE_TO_YE>(Ykg_ref, Ykg_ref);

    for (int i = 0; i < ne_Bprime; i++) {
      Info << "Ykg[" <<  mixBprime_.elementName(i) << "] = " << Ykg_ref[i] << nl;
    }

    for (int elementI=0; elementI<ne_Bprime; elementI++) {
      word nameYkg_ = "Ykg["+mixBprime_.elementName(elementI)+"]";
      Ykg_names.append(nameYkg_);
      volScalarField * Ykg_ = new volScalarField(meshLookupOrConstructScalar(mesh_, nameYkg_.c_str(),dimensionedScalar("0", dimless,  Ykg_ref[elementI]),"fixedValue"));
      Ykg_->store(Ykg_);
    }
    for (int speciesI=0; speciesI<ns_Bprime; speciesI++) {
      word nameXw_ = "Xw["+mixBprime_.speciesName(speciesI)+"]";
      Xw_names.append(nameXw_);
      volScalarField * Xw_ = new volScalarField(meshLookupOrConstructScalar(mesh_, nameXw_.c_str(),dimensionedScalar("0", dimless, scalar(0)),"fixedValue"));
      Xw_->store(Xw_);
    }
  }

  if (writeBprimeTable_) {
    writeBprimeTable();
    nameBprimeFile_ ="BPrimeTable-PATO" ;
    readBprimeTable_="yes";
    readBprimeTable();
  }
}

void Foam::BprimeBoundaryConditions::update()
{
  if (debug_) {
    Info << "--- update BoundaryMapping --- Foam::BprimeBoundaryConditions::update()" << endl;
  }
  updateBoundaryMapping();

  if (debug_) {
    Info << "--- update optionsBC --- Foam::BprimeBoundaryConditions::update()" << endl;
  }
  optionsBC_.update();

  if (debug_) {
    Info << "--- update Bprime --- Foam::BprimeBoundaryConditions::update()" << endl;
  }
  updateBprimeBC();
  if (debug_) {
    Info << "--- update temperature --- Foam::BprimeBoundaryConditions::update()" << endl;
  }
  updateTemperatureBC();

  if(dynamicMesh_) {
    if (debug_) {
      Info << "--- update motion --- Foam::BprimeBoundaryConditions::update()" << endl;
    }
    updateMotionBC();

    if(this->debug_) {
      Info << "--- end --- Foam::BprimeBoundaryConditions::update()" << endl;
    }
  }
}

void Foam::BprimeBoundaryConditions::updateBoundaryMapping()
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

void Foam::BprimeBoundaryConditions::updateMotionBC()
{
  // OpenFOAM 4.x
  if(mesh_.objectRegistry::foundObject<volVectorField>("cellMotionU")) {
    volVectorField& cellMotionU_ = const_cast<volVectorField&>(meshLookupOrConstructVector(mesh_, "cellMotionU"));
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
    recession_.boundaryFieldRef()[currentPatchID_][faceI] = mag(initialPosition_.boundaryField()[currentPatchID_][faceI] -  mesh_.Cf().boundaryField()[currentPatchID_][faceI]);
  }
}

void Foam::BprimeBoundaryConditions::updateTemperatureBC()
{
  forAll(T_.boundaryFieldRef()[currentPatchID_], faceI) {
    // not used yet
    scalar wallSpeciesDiffusion = 1;

    // Updated by BoundaryMapping
    scalar& h_r_BF = h_r.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& Tbackground_BF = Tbackground.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qRad_BF = qRad.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& heatOn_BF = heatOn.boundaryFieldRef()[currentPatchID_][faceI];

    // Already updated
    const scalar& rhoeUeCH_BF = rhoeUeCH.boundaryFieldRef()[currentPatchID_][faceI];
    const vector nf =
        - mesh_.Sf().boundaryField()[currentPatchID_][faceI]
        / mesh_.magSf().boundaryField()[currentPatchID_][faceI];
    const scalar invDx_BF = mesh_.deltaCoeffs().boundaryField()[currentPatchID_][faceI];
    const tmp<scalarField> Tint_tmp = T_.boundaryField()[currentPatchID_].patchInternalField();
    const scalarField& Tint_ = Tint_tmp();
    const scalar Tint_BF = Tint_[faceI];
    const scalar mDotGFace_BF = mDotGFace_.boundaryField()[currentPatchID_][faceI]&(-nf);
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

    if(heatOn_BF == 1) {
      qAdvPyro_BF = mDotGFace_BF * (h_g_BF - h_w_BF);
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

void Foam::BprimeBoundaryConditions::updateBprimeBC()
{
  if (debug_) {
    Info << "--- currentPatchID_=" << currentPatchID_ << " --- Foam::BprimeBoundaryConditions::updateBprimeBC()"  << endl;
    Info << "--- update Boundary Mapping for patch " << mesh_.boundaryMesh()[currentPatchID_].name() << " --- Foam::BprimeBoundaryConditions::updateBprimeBC()"  << endl;
  }

  forAll(BprimeCw_.boundaryField()[currentPatchID_], faceI) {

    // Species mole fraction at the wall
    scalarList p_Xw_(ns_Bprime);
    // Elements mass fraction inside the material
    scalarList p_Ykg_(ne_Bprime);
    // Elements mass fraction in the environment
    scalarList p_Yke_(ne_Bprime);

    if (debug_ && faceI==0) {
      Info << "--- updated by Boundary Mapping --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
    }
    // Updated by Boundary Mapping
    scalar& rhoeUeCH_BF = rhoeUeCH.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& lambda_BF = lambda.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& heatOn_BF = heatOn.boundaryFieldRef()[currentPatchID_][faceI];

    if (debug_ && faceI==0) {
      Info << "--- already updated --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
    }
    // Already updated
    const vector nf =
        - mesh_.Sf().boundaryField()[currentPatchID_][faceI]
        / mesh_.magSf().boundaryField()[currentPatchID_][faceI];
    const scalar mDotGFace_BF = mDotGFace_.boundaryField()[currentPatchID_][faceI]&(-nf);
    const scalar p_BF = p_.boundaryField()[currentPatchID_][faceI];
    const scalar T_BF = T_.boundaryField()[currentPatchID_][faceI];

    if (debug_ && faceI==0) {
      Info << "--- will be updated --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
    }
    // Will be updated
    scalar& mDotGw_BF = mDotGw_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& mDotCw_BF = mDotCw_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& blowingCorrection_BF = blowingCorrection.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& BprimeGw_BF = BprimeGw_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& BprimeCw_BF = BprimeCw_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& h_w_BF = h_w.boundaryFieldRef()[currentPatchID_][faceI];

    if (debug_ && faceI==0) {
      Info << "--- update pyrolysis gas flow rate --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
    }
    // pyrolysis gas flow rate
    mDotGw_BF = mDotGFace_BF;

    if (debug_ && faceI==0) {
      Info << "--- blowing condition --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
      Info << "--- heatOn_BF=" << heatOn_BF << " --- Foam::BprimeBoundaryConditions::updateBprimeBC()"<<endl;
      Info << "--- rhoeUeCH_BF" << rhoeUeCH_BF << endl;
    }
    if(heatOn_BF == 1) {
      if (rhoeUeCH_BF < 1e-6) {
        rhoeUeCH_BF=1e-6;
      }
      // Blowing correction
      bool condition_ = (rhoeUeCH_BF <= 1e-4);
      if (materialDict_.isDict("Pyrolysis")) {
        word typePyro = materialDict_.subDict("Pyrolysis").lookupOrDefault<word>("PyrolysisType", "noPyrolysisSolver<specifiedPyrolysisModel>");
        typePyro.replaceAll("PyrolysisSolver<specifiedPyrolysisModel>", "");
        if (typePyro != "no") {
          condition_ = (mDotGFace_BF <= 1e-6) || (rhoeUeCH_BF <= 1e-4);
        }
      }
      if (debug_ && faceI==0) {
        Info << "--- update blowing correction --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- rhoeUeCH=" << rhoeUeCH_BF << " --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- BprimeCw_BF=" << BprimeCw_BF << " --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- mDotGFace_BF=" << mDotGFace_BF << " --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
        Info << "--- condition_=" << condition_ << " --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
      }

      if (((mDotGFace_BF / (rhoeUeCH_BF * blowingCorrection_BF) + BprimeCw_BF )==0)||lambda_BF==0) {
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
                    + BprimeCw_BF
                )
            ) /
            (
                2. * lambda_BF *
                (
                    mDotGFace_BF / (rhoeUeCH_BF * blowingCorrection_BF)
                    + BprimeCw_BF
                )
            );
      }

      if (debug_ && faceI==0) {
        Info << "--- update BprimeG --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
      }
      // BprimeG
      BprimeGw_BF = max(mDotGFace_BF / (rhoeUeCH_BF * blowingCorrection_BF),0);
    }

    if(heatOn_BF == 0) {
      if (debug_ && faceI==0) {
        Info << "--- update (BprimeC, mDotC, hw) == 0 --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
      }
      BprimeCw_BF = 0;
      mDotCw_BF = 0;
      h_w_BF = 0;
      continue;
    }

    if (debug_ && faceI==0) {
      Info << "--- update BprimeC and hw --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
    }

    forAll(p_Yke_, elementI) {
      p_Yke_[elementI] = Yke_ref[elementI];
    }

    forAll(p_Ykg_, elementI) {
      volScalarField& Ykgw = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(Ykg_names[elementI]));
      p_Ykg_[elementI] = Ykgw.boundaryFieldRef()[currentPatchID_][faceI];
    }

    // BprimeC and h_w
    computeSurfaceMassBalance(p_Yke_, p_Ykg_, p_BF, T_BF, BprimeGw_BF, BprimeCw_BF, h_w_BF, p_Xw_);

    forAll(p_Xw_, speciesI) {
      volScalarField& Xw = const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(Xw_names[speciesI]));
      Xw.boundaryFieldRef()[currentPatchID_][faceI] = p_Xw_[speciesI];
    }

    if (debug_ && faceI==0) {
      Info << "--- update char ablation rate --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
    }
    // char ablation rate
    mDotCw_BF = BprimeCw_BF* rhoeUeCH_BF * blowingCorrection_BF;

    if (!recessionFlag) {
      BprimeCw_BF = 0;
      mDotCw_BF = 0;
    }

    if(!dynamicMesh_) {
      if (debug_ && faceI==0) {
        Info << "--- update (BprimeC, mDotC) == 0 --- Foam::BprimeBoundaryConditions::updateBprimeBC()" << endl;
      }
      BprimeCw_BF = 0;
      mDotCw_BF = 0;
    }
  }
}

Foam::wordList Foam::BprimeBoundaryConditions::neededFields()
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

void Foam::BprimeBoundaryConditions::writeBprimeTable()
{
  // Elements mass fraction inside the material
  scalarList p_Ykg_(ne_Bprime);
  // Elements mass fraction in the environment
  scalarList p_Yke_(ne_Bprime);

  forAll(p_Yke_, elementI) {
    p_Yke_[elementI] = Yke_ref[elementI];
  }

  forAll(p_Ykg_, elementI) {
    p_Ykg_[elementI] = Ykg_ref[elementI];
  }

  scalarList BprimeTablePw_(dict_.lookup("BprimeTablePw"));
  scalar minT(readScalar(dict_.lookup("BprimeTableTmin")));
  scalar maxT(readScalar(dict_.lookup("BprimeTableTmax")));
  scalar numberBc(readScalar(dict_.lookup("BprimeTableBcNumber")));
  scalarList BprimeTableBg_(dict_.lookup("BprimeTableBg"));

  if(minBc <= 0) {
    FatalErrorInFunction << "minBc has to be greater than 0." << exit(FatalError);
  }
  if(maxBc <= minBc) {
    FatalErrorInFunction << "maxBc <= minBc" << exit(FatalError);
  }
  if (minT >= maxT) {
    FatalErrorInFunction << "BprimeTableTmin >= BprimeTableTmax" << exit(FatalError);
  }

  if (numberBc <= 3) {
    FatalErrorInFunction << "BprimeTableBcNumber <= 3" << exit(FatalError);
  }

  // sort and reverse p and B'g
  sort(BprimeTablePw_); // atm
  sort(BprimeTableBg_); // [-]
  reverse(BprimeTablePw_);
  reverse(BprimeTableBg_);

  // B'c and T outputs
  scalarList BprimeTableBc_(numberBc);
  scalarList BprimeTableTw_(numberBc);
  scalarList BprimeTableHw_(numberBc);

  scalar minSlope = 1e-10;

  // high number of B'c and T for the inverse function
  scalar numberT = 1000;
  scalarList inverseTw_(numberT);
  scalarList inverseBc_(numberT);
  scalarList inverseHw_(numberT);
  scalar deltaT = (maxT - minT )/(numberT-1);
  forAll(inverseTw_, TI) {
    inverseTw_[TI]=minT+deltaT*TI; // K
  }

  word outputFileName = "BPrimeTable";
  Info << "Writing " << outputFileName << "-PATO and " << outputFileName << "-FIAT" <<endl;
  Info << "N_P = " << BprimeTablePw_.size() << ", N_T = " << BprimeTableTw_.size() << ", N_Bg = " << BprimeTableBg_.size() << endl;
  OFstream outputPATO_(outputFileName+"-PATO");
  OFstream outputFIAT_(outputFileName+"-FIAT");

  outputPATO_ << "// p(Pa)       B'g (-)         B'c (-)         T (K)           h_w (J/kg)" << endl;
  outputFIAT_ << "//p(atm) B'g(-) B'c(-) B'f(-) Tsurf(K) Dexp(-) Zhw(-) h_w(cal/g) IABL(-) mat(-)" << endl;

  scalar convHw = 4184; // cal/g to J/kg
  scalar ONEATM = 101325; // atm to Pa

  scalarList p_Xw(ns_Bprime);
  for(int i=0; i<ns_Bprime; i++) {
    p_Xw[i]=0;
  }
  scalar Tw = 0;
  scalar pw = 0;
  scalar Bg = 0;

  forAll(BprimeTablePw_, pwi) {
    forAll(BprimeTableBg_, Bgi) {
      forAll(inverseTw_, Twi) {
        for(int i=0; i<ns_Bprime; i++) {
          p_Xw[i]=0;
        }
        Tw = inverseTw_[Twi];
        pw = BprimeTablePw_[pwi]*ONEATM;
        Bg = BprimeTableBg_[Bgi];

        computeSurfaceMassBalance(p_Yke_, p_Ykg_, pw, Tw, Bg, inverseBc_[Twi], inverseHw_[Twi], p_Xw);

        if (Twi > 0) {
          if ((inverseBc_[Twi]- inverseBc_[Twi-1])/deltaT < minSlope) {
            if ((abs(inverseBc_[Twi]- inverseBc_[Twi-1])/deltaT) > 1e-6) {
              Info << "(pw,Bg,Tw) = (" << BprimeTablePw_[pwi] << "," <<  BprimeTableBg_[Bgi] << "," << inverseTw_[Twi] <<
                   ") \t (inverseBc_[Twi]- inverseBc_[Twi-1])/deltaT = " << (inverseBc_[Twi]- inverseBc_[Twi-1])/deltaT << endl;
            }
            inverseBc_[Twi] = inverseBc_[Twi-1] + deltaT * minSlope;
          }
        }

        if (inverseBc_[Twi] >= maxBc || Twi == numberT - 1) {
          scalarList inverseTw_crop;
          scalarList inverseBc_crop;
          for(int i = 0 ; i <= Twi ; i++) {
            inverseTw_crop.append(inverseTw_[i]);
            inverseBc_crop.append(inverseBc_[i]);
          }

          scalar deltaBc_log = (log(inverseBc_[Twi]) - log(inverseBc_[0]) )/ (numberBc-1);

          forAll(BprimeTableBc_, BcI) {
            BprimeTableBc_[BcI] = exp(log(inverseBc_[0])+deltaBc_log*BcI);
            BprimeTableTw_[BcI] = linearInterpolation(log(inverseBc_crop), inverseTw_crop,log(BprimeTableBc_[BcI]));
            BprimeTableHw_[BcI]=0;
            computeSurfaceMassBalance(p_Yke_, p_Ykg_, pw, BprimeTableTw_[BcI], Bg, BprimeTableBc_[BcI], BprimeTableHw_[BcI], p_Xw);

            outputPATO_ << pw << " " << Bg << " " << BprimeTableBc_[BcI] << " " <<  BprimeTableTw_[BcI] << " " << BprimeTableHw_[BcI] << endl;
            int IABL = 0;
            word mat = "char";

            if (BprimeTableBc_[BcI] > minBc) {
              IABL = 1;
              mat = "C(gr)*__";
            }

            double hw_cal_g = BprimeTableHw_[BcI]/convHw;

            outputFIAT_.setf(ios_base::fixed, ios_base::floatfield);
            outputFIAT_.precision(5);
            outputFIAT_ << setw(12)  << BprimeTablePw_[pwi];
            outputFIAT_ << setw(12)  << Bg;
            outputFIAT_ << setw(12)<< BprimeTableBc_[BcI];
            outputFIAT_ << setw(12)  << "0.00000";
            outputFIAT_ << setw(12)  <<  BprimeTableTw_[BcI] ;
            outputFIAT_ << setw(5) <<  "0" ;
            outputFIAT_.precision(3);
            outputFIAT_ << setw(10)  << hw_cal_g;
            outputFIAT_ << setw(10) << hw_cal_g;
            outputFIAT_ << setw(2) << IABL;
            outputFIAT_ << " " <<  mat << endl;
            outputFIAT_.unsetf(ios_base::right);

          }
          break;
        }

      }

    }
  }
}

void Foam::BprimeBoundaryConditions::readBprimeTable()
{
  // read and store the B' table into the RAM for faster access
  Info << "Reading the B' Table" << nl;
  fileName BprimeFileName(changeEnviVar(nameBprimeFile_));
  if(debug_) {
    Info << "--- BprimeFileName=" << BprimeFileName << endl;
  }
//  fileName BprimeFileName(fileName(nameBprimeFile_).expand());
  IFstream BprimeInputFile(BprimeFileName.c_str());  // opens an input file

  if (BprimeInputFile.good() == false) { // checks the input file is opened
    FatalErrorInFunction << "Bprime table not found. Either update the path or change the simulation options."
                         << exit(FatalError); // exit
  }

  IFstream BprimeInputFileTemp(BprimeFileName.c_str());  // open a second input file
  Info << "The Bprime table must have 5 columns: p, B'g, B'c, T, Hg, ..." << nl;

  int columnTableBp = 5;
  int rawTableBp = 0;
  scalar tempBp, rawTableFracBp, rawTableIntBp;
  int i_rawBp = 0;

  while (true) {
    BprimeInputFileTemp >> tempBp;
    if (BprimeInputFileTemp.eof() == 1)
      break;
    i_rawBp++;
  }
  rawTableFracBp =
      modf
      (
          static_cast<scalar>(i_rawBp) /
          static_cast<scalar>(columnTableBp),
          &rawTableIntBp
      );

  if (rawTableFracBp != 0) {
    FatalErrorInFunction << "... and it does not seem to have 5 columns." << nl
                         << rawTableIntBp << " lines and "
                         << rawTableFracBp * columnTableBp << " 'scalar' have been read."
                         << exit(FatalError);
  } else {
    rawTableBp = i_rawBp / columnTableBp;
    Info << "The Bprime table has been read (" << rawTableBp << " lines)."
         << nl;
  }

  RectangularMatrix<scalar> Table(rawTableBp, columnTableBp, 0);

  for (int x = 0; x < rawTableBp; x++) {
    for (int i = 0; i < columnTableBp; i++) {
      BprimeInputFile >> Table[x][i];
    }
  }

  // BprimeT table is an object of the class BprimeSI, which handles the B' table interpolations when called
  BprimeT_ptr = new BprimeTable(Table);

  Info << "The B' Table has been read" << nl;
}

void Foam::BprimeBoundaryConditions::computeSurfaceMassBalance(scalarList Yke_, scalarList Ykg_, scalar pw_, scalar Tw_, scalar Bg_, scalar& Bc_, scalar& hw_, scalarList& p_Xw_)
{
  if (debug_) {
    Info << "--- compute Surface Mass Balance" << endl;
  }

  double pw = pw_;
  double Tw = Tw_;
  double Bg = Bg_;
  double Bc = 0;
  double hw = 0;
  double * Yke = new double[ne_Bprime];
  double * Ykg = new double[ne_Bprime];
  for(int i = 0; i < ne_Bprime; i++) {
    Yke[i] = Yke_[i];
    Ykg[i] = Ykg_[i];
  }
  double * Xw = new double[ns_Bprime];
  for(int i=0; i<ns_Bprime; i++) {
    Xw[i]=0;
  }

  BprimeTable& BprimeT_ = *BprimeT_ptr;
  if(readBprimeTable_) {
    if (debug_) {
      Info << "--- update Bc from " << nameBprimeFile_ << endl;
    }
    Bc = BprimeT_.BprimeC(pw, Bg, Tw);
    if (debug_) {
      Info << "--- update hw from " << nameBprimeFile_ << endl;
    }
    hw = BprimeT_.hw();
  } else {
    if (debug_) {
      Info << "--- update Bc from Mutation++ with \"" << mixtureMutationBprime << "\" mixture" << endl;
    }
    Mutation::Mixture& mixBprime_ = *mixBprime_ptr;
    mixBprime_.surfaceMassBalance(Yke, Ykg, Tw, pw, Bg, Bc, hw, Xw);

    if (debug_) {
      Info << "--- update Bc for pure C" << endl;
    }
    // For char only, infinity amount of C
    if (Bc>=199.99 && Tw < 1700) {
      Bc=0;
    }
    if (Bc>=199.99 && Tw >= 1700) {
      Bc=maxBc;
    }
    if (Bc>=maxBc) {
      Bc=maxBc;
    }
    if (Bc < minBc) {
      Bc =  minBc;
    }
  }

  if (debug_) {
    Info << "--- update X_w" << endl;
  }
  p_Xw_.resize(ns_Bprime);
  for(int i=0; i<ns_Bprime; i++) {
    p_Xw_[i] = Xw[i];
  }
  Bc_=Bc;
  hw_=hw;
}

void Foam::BprimeBoundaryConditions::write(Ostream& os) const
{
  if (boundaryMapping_ptr) {
    os << "\t// --- start --- Boundary Mapping Inputs" << endl;
    boundaryMapping_ptr->write(os);
    os << "\t// --- end --- Boundary Mapping Inputs" << endl;
  }

  if (mixtureMutationBprime!="none") {
    os.writeKeyword("mixtureMutationBprime") << mixtureMutationBprime << token::END_STATEMENT << nl;
  }
  if (environmentDirectory!="none") {
    os.writeKeyword("environmentDirectory") << environmentDirectory << token::END_STATEMENT << nl;
  }
  if (moleFractionGasInMaterial_.size()>0) {
    os.writeKeyword("moleFractionGasInMaterial") << moleFractionGasInMaterial_ << token::END_STATEMENT << nl;
  }
  if (writeBprimeTable_) {
    os.writeKeyword("writeBprimeTable") << writeBprimeTable_ << token::END_STATEMENT << nl;
  }
  if (minBc != 1e-5) {
    os.writeKeyword("BprimeTableBcmin") << minBc << token::END_STATEMENT << nl;
  }
  if (maxBc != 100) {
    os.writeKeyword("BprimeTableBcmax") << maxBc << token::END_STATEMENT << nl;
  }
  if (nameBprimeFile_ != fileName(word::null).expand()) {
    os.writeKeyword("BprimeFile") << nameBprimeFile_ << token::END_STATEMENT << nl;
  }
  if(dict_.found("BprimeTablePw")) {
    scalarList BprimeTablePw_(dict_.lookup("BprimeTablePw"));
    os.writeKeyword("BprimeTablePw") << BprimeTablePw_ << token::END_STATEMENT << nl;
  }
  if(dict_.found("BprimeTableTmin")) {
    scalar minT(readScalar(dict_.lookup("BprimeTableTmin")));
    os.writeKeyword("BprimeTableTmin") << minT << token::END_STATEMENT << nl;
  }
  if(dict_.found("BprimeTableTmax")) {
    scalar maxT(readScalar(dict_.lookup("BprimeTableTmax")));
    os.writeKeyword("BprimeTableTmax") << maxT << token::END_STATEMENT << nl;
  }
  if(dict_.found("BprimeTableBcNumber")) {
    scalar numberBc(readScalar(dict_.lookup("BprimeTableBcNumber")));
    os.writeKeyword("BprimeTableBcNumber") << numberBc << token::END_STATEMENT << nl;
  }
  if(dict_.found("BprimeTableBg")) {
    scalarList BprimeTableBg_(dict_.lookup("BprimeTableBg"));
    os.writeKeyword("BprimeTableBg") << BprimeTableBg_ << token::END_STATEMENT << nl;
  }
}

// ************************************************************************* //

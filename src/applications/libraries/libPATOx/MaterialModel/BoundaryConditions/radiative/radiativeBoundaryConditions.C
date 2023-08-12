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
    along with OpenFOAM.  If radiativet, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "radiativeBoundaryConditions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiativeBoundaryConditions::radiativeBoundaryConditions
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
energyModel_(meshLookupOrConstructModel<simpleEnergyModel>(mesh_,regionName_,simpleEnergyModel::modelName)),
initialized_(false),
sigmaSB(::constant::physicoChemical::sigma),
neededFields_(neededFields())
{
  // Create new fields in Energy Model
  scalarFields_.insert("qRadEmission",energyModel_.createVolField<scalar>("qRadEmission",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("qRadAbsorption",energyModel_.createVolField<scalar>("qRadAbsorption",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("qConv",energyModel_.createVolField<scalar>("qConv",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("qCond",energyModel_.createVolField<scalar>("qCond",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("rhoeUeCH",energyModel_.createVolField<scalar>("rhoeUeCH", dimensionedScalar("0", dimMass/ pow3(dimLength)/ dimTime, scalar(0.0))));
  scalarFields_.insert("blowingCorrection",energyModel_.createVolField<scalar>("blowingCorrection",dimensionedScalar("1", dimless, scalar(1.0))));
  scalarFields_.insert("h_r",energyModel_.createVolField<scalar>("h_r",dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0.0))));
  scalarFields_.insert("Tbackground",energyModel_.createVolField<scalar>("Tbackground", dimensionedScalar("0", dimTemperature, scalar(0.0))));
  scalarFields_.insert("lambda",energyModel_.createVolField<scalar>("lambda", dimensionedScalar("0", dimless, scalar(0.0))));
  scalarFields_.insert("qRad",energyModel_.createVolField<scalar>("qRad", dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0))));
  scalarFields_.insert("chemistryOn",energyModel_.createVolField<scalar>("chemistryOn", dimensionedScalar("0", dimless, scalar(0.0))));
  scalarFields_.insert("emissivity",energyModel_.createVolField<scalar>("emissivity",dimensionedScalar("0", dimless, scalar(0))));
  scalarFields_.insert("absorptivity",energyModel_.createVolField<scalar>("absorptivity",dimensionedScalar("0", dimless, scalar(0))));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::radiativeBoundaryConditions::~radiativeBoundaryConditions()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::radiativeBoundaryConditions::initialize()
{
  word bold_on="\e[1m";
  word bold_off="\e[0m";
  Info << simpleModel::getTabLevel("| ") << bold_on << "Initialize Temperature BC (radiative)" <<  bold_off << endl;
  simpleModel::tabLevel_++;
  // References to fields from Energy Model
  scalarFields_.insert("Ta",energyModel_.refVolField<scalar>("Ta"));
  tensorFields_.insert("k",energyModel_.refVolField<tensor>("k"));
  // Boundary Mapping
  boundaryMapping_=simpleBoundaryMappingModel::New(mesh_, neededFields_,dict_);
  boundaryMapping_ptr=&boundaryMapping_();
  // Initialized flag
  initialized_=true;
  simpleModel::tabLevel_--;
}

void Foam::radiativeBoundaryConditions::update()
{
  if(!initialized_) {
    initialize();
  }
  if (debug_) {
    Info << "--- update BoundaryMapping --- Foam::radiativeBoundaryConditions::update()" << endl;
  }
  updateBoundaryMapping();

  if (debug_) {
    Info << "--- update temperature --- Foam::radiativeBoundaryConditions::update()" << endl;
  }
  updateTemperatureBC();
}

void Foam::radiativeBoundaryConditions::updateBoundaryMapping()
{
  // Update all the other fields (different than Ta) needed for radiativeBoundaryConditions (p, rhoeUeCH, ...)
  forAll(boundaryMapping_ptr->mappingFieldsName(), fieldI) {
    if (!isA<boundaryMappingFvPatchScalarField>(const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(boundaryMapping_ptr->mappingFieldsName()[fieldI])).boundaryField()[currentPatchID_])) {
      boundaryMapping_ptr->update(mesh_.time().value(),currentPatchID_,boundaryMapping_ptr->mappingFieldsName()[fieldI]);
    } else { // correct the fields with a different mappingBoundary patch type
      const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(boundaryMapping_ptr->mappingFieldsName()[fieldI])).correctBoundaryConditions();
    }
  }
}

void Foam::radiativeBoundaryConditions::updateTemperatureBC()
{
  // Reference to fields
  volScalarField& T_=scalarFields_["Ta"];
  volScalarField& h_r=scalarFields_["h_r"];
  volScalarField& Tbackground=scalarFields_["Tbackground"];
  volScalarField& qRad=scalarFields_["qRad"];
  volScalarField& rhoeUeCH_=scalarFields_["rhoeUeCH"];
  volScalarField& emissivity_=scalarFields_["emissivity"];
  volScalarField& absorptivity_=scalarFields_["absorptivity"];
  volTensorField& k_=tensorFields_["k"];
  volScalarField& qRadEmission_=scalarFields_["qRadEmission"];
  volScalarField& qRadAbsorption_=scalarFields_["qRadAbsorption"];
  volScalarField& qConv_=scalarFields_["qConv"];
  volScalarField& qCond_=scalarFields_["qCond"];

  forAll(T_.boundaryFieldRef()[currentPatchID_], faceI) {
    // Updated by BoundaryMapping
    scalar& h_r_BF = h_r.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& Tbackground_BF = Tbackground.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qRad_BF = qRad.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& emissivity_BF = emissivity_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& absorptivity_BF = absorptivity_.boundaryFieldRef()[currentPatchID_][faceI];
    emissivity_BF = 0.9;
    absorptivity_BF = 0.9;

    // Already updated
    const scalar& rhoeUeCH_BF = rhoeUeCH_.boundaryFieldRef()[currentPatchID_][faceI];
    const vector nf =
        - mesh_.Sf().boundaryField()[currentPatchID_][faceI]
        / mesh_.magSf().boundaryField()[currentPatchID_][faceI];
    const scalar invDx_BF = mesh_.deltaCoeffs().boundaryField()[currentPatchID_][faceI];
    const tmp<scalarField> Tint_tmp = T_.boundaryField()[currentPatchID_].patchInternalField();
    const scalarField& Tint_ = Tint_tmp();
    const scalar Tint_BF = Tint_[faceI];
    const tensor k_BF = k_.boundaryField()[currentPatchID_][faceI];
    const scalar kProj_BF = nf & k_BF & nf; // Projection of k on the surface normal

    // Will be updated
    scalar& qRadEmission_BF = qRadEmission_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qRadAbsorption_BF = qRadAbsorption_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qConv_BF = qConv_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qCond_BF = qCond_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& T_BF = T_.boundaryFieldRef()[currentPatchID_][faceI];

    // Surface energy balance: heat fluxes and temperature
    qRadEmission_BF = - emissivity_BF * sigmaSB.value() * (pow4(T_BF) - pow4(Tbackground_BF));
    qRadAbsorption_BF = absorptivity_BF * qRad_BF;
    qConv_BF = rhoeUeCH_BF  * h_r_BF ;
    T_BF =
        Tint_BF
        + (
            1. / (kProj_BF * invDx_BF) *
            (
                qConv_BF + qRadEmission_BF + qRadAbsorption_BF
            )
        );
    qCond_BF = (kProj_BF * invDx_BF)* (T_BF - Tint_BF);
  }



}

Foam::wordList Foam::radiativeBoundaryConditions::neededFields()
{
  wordList neededFields;
  neededFields.append("h_r");
  neededFields.append("Tbackground");
  neededFields.append("rhoeUeCH");
  neededFields.append("qRad");
  return neededFields;
}

void Foam::radiativeBoundaryConditions::write(Ostream& os) const
{
  os << "\t// --- start --- Boundary Mapping Inputs" << endl;
  boundaryMapping_ptr->write(os);
  os << "\t// --- end --- Boundary Mapping Inputs" << endl;
}





// ************************************************************************* //

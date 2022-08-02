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

#include "qRad_emission_absorptionBoundaryConditions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::qRad_emission_absorptionBoundaryConditions::qRad_emission_absorptionBoundaryConditions
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
qRadEmission_(meshLookupOrConstructScalar(mesh, "qRadEmission",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
qRadAbsorption_(meshLookupOrConstructScalar(mesh, "qRadAbsorption",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
qCond_(meshLookupOrConstructScalar(mesh, "qCond",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
Tbackground(meshLookupOrConstructScalar(mesh, "Tbackground", dimensionedScalar("0", dimTemperature, scalar(0.0)))),
qRad(meshLookupOrConstructScalar(mesh, "qRad", dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
emissivity_(meshLookupOrConstructScalar(mesh, "emissivity", dimensionedScalar("0", dimless, scalar(0)))),
absorptivity_(meshLookupOrConstructScalar(mesh, "absorptivity", dimensionedScalar("0", dimless, scalar(0)))),
k_(meshLookupOrConstructTensor(mesh, "k", dimensionedTensor("0", dimensionSet(1, 1, -3, -1, 0, 0, 0), tensor(1,0,0,0,1,0,0,0,1)))),
T_(meshLookupOrConstructScalar(mesh, "Ta")),
p_(meshLookupOrConstructScalar(mesh, "p")),
sigmaSB(::constant::physicoChemical::sigma),
qRad_factor_(dict_.lookupOrDefault<scalar>("qRad_factor",1)),
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::qRad_emission_absorptionBoundaryConditions::~qRad_emission_absorptionBoundaryConditions()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::qRad_emission_absorptionBoundaryConditions::update()
{
  if (debug_) {
    Info << "--- update BoundaryMapping --- Foam::qRad_emission_absorptionBoundaryConditions::update()" << endl;
  }
  updateBoundaryMapping();

  if (debug_) {
    Info << "--- update temperature --- Foam::qRad_emission_absorptionBoundaryConditions::update()" << endl;
  }
  updateTemperatureBC();
}

void Foam::qRad_emission_absorptionBoundaryConditions::updateBoundaryMapping()
{
  // Update all the other fields (different than Ta) needed for qRad_emission_absorptionBoundaryConditions (p, rhoeUeCH, ...)
  forAll(boundaryMapping_ptr->mappingFieldsName(), fieldI) {
    if (!isA<boundaryMappingFvPatchScalarField>(const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(boundaryMapping_ptr->mappingFieldsName()[fieldI])).boundaryField()[currentPatchID_])) {
      boundaryMapping_ptr->update(mesh_.time().value(),currentPatchID_,boundaryMapping_ptr->mappingFieldsName()[fieldI]);
    } else { // correct the fields with a different mappingBoundary patch type
      const_cast<volScalarField&>(mesh_.objectRegistry::lookupObject<volScalarField>(boundaryMapping_ptr->mappingFieldsName()[fieldI])).correctBoundaryConditions();
    }
  }
}

void Foam::qRad_emission_absorptionBoundaryConditions::updateTemperatureBC()
{
  forAll(T_.boundaryFieldRef()[currentPatchID_], faceI) {
    // Updated by BoundaryMapping
    scalar& Tbackground_BF = Tbackground.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qRad_BF = qRad.boundaryFieldRef()[currentPatchID_][faceI];

    // Already updated
    const vector nf =
        - mesh_.Sf().boundaryField()[currentPatchID_][faceI]
        / mesh_.magSf().boundaryField()[currentPatchID_][faceI];
    const scalar invDx_BF = mesh_.deltaCoeffs().boundaryField()[currentPatchID_][faceI];
    const tmp<scalarField> Tint_tmp = T_.boundaryField()[currentPatchID_].patchInternalField();
    const scalarField& Tint_ = Tint_tmp();
    const scalar Tint_BF = Tint_[faceI];
    const tensor k_BF = k_.boundaryField()[currentPatchID_][faceI];
    const scalar kProj_BF = nf & k_BF & nf; // Projection of k on the surface normal
    const scalar& emissivity_BF = emissivity_.boundaryFieldRef()[currentPatchID_][faceI];
    const scalar& absorptivity_BF = absorptivity_.boundaryFieldRef()[currentPatchID_][faceI];

    // Will be updated
    scalar& qRadEmission_BF = qRadEmission_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qRadAbsorption_BF = qRadAbsorption_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& qCond_BF = qCond_.boundaryFieldRef()[currentPatchID_][faceI];
    scalar& T_BF = T_.boundaryFieldRef()[currentPatchID_][faceI];

    // Surface energy balance: heat fluxes and temperature
    qRadEmission_BF = - emissivity_BF * sigmaSB.value() * (pow4(T_BF) - pow4(Tbackground_BF));
    qRadAbsorption_BF = qRad_factor_ * absorptivity_BF * qRad_BF;
    T_BF =
        Tint_BF
        + (
            1. / (kProj_BF * invDx_BF) *
            (
                qRadEmission_BF + qRadAbsorption_BF
            )
        );
    qCond_BF = (kProj_BF * invDx_BF)* (T_BF - Tint_BF);
  }



}

Foam::wordList Foam::qRad_emission_absorptionBoundaryConditions::neededFields()
{
  wordList neededFields;
  neededFields.append("Tbackground");
  neededFields.append("qRad");
  neededFields.append("p");
  return neededFields;
}

void Foam::qRad_emission_absorptionBoundaryConditions::write(Ostream& os) const
{
  os << "\t// --- start --- Boundary Mapping Inputs" << endl;
  boundaryMapping_ptr->write(os);
  os << "\t// --- end --- Boundary Mapping Inputs" << endl;
}





// ************************************************************************* //

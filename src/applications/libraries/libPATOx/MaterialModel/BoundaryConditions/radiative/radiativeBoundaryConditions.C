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
qConv_(meshLookupOrConstructScalar(mesh, "qConv",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
qCond_(meshLookupOrConstructScalar(mesh, "qCond",dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
rhoeUeCH(meshLookupOrConstructScalar(mesh, "rhoeUeCH", dimensionedScalar("0", dimMass/ pow3(dimLength)/ dimTime, scalar(0.0)))),
h_r(meshLookupOrConstructScalar(mesh, "h_r",dimensionedScalar("0", pow(dimLength,2)/pow(dimTime,2), scalar(0.0)))),
Tbackground(meshLookupOrConstructScalar(mesh, "Tbackground", dimensionedScalar("0", dimTemperature, scalar(0.0)))),
lambda(meshLookupOrConstructScalar(mesh, "lambda", dimensionedScalar("0", dimless, scalar(0.0)))),
qRad(meshLookupOrConstructScalar(mesh, "qRad", dimensionedScalar("0", dimMass/ pow3(dimTime)/ dimLength, scalar(0.0)))),
emissivity_(meshLookupOrConstructScalar(mesh, "emissivity", dimensionedScalar("0", dimless, scalar(0)))),
absorptivity_(meshLookupOrConstructScalar(mesh, "absorptivity", dimensionedScalar("0", dimless, scalar(0)))),
k_(meshLookupOrConstructTensor(mesh, "k", dimensionedTensor("0", dimensionSet(1, 1, -3, -1, 0, 0, 0), tensor(1,0,0,0,1,0,0,0,1)))),
T_(meshLookupOrConstructScalar(mesh, "Ta")),
sigmaSB(::constant::physicoChemical::sigma),
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


Foam::radiativeBoundaryConditions::~radiativeBoundaryConditions()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiativeBoundaryConditions::update()
{
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
    const scalar& rhoeUeCH_BF = rhoeUeCH.boundaryFieldRef()[currentPatchID_][faceI];
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

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
    along with OpenFOAM.  If pyro_recessiont, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pyro_recessionBoundaryConditions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pyro_recessionBoundaryConditions::pyro_recessionBoundaryConditions
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
dynamicMesh_(isA<dynamicFvMesh>(mesh)),
energyModel_(meshLookupOrConstructModel<simpleEnergyModel>(mesh_,regionName_,simpleEnergyModel::modelName)),
neededFields_( {"Ta","recessionRate"}),
               debug_(energyModel_.materialDict().lookupOrDefault<Switch>("debug","no"))
{
  // Create new fields in Energy Model
  scalarFields_.insert("recessionRate",energyModel_.createVolField<scalar>("recessionRate", dimensionedScalar("0", dimLength/dimTime, scalar(0.0))));
  scalarFields_.insert("recession",energyModel_.createVolField<scalar>("recession", dimensionedScalar("0", dimLength, scalar(0.0))));
  vectorFields_.insert("initialPosition",energyModel_.createVolField<vector>("initialPosition", dimensionedVector("0", dimLength, vector(0.0,0.0,0.0))));

  if(!dynamicMesh_) {
    FatalError << "This boundary works only with dynamicMesh." << exit(FatalError);
  }

  cellMotionU_ptr = &const_cast<volVectorField&>(mesh_.objectRegistry::lookupObject<volVectorField>("cellMotionU"));

  volVectorField& initialPosition_ = vectorFields_["initialPosition"];
  forAll(initialPosition_.boundaryField()[currentPatchID_], faceI) {
    initialPosition_.boundaryFieldRef()[currentPatchID_][faceI] =  mesh_.Cf().boundaryField()[currentPatchID_][faceI];
  }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::pyro_recessionBoundaryConditions::~pyro_recessionBoundaryConditions()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pyro_recessionBoundaryConditions::updateBoundaryMapping()
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

void Foam::pyro_recessionBoundaryConditions::update()
{
  // Reference to fields
  volScalarField& recessionRate_=scalarFields_["recessionRate"];
  volScalarField& recession_=scalarFields_["recession"];
  volVectorField& initialPosition_=vectorFields_["initialPosition"];

  if (debug_) {
    Info << "--- update BoundaryMapping --- Foam::pyro_recessionBoundaryConditions::update()" << endl;
  }
  updateBoundaryMapping();
  if (debug_) {
    Info << "--- update cellMotionU --- Foam::pyro_recessionBoundaryConditions::update()" << endl;
  }
  volVectorField& cellMotionU_ = *cellMotionU_ptr;
  forAll(cellMotionU_.boundaryField()[currentPatchID_], faceI) {
    const vector nf =
        - mesh_.Sf().boundaryField()[currentPatchID_][faceI]
        / mesh_.magSf().boundaryField()[currentPatchID_][faceI];
    vector&  cellMotionU_BF = cellMotionU_.boundaryFieldRef()[currentPatchID_][faceI];
    cellMotionU_BF= recessionRate_.boundaryField()[currentPatchID_][faceI] * nf;
    recession_.boundaryFieldRef()[currentPatchID_][faceI] = mag(initialPosition_.boundaryField()[currentPatchID_][faceI] -  mesh_.Cf().boundaryField()[currentPatchID_][faceI]);
  }

}

void Foam::pyro_recessionBoundaryConditions::write(Ostream& os) const
{
  if (boundaryMapping_ptr) {
    os << "\t// --- start --- Boundary Mapping Inputs" << endl;
    boundaryMapping_ptr->write(os);
    os << "\t// --- end --- Boundary Mapping Inputs" << endl;
  }
}

// ************************************************************************* //

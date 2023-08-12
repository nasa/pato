/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "tractionDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
  :
fixedGradientFvPatchVectorField(p, iF),
traction_(p.size(), Zero),
pressure_(p.size(), 0.0)
{
  fvPatchVectorField::operator=(patchInternalField());
  gradient() = Zero;
}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
  :
fixedGradientFvPatchVectorField(p, iF),
traction_("traction", dict, p.size()),
pressure_("pressure", dict, p.size())
{
  fvPatchVectorField::operator=(patchInternalField());
  gradient() = Zero;
}

tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
  :
fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
traction_(mapper(tdpvf.traction_)),
pressure_(mapper(tdpvf.pressure_))
{}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf
)
  :
fixedGradientFvPatchVectorField(tdpvf),
traction_(tdpvf.traction_),
pressure_(tdpvf.pressure_)
{}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
  :
fixedGradientFvPatchVectorField(tdpvf, iF),
traction_(tdpvf.traction_),
pressure_(tdpvf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
  fixedGradientFvPatchVectorField::autoMap(m);
  m(traction_, traction_);
  m(pressure_, pressure_);
}


void tractionDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
  fixedGradientFvPatchVectorField::rmap(ptf, addr);

  const tractionDisplacementFvPatchVectorField& dmptf =
      refCast<const tractionDisplacementFvPatchVectorField>(ptf);

  traction_.rmap(dmptf.traction_, addr);
  pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void tractionDisplacementFvPatchVectorField::updateCoeffs()
{
  if (updated()) {
    return;
  }

  const fvPatchField<scalar>& mu =
      patch().lookupPatchField<volScalarField, scalar>("mu_sM");
  const fvPatchField<scalar>& lambda =
      patch().lookupPatchField<volScalarField, scalar>("lambda_sM");
  vectorField n(patch().nf());
  const fvPatchField<tensor>& gradField =
      patch().lookupPatchField<volTensorField, tensor>("gradD");

  gradient() =
      (traction_ - n*pressure_)
      - (n & (mu * (twoSymm(gradField) - gradField) - (mu + lambda)*gradField))
      - n*lambda*tr(gradField);

  // Thermal contribution
  const fvPatchField<scalar>&  threeKalpha =
      patch().lookupPatchField<volScalarField, scalar>
      (
          "threeKalpha"
      );

  const fvPatchScalarField& T =
      patch().lookupPatchField<volScalarField, scalar>("Ta");

  const fvPatchScalarField& T_0 =
      patch().lookupPatchField<volScalarField, scalar>("T0");

  gradient() += n*threeKalpha*(T - T_0);

  // Pyrolysis contribution
  const fvPatchField<scalar>&  threeKxi =
      patch().lookupPatchField<volScalarField, scalar>
      (
          "threeKxi"
      );
  const fvPatchScalarField& tau =
      patch().lookupPatchField<volScalarField, scalar>("tau");

  gradient() -= n*threeKxi*(1 - tau);

  gradient() /= (2.0*mu + lambda);

  fixedGradientFvPatchVectorField::updateCoeffs();
}

// Evaluate the patch field
void tractionDisplacementFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
  if (!this->updated()) {
    this->updateCoeffs();
  }

  const fvPatchField<tensor>& gradField =
      patch().lookupPatchField<volTensorField, tensor>("gradD");

  vectorField n(patch().nf());
  vectorField delta(patch().delta());

  //- non-orthogonal correction vectors
  vectorField k(delta - n*(n&delta));

  vectorField::operator=
  (
      this->patchInternalField()
      + (k & gradField.patchInternalField())
      + gradient()/this->patch().deltaCoeffs()
  );

  fvPatchVectorField::evaluate();
}

void tractionDisplacementFvPatchVectorField::write(Ostream& os) const
{
  fvPatchVectorField::write(os);
  writeEntry(os, "traction", traction_);
  writeEntry(os, "pressure", pressure_);
  writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    tractionDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

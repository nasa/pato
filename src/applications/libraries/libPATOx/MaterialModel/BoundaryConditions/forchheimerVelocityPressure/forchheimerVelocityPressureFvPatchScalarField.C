/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2020 PATO
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "forchheimerVelocityPressureFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

forchheimerVelocityPressureFvPatchScalarField::forchheimerVelocityPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
  :
fixedGradientFvPatchScalarField(p, iF),
U_(0)
{}


forchheimerVelocityPressureFvPatchScalarField::forchheimerVelocityPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
  :
fixedGradientFvPatchScalarField(p, iF)
{
  U_ = vectorField("U", dict, p.size());

  if (dict.found("gradient")) {
    gradient() = scalarField("gradient", dict, p.size());
    fixedGradientFvPatchScalarField::updateCoeffs();
    fixedGradientFvPatchScalarField::evaluate();
  } else {
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
  }
}

forchheimerVelocityPressureFvPatchScalarField::forchheimerVelocityPressureFvPatchScalarField
(
    const forchheimerVelocityPressureFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
  :
fixedGradientFvPatchScalarField(psf, p, iF, mapper),
U_(psf.U_)
{}





forchheimerVelocityPressureFvPatchScalarField::forchheimerVelocityPressureFvPatchScalarField
(
    const forchheimerVelocityPressureFvPatchScalarField& psf
)
  :
fixedGradientFvPatchScalarField(psf),
U_(psf.U_)
{}


forchheimerVelocityPressureFvPatchScalarField::forchheimerVelocityPressureFvPatchScalarField
(
    const forchheimerVelocityPressureFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
  :
fixedGradientFvPatchScalarField(psf, iF),
U_(psf.U_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void forchheimerVelocityPressureFvPatchScalarField::updateCoeffs()
{
  if (updated()) {
    return;
  }

  const Field<vector> nf(patch().nf());

  const Field<scalar>& mu =
      patch().template lookupPatchField<volScalarField, scalar>("mu_g");

  const Field<tensor>& K =
      patch().template lookupPatchField<volTensorField, tensor>("K");

  const Field<tensor>& Beta =
      patch().template lookupPatchField<volTensorField, tensor>("Beta");

  const Field<scalar>& rho =
      patch().template lookupPatchField<volScalarField, scalar>("rho_g");

  const uniformDimensionedVectorField& g =
      this->db().template lookupObject<uniformDimensionedVectorField>("g");

  gradient() =
      - (mu/((K & nf) & nf) + rho*((Beta & nf) & nf)*mag(U_))*(U_ & nf)
      + rho*(g.value() & this->patch().Cf());

  fixedGradientFvPatchScalarField::updateCoeffs();
}


void forchheimerVelocityPressureFvPatchScalarField::write(Ostream& os) const
{
  fvPatchScalarField::write(os);
  writeEntry(os, "U", U_);
  writeEntry(os, "gradient", gradient());
  writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, forchheimerVelocityPressureFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

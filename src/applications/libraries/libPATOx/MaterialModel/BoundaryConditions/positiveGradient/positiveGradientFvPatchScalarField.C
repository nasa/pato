/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "positiveGradientFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

positiveGradientFvPatchScalarField::positiveGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
  :
fixedGradientFvPatchScalarField(p, iF)
{}


positiveGradientFvPatchScalarField::positiveGradientFvPatchScalarField
(
    const positiveGradientFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
  :
fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


positiveGradientFvPatchScalarField::positiveGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
  :
fixedGradientFvPatchScalarField(p, iF)
{
  if (dict.found("gradient")) {
    gradient() = scalarField("gradient", dict, p.size());
    fixedGradientFvPatchScalarField::updateCoeffs();
    fixedGradientFvPatchScalarField::evaluate();
  } else {
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
  }
}


positiveGradientFvPatchScalarField::positiveGradientFvPatchScalarField
(
    const positiveGradientFvPatchScalarField& wbppsf
)
  :
fixedGradientFvPatchScalarField(wbppsf)

{}


positiveGradientFvPatchScalarField::positiveGradientFvPatchScalarField
(
    const positiveGradientFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
  :
fixedGradientFvPatchScalarField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void positiveGradientFvPatchScalarField::updateCoeffs()
{
  if (updated()) {
    return;
  }

  const fvPatchScalarField& Ziw = *this;

  gradient() = (Ziw - Ziw.patchInternalField()) * this->patch().deltaCoeffs() * pos((Ziw - Ziw.patchInternalField()));

  fixedGradientFvPatchScalarField::updateCoeffs();
}


void positiveGradientFvPatchScalarField::write(Ostream& os) const
{
  fvPatchScalarField::write(os);
  writeEntry(os, "gradient", gradient());
  writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, positiveGradientFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

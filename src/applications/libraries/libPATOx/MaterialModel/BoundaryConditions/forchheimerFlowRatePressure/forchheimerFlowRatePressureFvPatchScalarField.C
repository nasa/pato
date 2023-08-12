/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2011-2019 OpenFOAM Foundation
                            | Copyright (C) 2020      PATO
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

#include "forchheimerFlowRatePressureFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "one.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

forchheimerFlowRatePressureFvPatchScalarField::forchheimerFlowRatePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
  :
fixedGradientFvPatchScalarField(p, iF),
U_(0),
flowRate_(),
volumetric_(false),
rhoName_("rho_g"),
rhoInlet_(0.0)
{}


forchheimerFlowRatePressureFvPatchScalarField::forchheimerFlowRatePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
  :
fixedGradientFvPatchScalarField(p, iF),
U_(0),
rhoInlet_(dict.lookupOrDefault<scalar>("rhoInlet", -vGreat))
{
  if (dict.found("volumetricFlowRate")) {
    volumetric_ = true;
    flowRate_  = Function1<scalar>::New("volumetricFlowRate", dict);
    rhoName_ = "rho_g";
  } else if (dict.found("massFlowRate")) {
    volumetric_ = false;
    flowRate_ = Function1<scalar>::New("massFlowRate", dict);
    rhoName_ = word(dict.lookupOrDefault<word>("rho", "rho_g"));
  } else {
    FatalIOErrorInFunction(dict)
        << "Please supply either 'volumetricFlowRate' or"
        << " 'massFlowRate' and 'rho'" << exit(FatalIOError);
  }

  if (dict.found("gradient")) {
    gradient() = scalarField("gradient", dict, p.size());
    fixedGradientFvPatchScalarField::updateCoeffs();
    fixedGradientFvPatchScalarField::evaluate();
  } else {
    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = 0.0;
  }
}

forchheimerFlowRatePressureFvPatchScalarField::forchheimerFlowRatePressureFvPatchScalarField
(
    const forchheimerFlowRatePressureFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
  :
fixedGradientFvPatchScalarField(psf, p, iF, mapper),
U_(psf.U_),
flowRate_(psf.flowRate_, false),
volumetric_(psf.volumetric_),
rhoName_(psf.rhoName_),
rhoInlet_(psf.rhoInlet_)
{}


forchheimerFlowRatePressureFvPatchScalarField::forchheimerFlowRatePressureFvPatchScalarField
(
    const forchheimerFlowRatePressureFvPatchScalarField& psf
)
  :
fixedGradientFvPatchScalarField(psf),
U_(psf.U_),
flowRate_(psf.flowRate_, false),
volumetric_(psf.volumetric_),
rhoName_(psf.rhoName_),
rhoInlet_(psf.rhoInlet_)
{}


forchheimerFlowRatePressureFvPatchScalarField::forchheimerFlowRatePressureFvPatchScalarField
(
    const forchheimerFlowRatePressureFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
  :
fixedGradientFvPatchScalarField(psf, iF),
U_(psf.U_),
flowRate_(psf.flowRate_, false),
volumetric_(psf.volumetric_),
rhoName_(psf.rhoName_),
rhoInlet_(psf.rhoInlet_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class RhoType>
vectorField forchheimerFlowRatePressureFvPatchScalarField::updateValues
(
    const RhoType& rho
)
{
  const scalar t = db().time().timeOutputValue();
  const vectorField n(patch().nf());
  const scalar avgU = -flowRate_->value(t)/gSum(rho*patch().magSf());
  const vectorField flux(avgU*n);

  return flux;
}

void forchheimerFlowRatePressureFvPatchScalarField::updateCoeffs()
{
  if (updated()) {
    return;
  }

  if (volumetric_ || rhoName_ == "none") {
    U_ = updateValues(one());
  } else {
    if (db().foundObject<volScalarField>(rhoName_)) {
      const fvPatchField<scalar>& rhop =
          patch().lookupPatchField<volScalarField, scalar>(rhoName_);
      U_ = updateValues(rhop);
    } else {
      // use constant density
      if (rhoInlet_ < 0) {
        FatalErrorInFunction
            << "Did not find registered density field " << rhoName_
            << " and no constant density 'rhoInlet' specified"
            << exit(FatalError);
      }

      U_ = updateValues(rhoInlet_);
    }
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


void forchheimerFlowRatePressureFvPatchScalarField::write(Ostream& os) const
{
  fvPatchField<scalar>::write(os);
  writeEntry(os, flowRate_());
  if (!volumetric_) {
    writeEntryIfDifferent<word>(os, "rho", "rho_g", rhoName_);
    writeEntryIfDifferent<scalar>(os, "rhoInlet", -vGreat, rhoInlet_);
  }
  writeEntry(os, "gradient", gradient());
  writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, forchheimerFlowRatePressureFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

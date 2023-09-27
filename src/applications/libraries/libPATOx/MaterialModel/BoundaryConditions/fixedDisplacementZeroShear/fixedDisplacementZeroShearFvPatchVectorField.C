/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    fixedDisplacementZeroShearFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
  :
directionMixedFvPatchVectorField(p, iF),
traction_(p.size(), Zero),
pressure_(p.size(), 0.0)
{}


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
  :
directionMixedFvPatchVectorField(ptf, p, iF, mapper)
{}


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
  :
directionMixedFvPatchVectorField(p, iF),
traction_(p.size(), Zero),
pressure_(p.size(), 0.0)
{

  this->refGrad() = vector::zero;

  tmp<vectorField> n_tmp = patch().nf();
  vectorField n = n_tmp();
  this->valueFraction() = sqr(n);

  if (dict.found("value")) {
    Field<vector>::operator=(vectorField("value", dict, p.size()));
  } else {
    FatalError << "value entry not found for patch " << patch().name()
               << exit(FatalError);
  }

  this->refValue() = *this;

  tmp<Field<vector>> normalValue_tmp = transform(valueFraction(), refValue());
  Field<vector> normalValue = normalValue_tmp();

  tmp<Field<vector>> gradValue_tmp =
      this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();
  Field<vector> gradValue = gradValue_tmp();

  tmp<Field<vector>> transformGradValue_tmp =
      transform(I - valueFraction(), gradValue);
  Field<vector> transformGradValue = transformGradValue_tmp();

  Info << " patch    " << this->patch().name() <<  endl;
  Info << " VALUEEEE BC  " <<refValue() <<  endl;
  Info << " normalValue  " <<normalValue <<  endl;
  Info << " gradValue  " <<gradValue <<  endl;
  Info << " transformGradValue  " <<transformGradValue <<  endl;

  Field<vector>::operator=(normalValue + transformGradValue);
}


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
  :
directionMixedFvPatchVectorField(ptf, iF),
traction_(ptf.traction_),
pressure_(ptf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedDisplacementZeroShearFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
  directionMixedFvPatchVectorField::autoMap(m);
  m(traction_, traction_);
  m(pressure_, pressure_);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedDisplacementZeroShearFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
  directionMixedFvPatchVectorField::rmap(ptf, addr);

  const fixedDisplacementZeroShearFvPatchVectorField& dmptf =
      refCast<const fixedDisplacementZeroShearFvPatchVectorField>(ptf);

  traction_.rmap(dmptf.traction_, addr);
  pressure_.rmap(dmptf.pressure_, addr);
}


void fixedDisplacementZeroShearFvPatchVectorField::updateCoeffs()
{
  if (this->updated()) {
    return;
  }

  // Create result
  tmp<vectorField> tgradient(new vectorField(patch().size(), vector::zero));
  vectorField gradient = tgradient();

  // Standard isotropic solvers

  // Lookup material properties from the solver
  const fvPatchScalarField& mu =
      patch().lookupPatchField<volScalarField, scalar>("mu_sM");

  const fvPatchScalarField& lambda =
      patch().lookupPatchField<volScalarField, scalar>("lambda_sM");

  tmp<vectorField> n_tmp = patch().nf();
  vectorField n = n_tmp();

  // gradient of the field
  const fvPatchTensorField& gradField =
      patch().lookupPatchField<volTensorField, tensor>("gradD");


  // Calculate the normal gradient based on the traction

  gradient =
      (traction_ - n*pressure_)
      - (n & (mu*gradField.T() - (mu + lambda)*gradField))
      - n*lambda*tr(gradField);

  // Thermal effects

  const fvPatchField<scalar>& threeKalpha =
      patch().lookupPatchField<volScalarField, scalar>
      (
          "threeKalpha"
      );

  const fvPatchScalarField& T =
      patch().lookupPatchField<volScalarField, scalar>("Ta");

  const fvPatchScalarField& T_0 =
      patch().lookupPatchField<volScalarField, scalar>("T0");

  gradient += n*threeKalpha*(T-T_0);

  // Pyrolysis effects

  const fvPatchField<scalar>& threeKxi =
      patch().lookupPatchField<volScalarField, scalar>
      (
          "threeKxi"
      );

  const fvPatchScalarField& tau =
      patch().lookupPatchField<volScalarField, scalar>("tau");

  const fvPatchScalarField& tau_0 =
      patch().lookupPatchField<volScalarField, scalar>("tau0");

  gradient -= n*threeKxi*(tau_0-tau);

  gradient /= (2.0*mu + lambda);

  directionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void fixedDisplacementZeroShearFvPatchVectorField::write(Ostream& os) const
{
  directionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedDisplacementZeroShearFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

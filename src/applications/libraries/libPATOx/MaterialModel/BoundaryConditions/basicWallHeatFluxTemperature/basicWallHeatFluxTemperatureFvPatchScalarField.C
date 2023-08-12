/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "basicWallHeatFluxTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
template<>
const char*
NamedEnum
<
basicWallHeatFluxTemperatureFvPatchScalarField::operationMode,
                                               3
>::names[] = {
  "power",
  "flux",
  "coefficient"
};
}

const Foam::NamedEnum
<
Foam::basicWallHeatFluxTemperatureFvPatchScalarField::operationMode,
     3
     > Foam::basicWallHeatFluxTemperatureFvPatchScalarField::operationModeNames;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicWallHeatFluxTemperatureFvPatchScalarField::
basicWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
  :
mixedFvPatchScalarField(p, iF),
mode_(fixedHeatFlux),
Q_(0),
Tinf_(),
relaxation_(1),
emissivity_(0),
qrRelaxation_(1),
qrName_("undefined-qr"),
kappa_(1)
{
  refValue() = 0;
  refGrad() = 0;
  valueFraction() = 1;
}


Foam::basicWallHeatFluxTemperatureFvPatchScalarField::
basicWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
  :
mixedFvPatchScalarField(p, iF),
mode_(operationModeNames.read(dict.lookup("mode"))),
Q_(0),
Tinf_(),
relaxation_(dict.lookupOrDefault<scalar>("relaxation", 1)),
emissivity_(dict.lookupOrDefault<scalar>("emissivity", 0)),
qrRelaxation_(dict.lookupOrDefault<scalar>("qrRelaxation", 1)),
qrName_(dict.lookupOrDefault<word>("qr", "none")),
kappa_(dict.lookupOrDefault<scalar>("kappa", 1))
{
  switch (mode_) {
  case fixedPower: {
    dict.lookup("Q") >> Q_;

    break;
  }
  case fixedHeatFlux: {
    q_ = scalarField("q", dict, p.size());

    break;
  }
  case fixedHeatTransferCoeff: {
    h_ = scalarField("h", dict, p.size());
    Tinf_ = Function1<scalar>::New("Tinf", dict);

    break;
  }
  }

  fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

  if (qrName_ != "none") {
    if (dict.found("qrPrevious")) {
      qrPrevious_ = scalarField("qrPrevious", dict, p.size());
    } else {
      qrPrevious_.setSize(p.size(), 0);
    }
  }

  if (dict.found("refValue")) {
    // Full restart
    refValue() = scalarField("refValue", dict, p.size());
    refGrad() = scalarField("refGradient", dict, p.size());
    valueFraction() = scalarField("valueFraction", dict, p.size());
  } else {
    // Start from user entered data. Assume fixedValue.
    refValue() = *this;
    refGrad() = 0;
    valueFraction() = 1;
  }
}


Foam::basicWallHeatFluxTemperatureFvPatchScalarField::
basicWallHeatFluxTemperatureFvPatchScalarField
(
    const basicWallHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
  :
mixedFvPatchScalarField(ptf, p, iF, mapper),
mode_(ptf.mode_),
Q_(ptf.Q_),
Tinf_(ptf.Tinf_, false),
relaxation_(ptf.relaxation_),
emissivity_(ptf.emissivity_),
qrRelaxation_(ptf.qrRelaxation_),
qrName_(ptf.qrName_)
{
  switch (mode_) {
  case fixedPower: {
    break;
  }
  case fixedHeatFlux: {
    mapper(q_, ptf.q_);
    break;
  }
  case fixedHeatTransferCoeff: {
    mapper(h_, ptf.h_);
    break;
  }
  }

  if (qrName_ != "none") {
    mapper(qrPrevious_, ptf.qrPrevious_);
  }
}


Foam::basicWallHeatFluxTemperatureFvPatchScalarField::
basicWallHeatFluxTemperatureFvPatchScalarField
(
    const basicWallHeatFluxTemperatureFvPatchScalarField& tppsf
)
  :
mixedFvPatchScalarField(tppsf),
mode_(tppsf.mode_),
Q_(tppsf.Q_),
q_(tppsf.q_),
h_(tppsf.h_),
Tinf_(tppsf.Tinf_, false),
relaxation_(tppsf.relaxation_),
emissivity_(tppsf.emissivity_),
qrPrevious_(tppsf.qrPrevious_),
qrRelaxation_(tppsf.qrRelaxation_),
qrName_(tppsf.qrName_)
{}


Foam::basicWallHeatFluxTemperatureFvPatchScalarField::
basicWallHeatFluxTemperatureFvPatchScalarField
(
    const basicWallHeatFluxTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
  :
mixedFvPatchScalarField(tppsf, iF),
mode_(tppsf.mode_),
Q_(tppsf.Q_),
q_(tppsf.q_),
h_(tppsf.h_),
Tinf_(tppsf.Tinf_, false),
relaxation_(tppsf.relaxation_),
emissivity_(tppsf.emissivity_),
qrPrevious_(tppsf.qrPrevious_),
qrRelaxation_(tppsf.qrRelaxation_),
qrName_(tppsf.qrName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::basicWallHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
  mixedFvPatchScalarField::autoMap(m);

  switch (mode_) {
  case fixedPower: {
    break;
  }
  case fixedHeatFlux: {
    m(q_, q_);

    break;
  }
  case fixedHeatTransferCoeff: {
    m(h_, h_);

    break;
  }
  }

  if (qrName_ != "none") {
    m(qrPrevious_, qrPrevious_);
  }
}


void Foam::basicWallHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
  mixedFvPatchScalarField::rmap(ptf, addr);

  const basicWallHeatFluxTemperatureFvPatchScalarField& tiptf =
      refCast<const basicWallHeatFluxTemperatureFvPatchScalarField>(ptf);

  switch (mode_) {
  case fixedPower: {
    break;
  }
  case fixedHeatFlux: {
    q_.rmap(tiptf.q_, addr);

    break;
  }
  case fixedHeatTransferCoeff: {
    h_.rmap(tiptf.h_, addr);

    break;
  }
  }

  if (qrName_ != "none") {
    qrPrevious_.rmap(tiptf.qrPrevious_, addr);
  }
}


void Foam::basicWallHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
  if (updated()) {
    return;
  }

  const scalarField& Tp(*this);

  // Store current valueFraction and refValue for relaxation
  const scalarField valueFraction0(valueFraction());
  const scalarField refValue0(refValue());

  scalarField qr(Tp.size(), 0);
  if (qrName_ != "none") {
    qr =
        qrRelaxation_
        *patch().lookupPatchField<volScalarField, scalar>(qrName_)
        + (1 - qrRelaxation_)*qrPrevious_;

    qrPrevious_ = qr;
  }

  switch (mode_) {
  case fixedPower: {
    refGrad() = (Q_/gSum(patch().magSf()) + qr)/kappa_;
    refValue() = Tp;
    valueFraction() = 0;

    break;
  }
  case fixedHeatFlux: {
    refGrad() = (q_ + qr)/kappa_;
    refValue() = Tp;
    valueFraction() = 0;

    break;
  }
  case fixedHeatTransferCoeff: {
    scalarField hp(h_);

    const scalar Tinf = Tinf_->value(this->db().time().timeOutputValue());
    scalarField hpTinf(hp*Tinf);

    if (emissivity_ > 0) {
      // Evaluate the radiative flux to the environment
      // from the surface temperature ...
      hp += emissivity_*sigma.value()*pow3(Tp);
      hpTinf += emissivity_*sigma.value()*pow4(Tinf);
    }

    const scalarField kappaDeltaCoeffs
    (
        this->kappa_*patch().deltaCoeffs()
    );

    refGrad() = 0;

    forAll(Tp, i) {
      if (qr[i] < 0) {
        const scalar hpmqr = hp[i] - qr[i]/Tp[i];

        refValue()[i] = hpTinf[i]/hpmqr;
        valueFraction()[i] = hpmqr/(hpmqr + kappaDeltaCoeffs[i]);
      } else {
        refValue()[i] = (hpTinf[i] + qr[i])/hp[i];
        valueFraction()[i] = hp[i]/(hp[i] + kappaDeltaCoeffs[i]);
      }
    }

    break;
  }
  }

  valueFraction() =
      relaxation_*valueFraction()
      + (1 - relaxation_)*valueFraction0;

  refValue() = relaxation_*refValue() + (1 - relaxation_)*refValue0;

  mixedFvPatchScalarField::updateCoeffs();

  if (debug) {
    const scalar Q = gSum(kappa_*patch().magSf()*snGrad());

    Info<< patch().boundaryMesh().mesh().name() << ':'
        << patch().name() << ':'
        << this->internalField().name() << " :"
        << " heat transfer rate:" << Q
        << " walltemperature "
        << " min:" << gMin(*this)
        << " max:" << gMax(*this)
        << " avg:" << gAverage(*this)
        << endl;
  }
}


void Foam::basicWallHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
  fvPatchScalarField::write(os);

  writeEntry(os, "mode", operationModeNames[mode_]);

  switch (mode_) {
  case fixedPower: {
    writeEntry(os, "Q", Q_);

    break;
  }
  case fixedHeatFlux: {
    writeEntry(os, "q", q_);

    break;
  }
  case fixedHeatTransferCoeff: {
    writeEntry(os, "h", h_);
    writeEntry(os, Tinf_());

    if (relaxation_ < 1) {
      writeEntry(os, "relaxation", relaxation_);
    }

    if (emissivity_ > 0) {
      writeEntry(os, "emissivity", emissivity_);
    }
    break;
  }
  }

  writeEntry(os, "qr", qrName_);

  if (qrName_ != "none") {
    writeEntry(os, "qrRelaxation", qrRelaxation_);

    writeEntry(os, "qrPrevious", qrPrevious_);
  }

  writeEntry(os, "refValue", refValue());
  writeEntry(os, "refGradient", refGrad());
  writeEntry(os, "valueFraction", valueFraction());
  writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    basicWallHeatFluxTemperatureFvPatchScalarField
);
}

// ************************************************************************* //

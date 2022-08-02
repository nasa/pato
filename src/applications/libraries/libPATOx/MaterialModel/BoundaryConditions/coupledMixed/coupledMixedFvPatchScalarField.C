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

#include "coupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
template<>
const char* Foam::NamedEnum
<
Foam::coupledMixedFvPatchScalarField::KMethodType,
     4
>::names[] = {
  "fluidThermo",
  "solidThermo",
  "directionalSolidThermo",
  "lookup"
};
}

const Foam::NamedEnum<Foam::coupledMixedFvPatchScalarField::KMethodType, 4> Foam::coupledMixedFvPatchScalarField::KMethodTypeNames_;

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledMixedFvPatchScalarField::
coupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
  :
mixedFvPatchScalarField(p, iF),
patch_(patch()),
method_(KMethodTypeNames_["undefined"]),
kappaName_("undefined"),
alphaAniName_("undefined-K"),
TnbrName_("undefined-Tnbr"),
thicknessLayers_(0),
kappaLayers_(0),
contactRes_(0)
{
  this->refValue() = 0.0;
  this->refGrad() = 0.0;
  this->valueFraction() = 1.0;
}


coupledMixedFvPatchScalarField::
coupledMixedFvPatchScalarField
(
    const coupledMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
  :
mixedFvPatchScalarField(ptf, p, iF, mapper),
patch_(patch()),
method_(ptf.method_),
kappaName_(ptf.kappaName_),
alphaAniName_(ptf.alphaAniName_),
TnbrName_(ptf.TnbrName_),
thicknessLayers_(ptf.thicknessLayers_),
kappaLayers_(ptf.kappaLayers_),
contactRes_(ptf.contactRes_)
{}


coupledMixedFvPatchScalarField::
coupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
  :
mixedFvPatchScalarField(p, iF),
patch_(patch()),
method_(KMethodTypeNames_.read(dict.lookup("kappaMethod"))),
kappaName_(dict.lookupOrDefault<word>("kappa", "none")),
alphaAniName_(dict.lookupOrDefault<word>("alphaAni","Anialpha")),
TnbrName_(dict.lookup("Tnbr")),
thicknessLayers_(0),
kappaLayers_(0),
contactRes_(0.0)
{
  switch (method_) {
  case mtDirectionalSolidThermo: {
    if (!dict.found("alphaAni")) {
      FatalIOErrorInFunction(dict)
          << "Did not find entry 'alphaAni'"
          " required for 'kappaMethod' "
          << KMethodTypeNames_[method_]
          << exit(FatalIOError);
    }

    break;
  }

  case mtLookup: {
    if (!dict.found("kappa")) {
      FatalIOErrorInFunction(dict)
          << "Did not find entry 'kappa'"
          " required for 'kappaMethod' "
          <<  KMethodTypeNames_[method_] << nl
          << "    Please set 'kappa' to the name of a volScalarField"
          " or volSymmTensorField"
          << exit(FatalIOError);
    }

    break;
  }

  default:
    break;
  }

  if (!isA<mappedPatchBase>(this->patch().patch())) {
    FatalErrorInFunction
        << "' not type '" << mappedPatchBase::typeName << "'"
        << "\n    for patch " << p.name()
        << " of field " << internalField().name()
        << " in file " << internalField().objectPath()
        << exit(FatalError);
  }

  if (dict.found("thicknessLayers")) {
    dict.lookup("thicknessLayers") >> thicknessLayers_;
    dict.lookup("kappaLayers") >> kappaLayers_;

    if (thicknessLayers_.size() > 0) {
      // Calculate effective thermal resistance by harmonic averaging
      forAll(thicknessLayers_, iLayer) {
        contactRes_ += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
      }
      contactRes_ = 1.0/contactRes_;
    }
  }

  fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

  if (dict.found("refValue")) {
    // Full restart
    refValue() = scalarField("refValue", dict, p.size());
    refGrad() = scalarField("refGradient", dict, p.size());
    valueFraction() = scalarField("valueFraction", dict, p.size());
  } else {
    // Start from user entered data. Assume fixedValue.
    refValue() = *this;
    refGrad() = 0.0;
    valueFraction() = 1.0;
  }
}


coupledMixedFvPatchScalarField::
coupledMixedFvPatchScalarField
(
    const coupledMixedFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
  :
mixedFvPatchScalarField(wtcsf, iF),
patch_(patch()),
method_(wtcsf.method_),
kappaName_(wtcsf.kappaName_),
alphaAniName_(wtcsf.alphaAniName_),
TnbrName_(wtcsf.TnbrName_),
thicknessLayers_(wtcsf.thicknessLayers_),
kappaLayers_(wtcsf.kappaLayers_),
contactRes_(wtcsf.contactRes_)
{}

coupledMixedFvPatchScalarField::~coupledMixedFvPatchScalarField()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coupledMixedFvPatchScalarField::updateCoeffs()
{
  if (updated()) {
    return;
  }

  // Since we're inside initEvaluate/evaluate there might be processor
  // comms underway. Change the tag we use.
  int oldTag = UPstream::msgType();
  UPstream::msgType() = oldTag+1;

  // Get the coupling information from the mappedPatchBase
  const mappedPatchBase& mpp =
      refCast<const mappedPatchBase>(patch().patch());
  const polyMesh& nbrMesh = mpp.sampleMesh();
  const label samplePatchi = mpp.samplePolyPatch().index();
  const fvPatch& nbrPatch =
      refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

  // Calculate the temperature by harmonic averaging
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  const coupledMixedFvPatchScalarField& nbrField =
      refCast
      <
      const coupledMixedFvPatchScalarField
      >
      (
          nbrPatch.lookupPatchField<volScalarField, scalar>
          (
              TnbrName_
          )
      );

  // Swap to obtain full local values of neighbour internal field
  tmp<scalarField> nbrIntFld(new scalarField(nbrField.size(), 0.0));
  tmp<scalarField> nbrKDelta(new scalarField(nbrField.size(), 0.0));

  if (contactRes_ == 0.0) {
    nbrIntFld.ref() = nbrField.patchInternalField();
    nbrKDelta.ref() = nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs();
  } else {
    nbrIntFld.ref() = nbrField;
    nbrKDelta.ref() = contactRes_;
  }

  mpp.distribute(nbrIntFld.ref());
  mpp.distribute(nbrKDelta.ref());

  tmp<scalarField> myKDelta = kappa(*this)*patch().deltaCoeffs();

  // Both sides agree on
  // - temperature : (myKDelta*fld + nbrKDelta*nbrFld)/(myKDelta+nbrKDelta)
  // - gradient    : (temperature-fld)*delta
  // We've got a degree of freedom in how to implement this in a mixed bc.
  // (what gradient, what fixedValue and mixing coefficient)
  // Two reasonable choices:
  // 1. specify above temperature on one side (preferentially the high side)
  //    and above gradient on the other. So this will switch between pure
  //    fixedvalue and pure fixedgradient
  // 2. specify gradient and temperature such that the equations are the
  //    same on both sides. This leads to the choice of
  //    - refGradient = zero gradient
  //    - refValue = neighbour value
  //    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

  this->refValue() = nbrIntFld();
  this->refGrad() = 0.0;
  this->valueFraction() = nbrKDelta()/(nbrKDelta() + myKDelta());

  mixedFvPatchScalarField::updateCoeffs();

  if (debug) {
    scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

    Info<< patch().boundaryMesh().mesh().name() << ':'
        << patch().name() << ':'
        << this->internalField().name() << " <- "
        << nbrMesh.name() << ':'
        << nbrPatch.name() << ':'
        << this->internalField().name() << " :"
        << " heat transfer rate:" << Q
        << " walltemperature "
        << " min:" << gMin(*this)
        << " max:" << gMax(*this)
        << " avg:" << gAverage(*this)
        << endl;
  }

  // Restore tag
  UPstream::msgType() = oldTag;

// mixedFvPatchField.C
//    Field<Type>::operator=
//    (
//        valueFraction_*refValue_
//      +
//        (1.0 - valueFraction_)*
//        (
//            this->patchInternalField()
//          + refGrad_/this->patch().deltaCoeffs()
//        )
//    );

}


void coupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
  mixedFvPatchScalarField::write(os);
  os.writeKeyword("Tnbr")<< TnbrName_
                         << token::END_STATEMENT << nl;
  writeEntry(os, "thicknessLayers", thicknessLayers_);
  writeEntry(os, "kappaLayers", kappaLayers_);

  os.writeKeyword("kappaMethod") << KMethodTypeNames_[method_]
                                 << token::END_STATEMENT << nl;
  os.writeKeyword("kappa") << kappaName_ << token::END_STATEMENT << nl;
  os.writeKeyword("alphaAni") << alphaAniName_ << token::END_STATEMENT << nl;
}

tmp<scalarField> coupledMixedFvPatchScalarField::kappa
(
    const scalarField& Tp
) const
{
  const fvMesh& mesh = patch_.boundaryMesh().mesh();
  const label patchi = patch_.index();

  switch (method_) {
  case mtFluidThermo: {
    typedef compressible::turbulenceModel turbulenceModel;

    word turbName(turbulenceModel::propertiesName);

    if
    (
        mesh.foundObject<turbulenceModel>(turbName)
    ) {
      const turbulenceModel& turbModel =
          mesh.lookupObject<turbulenceModel>(turbName);

      return turbModel.kappaEff(patchi);
    } else if (mesh.foundObject<fluidThermo>(basicThermo::dictName)) {
      const fluidThermo& thermo =
          mesh.lookupObject<fluidThermo>(basicThermo::dictName);

      return thermo.kappa(patchi);
    } else {
      FatalErrorInFunction
          << "kappaMethod defined to employ "
          << KMethodTypeNames_[method_]
          << " method, but thermo package not available"
          << exit(FatalError);
    }

    break;
  }

  case mtSolidThermo: {
    const solidThermo& thermo =
        mesh.lookupObject<solidThermo>(basicThermo::dictName);

    return thermo.kappa(patchi);
    break;
  }

  case mtDirectionalSolidThermo: {
    const solidThermo& thermo =
        mesh.lookupObject<solidThermo>(basicThermo::dictName);

    const symmTensorField& alphaAni =
        patch_.lookupPatchField<volSymmTensorField, scalar>
        (
            alphaAniName_
        );

    const scalarField& pp = thermo.p().boundaryField()[patchi];

    const symmTensorField kappa(alphaAni*thermo.Cp(pp, Tp, patchi));

    const vectorField n(patch_.nf());

    return n & kappa & n;
  }

  case mtLookup: {
    if (mesh.foundObject<volScalarField>(kappaName_)) {
      return patch_.lookupPatchField<volScalarField, scalar>
             (
                 kappaName_
             );
    } else if (mesh.foundObject<volSymmTensorField>(kappaName_)) {
      const symmTensorField& KWall =
          patch_.lookupPatchField<volSymmTensorField, scalar>
          (
              kappaName_
          );

      const vectorField n(patch_.nf());

      return n & KWall & n;
    } else {
      FatalErrorInFunction
          << "Did not find field " << kappaName_
          << " on mesh " << mesh.name() << " patch " << patch_.name()
          << nl
          << "    Please set 'kappa' to the name of a volScalarField"
          " or volSymmTensorField."
          << exit(FatalError);
    }

    break;
  }

  default: {
    FatalErrorInFunction
        << "Unimplemented method " << KMethodTypeNames_[method_] << nl
        << "    Please set 'kappaMethod' to one of "
        << KMethodTypeNames_.toc()
        << " and 'kappa' to the name of the volScalar"
        << " or volSymmTensor field (if kappa=lookup)"
        << exit(FatalError);
  }
  }

  return scalarField(0);
}


void coupledMixedFvPatchScalarField::setTnbrName(word name)
{
  TnbrName_=name;
}

void coupledMixedFvPatchScalarField::setKappa(word name)
{
  kappaName_=name;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    coupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //

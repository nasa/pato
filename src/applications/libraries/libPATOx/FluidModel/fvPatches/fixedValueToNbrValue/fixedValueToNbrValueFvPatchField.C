/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "fixedValueToNbrValueFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class Type>
Foam::fixedValueToNbrValueFvPatchField<Type>::fixedValueToNbrValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
  :
mixedFvPatchField<Type>(p, iF),
nbrName_("undefined-nbr")
{
  this->refValue() = pTraits<Type>::zero; // is loaded with <type>, vector
  this->refGrad() = pTraits<Type>::zero; // is loaded with <type>, vector
  this->valueFraction() = 0.0; // is loaded with scalar, scalar
}


template<class Type>
Foam::fixedValueToNbrValueFvPatchField<Type>::fixedValueToNbrValueFvPatchField
(
    const fixedValueToNbrValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
  :
mixedFvPatchField<Type>(ptf, p, iF, mapper),
nbrName_(ptf.nbrName_)
{}


template<class Type>
Foam::fixedValueToNbrValueFvPatchField<Type>::fixedValueToNbrValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
  :
mixedFvPatchField<Type>(p, iF),
nbrName_(dict.lookup("nbr"))
{
  /*
     if (!isA<mappedPatchBase>(this->patch().patch()))
     {
         FatalErrorIn
         (
             "fixedValueToNbrValueFvPatchField::"
             "fixedValueToNbrValueFvPatchField\n"
             "(\n"
             "    const fvPatch& p,\n"
             "    const DimensionedField<Type, volMesh>& iF,\n"
             "    const dictionary& dict\n"
             ")\n"
         )   << "\n    patch type '" << p.type()
             << "' not type '" << mappedPatchBase::typeName << "'"
             << "\n    for patch " << p.name()
             << " of field " << dimensionedInternalField().name()
             << " in file " << dimensionedInternalField().objectPath()
             << exit(FatalError);
     }*/

  if (dict.found("value")) {
    fvPatchField<Type>::operator=(Field<Type>("value", dict, p.size()));
  } else {
    fvPatchField<Type>::operator=(this->patchInternalField());
  }

  if (dict.found("refGradient")) {
    this->refGrad() = Field<Type>("refGradient", dict, p.size()); // in parent is <type>
  } else {
    this->refGrad() = pTraits<Type>::zero; // in parent is <type>
  }

  if (dict.found("valueFraction")) {
    this->valueFraction() = scalarField("valueFraction", dict, p.size()); // in parent is scalar
  } else {
    this->valueFraction() = 1.0; // in parent is <type> // 0.0 fixes to own value, 1.0 to nbr value
  }

  //this->mu() = read mu from porous medium to be used in fluid region;

}


template<class Type>
Foam::fixedValueToNbrValueFvPatchField<Type>::fixedValueToNbrValueFvPatchField
(
    const fixedValueToNbrValueFvPatchField<Type>& wtcsf,
    const DimensionedField<Type, volMesh>& iF
)
  :
mixedFvPatchField<Type>(wtcsf, iF),
nbrName_(wtcsf.nbrName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type>
void Foam::fixedValueToNbrValueFvPatchField<Type>::updateCoeffs()
{
  if (this->updated()) {
    return;
  }

  // Since we're inside initEvaluate/evaluate there might be processor
  // comms underway. Change the tag we use.
  int oldTag = UPstream::msgType();
  UPstream::msgType() = oldTag+1;

  // Get the coupling information from the mappedPatchBase
  const mappedPatchBase& mpp =  refCast<const mappedPatchBase>(this->patch().patch());

  const polyMesh& nbrMesh = mpp.sampleMesh();

  const label samplePatchI = mpp.samplePolyPatch().index();
  const fvPatch& nbrPatch = refCast<const fvMesh>(nbrMesh).boundary()[samplePatchI];

  const fvPatchField<Type>& nbrField =
      nbrPatch.lookupPatchField<GeometricField<Type, fvPatchField, volMesh>, Type>(nbrName_);

  /*    const vectorCoupledMixedFvPatchVectorField& nbrField =
          refCast<const vectorCoupledMixedFvPatchVectorField>
          (nbrPatch.lookupPatchField<volVectorField, vector>(nbrName_));  // this segment makes it impossible to have different conditions on different boundaries */

  // Swap to obtain full local values of neighbour internal field
  tmp<Field<Type>> nbrIntFld
                (
                    new Field<Type>(nbrField.size(), pTraits<Type>::zero)
                );

  nbrIntFld.ref() = nbrField.patchInternalField();
  mpp.distribute(nbrIntFld.ref());


  // Both sides agree on
  // - pressure : (myKDelta*fld + nbrKDelta*nbrFld)/(myKDelta+nbrKDelta)
  // - gradient    : (pressure-fld)*delta
  // We've got a degree of freedom in how to implement this in a mixed bc.
  // (what gradient, what fixedValue and mixing coefficient)
  // Two reasonable choices:
  // 1. specify above pressure on one side (preferentially the high side)
  //    and above gradient on the other. So this will switch between pure
  //    fixedvalue and pure fixedgradient
  // 2. specify gradient and pressure such that the equations are the
  //    same on both sides. This leads to the choice of
  //    - refGradient = zero gradient
  //    - refValue = neighbour value
  //    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())


  this->refValue() = nbrIntFld(); // multiply by mu from porous medium
  this->refGrad() = pTraits<Type>::zero;
  this->valueFraction() = 1.0;
// instead of changing BC, I can probably read muv in and define muv as mu*v inside of the solver and set it to write. However, the units dont appear to make sense anymore if I do this. We can also substitute kappa for mu_porous and mu_fluid, since this ratio matters, see DIN-A4

  mixedFvPatchField<Type>::updateCoeffs();
  /*
      if (debug)
      {
          scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

          Info<< patch().boundaryMesh().mesh().name() << ':'
              << patch().name() << ':'
              << this->dimensionedInternalField().name() << " <- "
              << nbrMesh.name() << ':'
              << nbrPatch.name() << ':'
              << this->dimensionedInternalField().name() << " :"
              << " heat transfer rate:" << Q
              << " wallpressure "
              << " min:" << gMin(*this)
              << " max:" << gMax(*this)
              << " avg:" << gAverage(*this)
              << endl;
      }
  */
  // Restore tag
  UPstream::msgType() = oldTag;
}

template<class Type>
void Foam::fixedValueToNbrValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
  mixedFvPatchField<Type>::write(os);
  os.writeKeyword("nbr")<< nbrName_
                        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
template<class Type>
Foam::makePatchTypeField // what does this do?
(
    fvPatchField<Type>,
    fixedValueToNbrValueFvPatchField<Type> // write a typedeff or so? how can I solve this?
);
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //

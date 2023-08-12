/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "boundaryMappingFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryMappingFvPatchScalarField::boundaryMappingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
  :
fixedValueFvPatchScalarField(p, iF),
mesh_(patch().boundaryMesh().mesh()),
dict_(initDict()),
patchName_(word::null),
boundaryMapping_(simpleBoundaryMappingModel::New(
    mesh_,
    neededFields_,
    dict_)
                ),
boundaryMapping_ptr(&boundaryMapping_())
{
  FatalErrorInFunction << "boundaryMappingFvPatchScalarField(p,iF) not implemented." << exit(FatalError);
}


Foam::boundaryMappingFvPatchScalarField::boundaryMappingFvPatchScalarField
(
    const boundaryMappingFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
  :
fixedValueFvPatchScalarField(ptf, p, iF, mapper),
mesh_(ptf.mesh_),
dict_(ptf.dict_),
mappingFields_(ptf.mappingFields_),
patchName_(ptf.patchName_),
boundaryMapping_(ptf.boundaryMapping_),
boundaryMapping_ptr(ptf.boundaryMapping_ptr)
{
}


Foam::boundaryMappingFvPatchScalarField::boundaryMappingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
  :
fixedValueFvPatchScalarField(p, iF, dict),
mesh_(patch().boundaryMesh().mesh()),
dict_(dict),
mappingFields_(dict_.lookup("mappingFields")),
patchName_(initPatchName()),
boundaryMapping_(
    simpleBoundaryMappingModel::New(
        mesh_,
        neededFields_,
        dict_)
),
boundaryMapping_ptr(&boundaryMapping_())
{
}


Foam::boundaryMappingFvPatchScalarField::boundaryMappingFvPatchScalarField
(
    const boundaryMappingFvPatchScalarField& frpsf
)
  :
fixedValueFvPatchScalarField(frpsf),
mesh_(frpsf.mesh_),
dict_(frpsf.dict_),
mappingFields_(frpsf.mappingFields_),
patchName_(frpsf.patchName_),
boundaryMapping_(frpsf.boundaryMapping_),
boundaryMapping_ptr(frpsf.boundaryMapping_ptr)
{
}


Foam::boundaryMappingFvPatchScalarField::boundaryMappingFvPatchScalarField
(
    const boundaryMappingFvPatchScalarField& frpsf,
    const DimensionedField<scalar, volMesh>& iF
)
  :
fixedValueFvPatchScalarField(frpsf, iF),
mesh_(frpsf.mesh_),
dict_(frpsf.dict_),
mappingFields_(frpsf.mappingFields_),
patchName_(frpsf.patchName_),
boundaryMapping_(frpsf.boundaryMapping_),
boundaryMapping_ptr(frpsf.boundaryMapping_ptr)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::boundaryMappingFvPatchScalarField::initPatchName()
{
  if(mappingFields_.size()!=1) {
    FatalErrorInFunction << "boundaryMapping Boundary Conditions: \"mappingFields\" must have a size of 1." << exit(FatalError);
  }
  word patchName_ =mappingFields_[0].first();
  wordList names_;
  names_.append(patchName_);
  foundFieldsInMesh(mesh_,names_);
  return patchName_;
}

Foam::dictionary Foam::boundaryMappingFvPatchScalarField::initDict()
{
  dictionary dict_;
  dict_.add("mappingType","constant");
  fileName file_ = "$FOAM_CASE/"+db().time().constant();
  dict_.add("mappingFileName",file_);
  List<Tuple2<word,scalar>> tuple_;
  tuple_.append(Tuple2<word,scalar>("p",1));
  dict_.add("mappingFields",tuple_);
  return dict_;
}

void Foam::boundaryMappingFvPatchScalarField::updateCoeffs()
{
  if (updated()) {
    return;
  }
//  scalarField& pw = *this;
  const label patchID = patch().index();
  const fvMesh& mesh_  =  patch().boundaryMesh().mesh();
  const Time& runTime_ = mesh_.time();
  boundaryMapping_ptr->update(runTime_.value(),patchID,patchName_);
  fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::boundaryMappingFvPatchScalarField::write(Ostream& os) const
{
  fvPatchScalarField::write(os);
  boundaryMapping_ptr->write(os);
  writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    boundaryMappingFvPatchScalarField
);
}

// ************************************************************************* //

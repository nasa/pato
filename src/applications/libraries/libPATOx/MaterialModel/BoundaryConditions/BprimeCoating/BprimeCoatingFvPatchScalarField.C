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

#include "BprimeCoatingFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BprimeCoatingFvPatchScalarField::BprimeCoatingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
  :
fixedValueFvPatchScalarField(p, iF),
debug_("no"),
mesh(patch().boundaryMesh().mesh()),
phaseName(word::null),
dictName(mesh.name()),
dict_(initDict()),
BprimeCoatingBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict_
)
{
  FatalErrorInFunction << "BprimeCoatingFvPatchScalarField(p,iF) not implemented." << exit(FatalError);
}


Foam::BprimeCoatingFvPatchScalarField::BprimeCoatingFvPatchScalarField
(
    const BprimeCoatingFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
  :
fixedValueFvPatchScalarField(ptf, p, iF, mapper),
debug_(ptf.debug_),
mesh(ptf.mesh),
phaseName(ptf.phaseName),
dictName(ptf.dictName),
dict_(ptf.dict_),
BprimeCoatingBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict_
)
{
}


Foam::BprimeCoatingFvPatchScalarField::BprimeCoatingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
  :
fixedValueFvPatchScalarField(p, iF, dict),
debug_(dict.lookupOrDefault<Switch>("debug","no")),
mesh(patch().boundaryMesh().mesh()),
phaseName(word::null),
dictName(mesh.name()),
BprimeCoatingBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict
)
{
}


Foam::BprimeCoatingFvPatchScalarField::BprimeCoatingFvPatchScalarField
(
    const BprimeCoatingFvPatchScalarField& frpsf
)
  :
fixedValueFvPatchScalarField(frpsf),
debug_(frpsf.debug_),
mesh(frpsf.mesh),
phaseName(frpsf.phaseName),
dictName(frpsf.dictName),
dict_(frpsf.dict_),
BprimeCoatingBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict_
)
{
}


Foam::BprimeCoatingFvPatchScalarField::BprimeCoatingFvPatchScalarField
(
    const BprimeCoatingFvPatchScalarField& frpsf,
    const DimensionedField<scalar, volMesh>& iF
)
  :
fixedValueFvPatchScalarField(frpsf, iF),
debug_(frpsf.debug_),
mesh(frpsf.mesh),
phaseName(frpsf.phaseName),
dictName(frpsf.dictName),
dict_(frpsf.dict_),
BprimeCoatingBoundaryConditions_(frpsf.BprimeCoatingBoundaryConditions_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::dictionary Foam::BprimeCoatingFvPatchScalarField::initDict()
{
  dictionary dict_;
  dict_.add("mixtureMutationBprime","tacot26");
  fileName envi_("$PATO_DIR/data/Environments/RawData/Earth");
  dict_.add("environmentDirectory",envi_);
  dict_.add("movingMesh","yes");
  dict_.add("mappingType","constant");
  fileName mapp_("constant/porousMat/BoundaryConditions");
  dict_.add("mappingFileName",mapp_);
  List<Tuple2<word,scalar>> tuple_;
  tuple_.append(Tuple2<word,scalar>("p",1));
  tuple_.append(Tuple2<word,scalar>("rhoeUeCH",4));
  tuple_.append(Tuple2<word,scalar>("h_r",5));
  tuple_.append(Tuple2<word,scalar>("qRad",6));
  tuple_.append(Tuple2<word,scalar>("heatOn",8));
  tuple_.append(Tuple2<word,scalar>("lambda",9));
  tuple_.append(Tuple2<word,scalar>("Tbackground",10));
  dict_.add("mappingFields",tuple_);
  return dict_;
}

void Foam::BprimeCoatingFvPatchScalarField::updateCoeffs()
{
  if (updated()) {
    return;
  }

//  scalarField& Tw = *this;
  if(debug_) {
    Info << "--- BprimeCoatingBoundaryConditions_.update(); --- Foam::BprimeCoatingFvPatchScalarField::updateCoeffs()" << endl;
  }
  BprimeCoatingBoundaryConditions_.update();

  if(debug_) {
    Info << "--- fixedValueFvPatchScalarField::updateCoeffs(); --- Foam::BprimeCoatingFvPatchScalarField::updateCoeffs()" << endl;
  }
  fixedValueFvPatchScalarField::updateCoeffs();
  if(debug_) {
    Info << "--- end --- Foam::BprimeCoatingFvPatchScalarField::updateCoeffs()" << endl;
  }
}


void Foam::BprimeCoatingFvPatchScalarField::write(Ostream& os) const
{
  fvPatchScalarField::write(os);
  BprimeCoatingBoundaryConditions_.write(os);
  writeEntry(os, "value", *this);
  //    writeEntryIfDifferent<scalar>(os, "inclinationAngle", -VGREAT, inclinationAngle_);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    BprimeCoatingFvPatchScalarField
);
}

// ************************************************************************* //

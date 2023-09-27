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

#include "qRad_emission_absorptionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::qRad_emission_absorptionFvPatchScalarField::qRad_emission_absorptionFvPatchScalarField
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
qRad_emission_absorptionBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict_
)
{
  FatalErrorInFunction << "qRad_emission_absorptionFvPatchScalarField(p,iF) not implemented." << exit(FatalError);
}


Foam::qRad_emission_absorptionFvPatchScalarField::qRad_emission_absorptionFvPatchScalarField
(
    const qRad_emission_absorptionFvPatchScalarField& ptf,
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
qRad_emission_absorptionBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict_
)
{
}


Foam::qRad_emission_absorptionFvPatchScalarField::qRad_emission_absorptionFvPatchScalarField
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
qRad_emission_absorptionBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict
)
{
}


Foam::qRad_emission_absorptionFvPatchScalarField::qRad_emission_absorptionFvPatchScalarField
(
    const qRad_emission_absorptionFvPatchScalarField& frpsf
)
  :
fixedValueFvPatchScalarField(frpsf),
debug_(frpsf.debug_),
mesh(frpsf.mesh),
phaseName(frpsf.phaseName),
dictName(frpsf.dictName),
dict_(frpsf.dict_),
qRad_emission_absorptionBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict_
)
{
}


Foam::qRad_emission_absorptionFvPatchScalarField::qRad_emission_absorptionFvPatchScalarField
(
    const qRad_emission_absorptionFvPatchScalarField& frpsf,
    const DimensionedField<scalar, volMesh>& iF
)
  :
fixedValueFvPatchScalarField(frpsf, iF),
debug_(frpsf.debug_),
mesh(frpsf.mesh),
phaseName(frpsf.phaseName),
dictName(frpsf.dictName),
dict_(frpsf.dict_),
qRad_emission_absorptionBoundaryConditions_(frpsf.qRad_emission_absorptionBoundaryConditions_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::dictionary Foam::qRad_emission_absorptionFvPatchScalarField::initDict()
{
  dictionary dict_;
  dict_.add("mappingType","constant");
  fileName mapp_("constant/porousMat/BoundaryConditions");
  dict_.add("mappingFileName",mapp_);
  List<Tuple2<word,scalar>> tuple_;
  tuple_.append(Tuple2<word,scalar>("p",1));
  tuple_.append(Tuple2<word,scalar>("rhoeUeCH",4));
  tuple_.append(Tuple2<word,scalar>("h_r",5));
  tuple_.append(Tuple2<word,scalar>("qRad",6));
  tuple_.append(Tuple2<word,scalar>("Tbackground",10));
  dict_.add("mappingFields",tuple_);
  return dict_;
}

void Foam::qRad_emission_absorptionFvPatchScalarField::updateCoeffs()
{
  if (updated()) {
    return;
  }

  if(debug_) {
    Info << "--- qRad_emission_absorptionBoundaryConditions_.update(); --- Foam::qRad_emission_absorptionFvPatchScalarField::updateCoeffs()" << endl;
  }
  qRad_emission_absorptionBoundaryConditions_.update();

  if(debug_) {
    Info << "--- fixedValueFvPatchScalarField::updateCoeffs(); --- Foam::qRad_emission_absorptionFvPatchScalarField::updateCoeffs()" << endl;
  }
  fixedValueFvPatchScalarField::updateCoeffs();
  if(debug_) {
    Info << "--- end --- Foam::qRad_emission_absorptionFvPatchScalarField::updateCoeffs()" << endl;
  }
}


void Foam::qRad_emission_absorptionFvPatchScalarField::write(Ostream& os) const
{
  fvPatchScalarField::write(os);
  qRad_emission_absorptionBoundaryConditions_.write(os);
  writeEntry(os, "value", *this);
  //    writeEntryIfDifferent<scalar>(os, "inclinationAngle", -VGREAT, inclinationAngle_);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    qRad_emission_absorptionFvPatchScalarField
);
}

// ************************************************************************* //

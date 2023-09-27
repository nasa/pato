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

#include "pyro_recessionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pyro_recessionFvPatchScalarField::pyro_recessionFvPatchScalarField
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
pyro_recessionBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict_
)
{
  FatalErrorInFunction << "pyro_recessionFvPatchScalarField(p,iF) not implemented." << exit(FatalError);
}


Foam::pyro_recessionFvPatchScalarField::pyro_recessionFvPatchScalarField
(
    const pyro_recessionFvPatchScalarField& ptf,
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
pyro_recessionBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict_
)
{
}


Foam::pyro_recessionFvPatchScalarField::pyro_recessionFvPatchScalarField
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
pyro_recessionBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict
)
{
}


Foam::pyro_recessionFvPatchScalarField::pyro_recessionFvPatchScalarField
(
    const pyro_recessionFvPatchScalarField& frpsf
)
  :
fixedValueFvPatchScalarField(frpsf),
debug_(frpsf.debug_),
mesh(frpsf.mesh),
phaseName(frpsf.phaseName),
dictName(frpsf.dictName),
dict_(frpsf.dict_),
pyro_recessionBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict_
)
{
}


Foam::pyro_recessionFvPatchScalarField::pyro_recessionFvPatchScalarField
(
    const pyro_recessionFvPatchScalarField& frpsf,
    const DimensionedField<scalar, volMesh>& iF
)
  :
fixedValueFvPatchScalarField(frpsf, iF),
debug_(frpsf.debug_),
mesh(frpsf.mesh),
phaseName(frpsf.phaseName),
dictName(frpsf.dictName),
dict_(frpsf.dict_),
pyro_recessionBoundaryConditions_(frpsf.pyro_recessionBoundaryConditions_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::dictionary Foam::pyro_recessionFvPatchScalarField::initDict()
{
  dictionary dict_;
  FatalError << "Not implemented" << exit(FatalError);
  return dict_;
}

void Foam::pyro_recessionFvPatchScalarField::updateCoeffs()
{
  if (updated()) {
    return;
  }

  if(debug_) {
    Info << "--- pyro_recessionBoundaryConditions_.update(); --- Foam::pyro_recessionFvPatchScalarField::updateCoeffs()" << endl;
  }
  pyro_recessionBoundaryConditions_.update();

  if(debug_) {
    Info << "--- fixedValueFvPatchScalarField::updateCoeffs(); --- Foam::pyro_recessionFvPatchScalarField::updateCoeffs()" << endl;
  }
  fixedValueFvPatchScalarField::updateCoeffs();
  if(debug_) {
    Info << "--- end --- Foam::pyro_recessionFvPatchScalarField::updateCoeffs()" << endl;
  }
}


void Foam::pyro_recessionFvPatchScalarField::write(Ostream& os) const
{
  fvPatchScalarField::write(os);
  pyro_recessionBoundaryConditions_.write(os);
  writeEntry(os, "value", *this);
  //    writeEntryIfDifferent<scalar>(os, "inclinationAngle", -VGREAT, inclinationAngle_);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchScalarField,
    pyro_recessionFvPatchScalarField
);
}

// ************************************************************************* //

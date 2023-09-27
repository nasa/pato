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

#include "erosionModelFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::erosionModelFvPatchVectorField::erosionModelFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
  :
fixedValueFvPatchVectorField(p, iF),
debug_("no"),
mesh(patch().boundaryMesh().mesh()),
phaseName(word::null),
dictName(mesh.name()),
dict_(initDict()),
erosionModelBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict_
)
{
  FatalErrorInFunction << "erosionModelFvPatchVectorField(p,iF) not implemented." << exit(FatalError);
}


Foam::erosionModelFvPatchVectorField::erosionModelFvPatchVectorField
(
    const erosionModelFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
  :
fixedValueFvPatchVectorField(ptf, p, iF, mapper),
debug_(ptf.debug_),
mesh(ptf.mesh),
phaseName(ptf.phaseName),
dictName(ptf.dictName),
dict_(ptf.dict_),
erosionModelBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict_
)
{
}


Foam::erosionModelFvPatchVectorField::erosionModelFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
  :
fixedValueFvPatchVectorField(p, iF, dict),
debug_(dict.lookupOrDefault<Switch>("debug","no")),
mesh(patch().boundaryMesh().mesh()),
phaseName(word::null),
dictName(mesh.name()),
erosionModelBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict
)
{
}


Foam::erosionModelFvPatchVectorField::erosionModelFvPatchVectorField
(
    const erosionModelFvPatchVectorField& frpsf
)
  :
fixedValueFvPatchVectorField(frpsf),
debug_(frpsf.debug_),
mesh(frpsf.mesh),
phaseName(frpsf.phaseName),
dictName(frpsf.dictName),
dict_(frpsf.dict_),
erosionModelBoundaryConditions_(
    mesh,
    phaseName,
    dictName,
    patch().index(),
    dict_
)
{
}


Foam::erosionModelFvPatchVectorField::erosionModelFvPatchVectorField
(
    const erosionModelFvPatchVectorField& frpsf,
    const DimensionedField<vector, volMesh>& iF
)
  :
fixedValueFvPatchVectorField(frpsf, iF),
debug_(frpsf.debug_),
mesh(frpsf.mesh),
phaseName(frpsf.phaseName),
dictName(frpsf.dictName),
dict_(frpsf.dict_),
erosionModelBoundaryConditions_(frpsf.erosionModelBoundaryConditions_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::dictionary Foam::erosionModelFvPatchVectorField::initDict()
{
  dictionary dict_;
  dict_.add("physicsBasedErosionModel","1");
  return dict_;
}

void Foam::erosionModelFvPatchVectorField::updateCoeffs()
{
  if (updated()) {
    return;
  }

  if(debug_) {
    Info << "--- erosionModelBoundaryConditions_.update(); --- Foam::erosionModelFvPatchVectorField::updateCoeffs()" << endl;
  }
  erosionModelBoundaryConditions_.update();

  if(debug_) {
    Info << "--- fixedValueFvPatchVectorField::updateCoeffs(); --- Foam::erosionModelFvPatchVectorField::updateCoeffs()" << endl;
  }
  fixedValueFvPatchVectorField::updateCoeffs();
  if(debug_) {
    Info << "--- end --- Foam::erosionModelFvPatchVectorField::updateCoeffs()" << endl;
  }
}


void Foam::erosionModelFvPatchVectorField::write(Ostream& os) const
{
  fvPatchVectorField::write(os);
  erosionModelBoundaryConditions_.write(os);
  writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField
(
    fvPatchVectorField,
    erosionModelFvPatchVectorField
);
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2020 PATO
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

Class
    Foam::forchheimerVelocityPressureFvPatchScalarField

Description
    This boundary condition provides a Darcy's velocity condition for pressure.
    It is a fixedGradient condition whereas pressure is calculated from the
    projection of the velocity vector and permeability tensor on the direction
    of the flow.

Usage
    \table
        Property | Description                 | Required | Default value
        U        | Darcy's velocity            | no       | U
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchname>
    {
        type forchheimerVelocityPressure;
        U    uniform (0.5 0 0);
    }
    \endvaerbatim

SourceFiles
    forchheimerVelocityPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef forchheimerVelocityPressureFvPatchScalarFields_H
#define forchheimerVelocityPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class forchheimerVelocityPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class forchheimerVelocityPressureFvPatchScalarField
:
public fixedGradientFvPatchScalarField
{
  // Private data

  //- Darcy's velocity field
  vectorField U_;

 public:

  //- Runtime type information
  TypeName("forchheimerVelocityPressure");


  // Constructors

  //- Construct from patch and internal field
  forchheimerVelocityPressureFvPatchScalarField
  (
      const fvPatch&,
      const DimensionedField<scalar, volMesh>&
  );

  //- Construct from patch, internal field and dictionary
  forchheimerVelocityPressureFvPatchScalarField
  (
      const fvPatch&,
      const DimensionedField<scalar, volMesh>&,
      const dictionary&
  );

  //- Construct by mapping given forchheimerVelocityPressureFvPatchScalarField onto
  //  a new patch
  forchheimerVelocityPressureFvPatchScalarField
  (
      const forchheimerVelocityPressureFvPatchScalarField&,
      const fvPatch&,
      const DimensionedField<scalar, volMesh>&,
      const fvPatchFieldMapper&
  );

  //- Copy constructor
  forchheimerVelocityPressureFvPatchScalarField
  (
      const forchheimerVelocityPressureFvPatchScalarField&
  );

  //- Construct and return a clone
  virtual tmp<fvPatchScalarField> clone() const {
    return tmp<fvPatchScalarField>
           (
               new forchheimerVelocityPressureFvPatchScalarField(*this)
           );
  }

  //- Copy constructor setting internal field reference
  forchheimerVelocityPressureFvPatchScalarField
  (
      const forchheimerVelocityPressureFvPatchScalarField&,
      const DimensionedField<scalar, volMesh>&
  );

  //- Construct and return a clone setting internal field reference
  virtual tmp<fvPatchScalarField> clone
  (
      const DimensionedField<scalar, volMesh>& iF
  ) const {
    return tmp<fvPatchScalarField>
           (
               new forchheimerVelocityPressureFvPatchScalarField(*this, iF)
           );
  }


  // Member functions

  //- Update the coefficients associated with the patch field
  virtual void updateCoeffs();

  //- Write
  virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

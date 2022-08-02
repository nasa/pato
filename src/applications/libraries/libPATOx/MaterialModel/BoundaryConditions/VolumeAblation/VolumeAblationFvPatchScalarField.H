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

Class
    Foam::VolumeAblationFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    Foam::VolumeAblationFvPatchScalarField

    This boundary condition provides a fixed boundary conditions from
    the surface energy balance.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
         type            surfaceEnergyBalance;
	     ...
    }
    \endverbatim

SourceFile
    VolumeAblationFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef VolumeAblationFvPatchScalarField_H
#define VolumeAblationFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "fvCFD.H"
#include "VolumeAblationBoundaryConditions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//class BprimeBoundaryConditions;
/*---------------------------------------------------------------------------*\
                 Class VolumeAblationFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class VolumeAblationFvPatchScalarField
:
public fixedValueFvPatchScalarField
{
  //- debug switch
  Switch debug_;

  //- Reference to the mesh database
  const fvMesh& mesh;

  //- Phase name
  const word phaseName;

  //- Region name
  const word dictName;

  //- Dictionary
  dictionary dict_;

 public:

  //- Volume Ablation boundary conditions object
  VolumeAblationBoundaryConditions VolumeAblationBoundaryConditions_;

  //- Runtime type information
  TypeName("VolumeAblation");


  // Constructors

  //- Construct from patch and internal field
  VolumeAblationFvPatchScalarField
  (
      const fvPatch&,
      const DimensionedField<scalar, volMesh>&
  );

  //- Construct from patch, internal field and dictionary
  VolumeAblationFvPatchScalarField
  (
      const fvPatch&,
      const DimensionedField<scalar, volMesh>&,
      const dictionary&
  );

  //- Construct by mapping given VolumeAblationFvPatchScalarField
  // onto a new patch
  VolumeAblationFvPatchScalarField
  (
      const VolumeAblationFvPatchScalarField&,
      const fvPatch&,
      const DimensionedField<scalar, volMesh>&,
      const fvPatchFieldMapper&
  );

  //- Construct as copy
  VolumeAblationFvPatchScalarField
  (
      const VolumeAblationFvPatchScalarField&
  );

  //- Construct and return a clone
  virtual tmp<fvPatchScalarField> clone() const {
    return tmp<fvPatchScalarField>
           (
               new VolumeAblationFvPatchScalarField(*this)
           );
  }

  //- Construct as copy setting internal field reference
  VolumeAblationFvPatchScalarField
  (
      const VolumeAblationFvPatchScalarField&,
      const DimensionedField<scalar, volMesh>&
  );

  //- Construct and return a clone setting internal field reference
  virtual tmp<fvPatchScalarField> clone
  (
      const DimensionedField<scalar, volMesh>& iF
  ) const {
    return tmp<fvPatchScalarField>
           (
               new VolumeAblationFvPatchScalarField(*this, iF)
           );
  }


  // Member functions

  // Evaluation functions

  //- Update the coefficients associated with the patch field
  virtual void updateCoeffs();

  //- Write
  virtual void write(Ostream&) const;

  //- Initialize
  dictionary initDict();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

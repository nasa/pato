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
    Foam::erosionModelFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Foam::erosionModelFvPatchVectorField

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
    erosionModelFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef erosionModelFvPatchVectorField_H
#define erosionModelFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "fvCFD.H"
#include "erosionModelBoundaryConditions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class erosionModelFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class erosionModelFvPatchVectorField
:
public fixedValueFvPatchVectorField
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
  erosionModelBoundaryConditions erosionModelBoundaryConditions_;

  //- Runtime type information
  TypeName("erosionModel");


  // Constructors

  //- Construct from patch and internal field
  erosionModelFvPatchVectorField
  (
      const fvPatch&,
      const DimensionedField<vector, volMesh>&
  );

  //- Construct from patch, internal field and dictionary
  erosionModelFvPatchVectorField
  (
      const fvPatch&,
      const DimensionedField<vector, volMesh>&,
      const dictionary&
  );

  //- Construct by mapping given erosionModelFvPatchVectorField
  // onto a new patch
  erosionModelFvPatchVectorField
  (
      const erosionModelFvPatchVectorField&,
      const fvPatch&,
      const DimensionedField<vector, volMesh>&,
      const fvPatchFieldMapper&
  );

  //- Construct as copy
  erosionModelFvPatchVectorField
  (
      const erosionModelFvPatchVectorField&
  );

  //- Construct and return a clone
  virtual tmp<fvPatchVectorField> clone() const {
    return tmp<fvPatchVectorField>
           (
               new erosionModelFvPatchVectorField(*this)
           );
  }

  //- Construct as copy setting internal field reference
  erosionModelFvPatchVectorField
  (
      const erosionModelFvPatchVectorField&,
      const DimensionedField<vector, volMesh>&
  );

  //- Construct and return a clone setting internal field reference
  virtual tmp<fvPatchVectorField> clone
  (
      const DimensionedField<vector, volMesh>& iF
  ) const {
    return tmp<fvPatchVectorField>
           (
               new erosionModelFvPatchVectorField(*this, iF)
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

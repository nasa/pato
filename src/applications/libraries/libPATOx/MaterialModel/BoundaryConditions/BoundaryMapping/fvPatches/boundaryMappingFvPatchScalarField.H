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
    Foam::boundaryMappingFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    Foam::boundaryMappingFvPatchScalarField

    This boundary condition provides a fixed pressure boundary conditions from the environment table.


Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            boundaryMapping;
	    mappingType constant;
	    mappingFileName "PATH_to_table";
        mappingFields ((p 1));
    }
    \endverbatim

SourceFile
    boundaryMappingFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryMappingFvPatchScalarField_H
#define boundaryMappingFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "fvCFD.H"
#include "PATOx.H"
#include "simpleBoundaryMappingModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class boundaryMappingFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class boundaryMappingFvPatchScalarField
:
public fixedValueFvPatchScalarField
{
 private:
  //- Reference to the mesh database
  const fvMesh& mesh_;

  //- Needed fields for boundaryMapping object
  const wordList neededFields_;

  //- Dictionary
  const dictionary dict_;

  //- mapping fields for boundaryMapping object
  const List<Tuple2<word,word> > mappingFields_;

  //- Name of the patch
  const word patchName_;

 public:

  //- Boundary Mapping object
  autoPtr<simpleBoundaryMappingModel> boundaryMapping_;
  simpleBoundaryMappingModel * boundaryMapping_ptr;

  //- Runtime type information
  TypeName("boundaryMapping");

  // Constructors

  //- Construct from patch and internal field
  boundaryMappingFvPatchScalarField
  (
      const fvPatch&,
      const DimensionedField<scalar, volMesh>&
  );

  //- Construct from patch, internal field and dictionary
  boundaryMappingFvPatchScalarField
  (
      const fvPatch&,
      const DimensionedField<scalar, volMesh>&,
      const dictionary&
  );

  //- Construct by mapping given boundaryMappingFvPatchScalarField
  // onto a new patch
  boundaryMappingFvPatchScalarField
  (
      const boundaryMappingFvPatchScalarField&,
      const fvPatch&,
      const DimensionedField<scalar, volMesh>&,
      const fvPatchFieldMapper&
  );

  //- Construct as copy
  boundaryMappingFvPatchScalarField
  (
      const boundaryMappingFvPatchScalarField&
  );

  //- Construct and return a clone
  virtual tmp<fvPatchScalarField> clone() const {
    return tmp<fvPatchScalarField>
           (
               new boundaryMappingFvPatchScalarField(*this)
           );
  }

  //- Construct as copy setting internal field reference
  boundaryMappingFvPatchScalarField
  (
      const boundaryMappingFvPatchScalarField&,
      const DimensionedField<scalar, volMesh>&
  );

  //- Construct and return a clone setting internal field reference
  virtual tmp<fvPatchScalarField> clone
  (
      const DimensionedField<scalar, volMesh>& iF
  ) const {
    return tmp<fvPatchScalarField>
           (
               new boundaryMappingFvPatchScalarField(*this, iF)
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


  //- Initialize patchName
  word initPatchName();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

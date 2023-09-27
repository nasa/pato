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
    Foam::twoDaxiBoundaryMappingModel

Description

SourceFiles
    twoDaxiBoundaryMappingModel.H
    twoDaxiBoundaryMappingModel.C

\*---------------------------------------------------------------------------*/

#ifndef twoDaxiBoundaryMappingModel_H
#define twoDaxiBoundaryMappingModel_H

#include "simpleBoundaryMappingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class twoDaxiBoundaryMappingModel Declaration
\*---------------------------------------------------------------------------*/

class twoDaxiBoundaryMappingModel
:
public simpleBoundaryMappingModel
{

 protected:

  //- axis of the mapping file(s) (2D-axi)
  word mappingAxis_;

  //- symmetry of the mapping field (3D-tecplot)
  word mappingSymmetry_;

  //- times of the mapping file(s)
  scalarList timesMappingFileData_;

  //- list of data from the mapping file(s)
  List<List<scalarList> > listMappingFileData_;

 public:

  //- Runtime type information
  TypeName("2D-axi");

  // Constructors

  //- Construct from mesh and phase name
  twoDaxiBoundaryMappingModel(      const fvMesh& mesh,
                                    const wordList& neededFields,
                                    const dictionary dict);

  //- Destructor
  virtual ~twoDaxiBoundaryMappingModel();

  //- Member Functions

  //- Update
  virtual void update(scalar timeValue, label patchID, word fieldName);

  //- write in os
  virtual void write(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

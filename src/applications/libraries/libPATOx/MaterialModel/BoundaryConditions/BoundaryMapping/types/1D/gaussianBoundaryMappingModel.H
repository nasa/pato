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
    Foam::gaussianBoundaryMappingModel

Description

SourceFiles
    gaussianBoundaryMappingModel.H
    gaussianBoundaryMappingModel.C

\*---------------------------------------------------------------------------*/

#ifndef gaussianBoundaryMappingModel_H
#define gaussianBoundaryMappingModel_H

#include "simpleBoundaryMappingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class gaussianBoundaryMappingModel Declaration
\*---------------------------------------------------------------------------*/

class gaussianBoundaryMappingModel
:
public simpleBoundaryMappingModel
{

 protected:

  //- Field coefficients
  List<scalarList> fieldCoeffs_;

 public:

  //- Runtime type information
  TypeName("gaussian");

  // Constructors

  //- Construct from mesh and phase name
  gaussianBoundaryMappingModel
  (
      const fvMesh& mesh,
      const wordList& neededFields,
      const dictionary dict
  );

  //- Destructor
  virtual ~gaussianBoundaryMappingModel();

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

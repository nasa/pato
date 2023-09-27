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
    Foam::twoDaxiInterpolateNearestBoundaryMappingModel

Description

SourceFiles
    twoDaxiInterpolateNearestBoundaryMappingModel.H
    twoDaxiInterpolateNearestBoundaryMappingModel.C

\*---------------------------------------------------------------------------*/

#ifndef twoDaxiInterpolateNearestBoundaryMappingModel_H
#define twoDaxiInterpolateNearestBoundaryMappingModel_H

#include "simpleBoundaryMappingModel.H"
#include "mathFunctions.H"
#include "IOFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class twoDaxiInterpolateNearestBoundaryMappingModel Declaration
\*---------------------------------------------------------------------------*/

class twoDaxiInterpolateNearestBoundaryMappingModel
:
public simpleBoundaryMappingModel
{

 protected:

  //- Update if mag(time - updatedTime) > updateTolTime, default value = 0.1 s
  scalar updateTolTime_;

  //- times of the mapping file(s)
  scalarList timesMappingFileData_;

  //- list of data from the mapping file(s)
  List<List<scalarList> > listMappingFileData_;

  //- Mutiplier coefficient
  scalar multiplier_;

 public:

  //- Runtime type information
  TypeName("2D-axi_InterpolateNearest");

  // Constructors

  //- Construct from mesh and phase name
  twoDaxiInterpolateNearestBoundaryMappingModel(const fvMesh& mesh,
      const wordList& neededFields,
      const dictionary dict);

  //- Destructor
  virtual ~twoDaxiInterpolateNearestBoundaryMappingModel();

  //- Member Functions

  //- Update
  virtual void update(scalar timeValue, label patchID, word fieldName);

  //- Interpolate from the nearest points (3D-tecplot)
  scalar interpolateNearest(label patchID, label faceI, label fieldI, List<scalarList>& mappingFileData);

  //- write in os
  virtual void write(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

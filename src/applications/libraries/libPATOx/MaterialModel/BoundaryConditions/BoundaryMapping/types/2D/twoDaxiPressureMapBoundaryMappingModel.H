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
    Foam::twoDaxiPressureMapBoundaryMappingModel

Description

SourceFiles
    twoDaxiPressureMapBoundaryMappingModel.H
    twoDaxiPressureMapBoundaryMappingModel.C

\*---------------------------------------------------------------------------*/

#ifndef twoDaxiPressureMapBoundaryMappingModel_H
#define twoDaxiPressureMapBoundaryMappingModel_H

#include "simpleBoundaryMappingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class twoDaxiPressureMapBoundaryMappingModel Declaration
\*---------------------------------------------------------------------------*/

class twoDaxiPressureMapBoundaryMappingModel
:
public simpleBoundaryMappingModel
{

 protected:

  //- dynamic pressure field name - 2D axi pressureMap only
  word dynamicPressureFieldName_;

  //- dynamic pressure field - 2D axi pressureMap only
  volScalarField * p_dyn_ptr;

  //- list of data from the mapping file(s)
  List<List<scalarList> > listMappingFileData_;

 public:

  //- Runtime type information
  TypeName("2D-axi_pressureMap");

  // Constructors

  //- Construct from mesh and phase name
  twoDaxiPressureMapBoundaryMappingModel(      const fvMesh& mesh,
      const wordList& neededFields,
      const dictionary dict);

  //- Destructor
  virtual ~twoDaxiPressureMapBoundaryMappingModel();

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

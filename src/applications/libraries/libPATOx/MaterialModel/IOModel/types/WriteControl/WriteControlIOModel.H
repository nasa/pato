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
    Foam::WriteControlIOModel

Description

SourceFiles
    WriteControlIOModel.H
    WriteControlIOModel.C

\*---------------------------------------------------------------------------*/

#ifndef WriteControlIOModel_H
#define WriteControlIOModel_H

#include "simpleIOModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class WriteControlIOModel Declaration
\*---------------------------------------------------------------------------*/

class WriteControlIOModel
:
public simpleIOModel
{

 protected:

  //- Mesh reference
  const fvMesh& mesh_;

  //- Phase name
  const word phaseName_;

  //- Dictionary Name
  const word& dictName_;

  //- Output List
  const wordList& outputList_;

  //- Init flag
  const bool initOutput_;

  //- List of times to write the output
  scalarList writeControlTimes;

  //- List of the write control flags: if yes, already wrote the output
  List<Switch> writeControlFlags;

 public:

  //- Runtime type information
  TypeName("WriteControl");

  // Constructors

  //- Construct from mesh and phase name
  WriteControlIOModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~WriteControlIOModel();

  //- Member Functions

  //- Update
  virtual void update();

  //- Write the output
  void writeOutput();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//#ifdef NoRepository
//#include "WriteControlIOModel.C"
//#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
#include "simpleMaterialPropertiesModel.H"

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

  //- Output List
  const wordList& outputList_;

  //- MaterialProperties Model
  simpleMaterialPropertiesModel& materialPropertiesModel_;

  //- Initialize the material properties
  Switch initMatProp_;

  //- Init flag
  const bool initOutput_;

  //- List of user specified times at which to write
  scalarList userWriteTimes;

  //- List of flags for user specified write times: if yes, already wrote the output
  List<Switch> userWriteFlags;

  //- List of user specified times to output for plotting
  scalarList userProbeTimes;

  //- List of flags for user specified probing times: if yes, already wrote the output
  List<Switch> userProbeFlags;

 public:

  //- Runtime type information
  TypeName("WriteControl");

  // Constructors

  //- Construct from mesh and phase name
  WriteControlIOModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~WriteControlIOModel();

  // Member Functions

  //- Initialize the material properties
  Switch initMatProp();

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

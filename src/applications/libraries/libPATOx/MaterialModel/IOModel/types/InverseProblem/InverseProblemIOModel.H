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
    along with OpenFOAM.  If InverseProblemt, see <http://www.gnu.org/licenses/>.

Class
    Foam::InverseProblemIOModel

Description

SourceFiles
    InverseProblemIOModel.H
    InverseProblemIOModel.C

\*---------------------------------------------------------------------------*/

#ifndef InverseProblemIOModel_H
#define InverseProblemIOModel_H

#include "noIOModel.H"
#include "mathFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class InverseProblemIOModel Declaration
\*---------------------------------------------------------------------------*/

class InverseProblemIOModel
:
public noIOModel
{

 protected:

  //- Mesh reference
  const fvMesh& mesh_;

  //- Dictionary name
  const word& dictName_;

  //- List of list null
  const List<wordList> nullListList_;

  //- runTimeModifiable in controlDict
  Switch runTimeModifiable_;

  //- writeInterval in controlDict
  scalar writeInterval_;

  //- writeControl in controlDict
  word writeControl_;

  //- Temperature [K]
  volScalarField& Ta_;

  // Output Model: READ
  //- List of the fields to write in the time folders
  const wordList& outputList_;
  //- List of the files to read
  const List<fileName>& readFilesList_;
  //- List of the data in the read files
  const List<List<scalarList> >& filesData_;

  // OutputModel: COMPARE
  //- list of the (file fields) to compare
  const List<Tuple2<Tuple2<word,scalar>,Tuple2<word,vector> > > compareList_;
  //- list of the fields name to compare
  wordList compareFields_;
  //- list of the fields probing coordinates
  List<vector> compareFieldsCoord_;
  //- list of the read file column to compare
  scalarList compareReadFileColumn_;
  //- list of the read file to compare
  wordList compareReadFile_;
  //- list of the read file data to compare
  List<List<scalarList> > compareReadFileData_;
  //- deltaTime interval between the results
  scalar compareTimeInterval_;
  //- Dakota method (MOGA or nl2sol)
  const word dakotaMethod_;
  //- Dakota results
  scalarList dakotaResults_;
  // Dakota results file name
  fileName dakotaResultsFileName_;
  //- List of the Sample Functions which are compared
  PtrList<SampleFunction> compareSampleFunctions_;
  //- indexes of sampling coordinates
  List<int> indexSamplingCoord_;
  //- indexes of sampling fields
  List<int> indexSamplingFields_;
  //- output delta time
  scalar outputDeltaTime_;
  //- output time
  scalar outputTime_;

  // OutputModel: INITIALIZE
  //- Switch to initialize the basic output model
  const bool initOutput_;

 public:

  //- Runtime type information
  TypeName("InverseProblem");

  // Constructors

  //- Construct from mesh and phase name
  InverseProblemIOModel(const fvMesh& mesh,const word& dictName);

//  //- Destructor
  virtual ~InverseProblemIOModel();

  //- Member Functions

  //- Update
  virtual void update();

  //- Update the dakota results
  inline void updateDakotaResults();

  //- Read the input files
  inline void readInput();

  //- Initialize Dakota input
  inline void dakotaInit();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::simpleIOModel

Description

SourceFiles
    simpleIOModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleIOModel_H
#define simpleIOModel_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "IOFunctions.H"
#include "SampleFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleIOModel Declaration
\*---------------------------------------------------------------------------*/

class simpleIOModel
:
public IOdictionary
{

 private:

  // Private Member Functions

 protected:

  // Protected data

  //- Mesh database
  const fvMesh& mesh_;

  //- Phase name
  const word phaseName_;

  //- Dictionary name
  const word& dictName_;

  //- Material dictionary
  IOdictionary materialDict_;

  //- Debug Switch
  const Switch debug_;

  //- Dynamic mesh flag
  const bool dynamicMesh_;

  //-  lduMatrix::debug
  label infoDebug_;

  //- Material dictionary path
  fileName materialDictPath_;

  //- empty list
  const wordList nullList;

  //- empty list
  const List<fileName> nullListFileName;

  //- list of the fields to write
  wordList outputList_;

  //- Name of the probing dictionary
  wordList probingDictNames_;

  // List of the Sample Functions which print the output folder
  PtrList<SampleFunction> sampleFunctions_;

  //- list of the files to read
  List<fileName> readFilesList_;

  //- list of the data in the read files
  List<List<scalarList> > filesData_;

  // Protected Member Functions

  //- Correct function - updates due to mesh changes
  void correct();
 public:

  //- Runtime type information
  TypeName("simpleIOModel");


  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleIOModel,
      fvMesh,
      (const fvMesh& mesh,const word& dictName),
      (mesh, dictName)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleIOModel
  (
      const fvMesh& mesh,
      const word& dictName
  );


  //- Selector
  static autoPtr<simpleIOModel> New
  (
      const fvMesh& mesh,
      const word& dictName
  );

  //- Destructor
  virtual ~simpleIOModel();

  //- Update
  virtual void update() = 0;

  //- Return the Info debug label
  inline const label& infoDebug() const;

  //- Read the list of files
  inline const List<fileName>& readFilesList() const;

  //- Data in the files
  inline const List<List<scalarList> >& filesData() const;

  //- Initialize the output model
  inline bool initOutput();


};

#include "simpleIOModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

#define namespaceFoam
#include "simpleModel.H"
#include "SampleFunction.H"
#include "simpleMaterialsModel.H"

#if defined(FOAM_EXTEND)
#undef FALSE
#undef TRUE
#define namespaceFoam_extend_
#include "declare_macros_foam_extend.H" // Redeclare Foam Extend 4.1 macros
#include "Foam_extend_fvCFD.H"
#include "Foam_extend_dynamicFvMesh.H"
#include "declare_macros_openfoam.H" // Redeclare OpenFOAM 7 macros
#endif
using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleIOModel Declaration
\*---------------------------------------------------------------------------*/

class simpleIOModel
:
public simpleModel
{

 private:

  // Private Member Functions

 protected:

  // Protected data

#if defined(FOAM_EXTEND)
  //- Foam Extend mesh database
  const Foam_extend_::fvMesh& foam_extend_mesh_;
#endif

  //- Phase name
  const word phaseName_;

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

  //- Write fields flag in controlDict
  const Switch writeFields_;

#if defined(FOAM_EXTEND)
  //- List of the Foam Extend fields to write
  wordList foam_extend_write_fields_;
#endif

  // Protected Member Functions

  //- Correct function - updates due to mesh changes
  void correct();

 public:

  //- Runtime type information
  TypeName("simpleIOModel");

  //- Model name
  static const word modelName;

  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleIOModel,
      fvMesh,
      (const fvMesh& mesh,const word& regionName),
      (mesh, regionName)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleIOModel
  (
      const fvMesh& mesh,
      const word& regionName
  );


  //- Selector
  static autoPtr<simpleIOModel> New
  (
      const fvMesh& mesh,
      const word& regionName
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

  //- Write fields flag in controlDict
  inline const Switch& writeFields() const;

  //- Initialize the output model
  virtual bool initOutput();

#if defined(FOAM_EXTEND)
  //- Initialize the Foam Extend fields to write
  virtual void initFoamExtendFields();

  //- Update the Foam Extend fields to write
  virtual void updateFoamExtendFields();
#endif


};

#include "simpleIOModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

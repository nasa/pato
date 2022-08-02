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
    Foam::SampleFunction

Description
    A base class for IO model

SourceFiles
    SampleFunction.C

\*---------------------------------------------------------------------------*/

#ifndef SampleFunction_H
#define SampleFunction_H

#include "fvCFD.H"
#include "functionObjectListTest.H"
#include "dictionaryEntry.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "meshSearch.H"
#include "IOFunctions.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SampleFunction Declaration
\*---------------------------------------------------------------------------*/

class SampleFunction
{
 protected:
  //- Reference to the mesh database
  const fvMesh& mesh_;

  //- Mesh searching enigne
  autoPtr<meshSearch> meshSearchPtr_;

  //- Run Time
  const Time& runTime;

  //- Name of the probing dictionary
  const word& probingDictName_;

  //- List of the sub dictionaries name
  PtrList<word> subProbingDictTOC_;

  //- List of the sub dictionaries
  List<dictionary> listSubProbingDict_;

  //- list of the points to probe per sub dictionary in the probing dictionary
  PtrList<List<vector> > probingPoints_;

  //- list of the scalar fields to probe
  PtrList<word> probingFields_scalar_;

  //- list of the vector fields to probe
  PtrList<word> probingFields_vector_;

  //- list of the tensor fields to probe
  PtrList<word> probingFields_tensor_;

  //- list of the scalar fields to probe
  PtrList<PtrList<OFstream> > os_out_scalar_;

  //- list of the vector fields to probe
  PtrList<PtrList<PtrList<OFstream> > > os_out_vector_;

  //- list of the tensor fields to probe
  PtrList<PtrList<PtrList<OFstream> > > os_out_tensor_;

  //- List of the function object
  autoPtr<functionObjectListTest> functionPtr_;

  // Reference to the fields list
  autoPtr<List<List<Field<scalar> > > > scalarMasterFields_ref;
  autoPtr<List<List<Field<vector> > > > vectorMasterFields_ref;
  autoPtr<List<List<Field<sphericalTensor> > > > sphericalTensorMasterFields_ref;
  autoPtr<List<List<Field<symmTensor> > > > symmTensorMasterFields_ref;
  autoPtr<List<List<Field<tensor> > > > tensorMasterFields_ref;

  // Path to the dictionary
  fileName dictPath_;

  // Initial cell centers
  volVectorField C_initial;

  //- dynamic mesh flag
  const bool dynamicMesh_;

  // function dictionary
  dictionary functionDict_;

  //- processor rank
  const int rank_;

  //- First iteration flag
  Switch firstIteration_;

 public:

  // Constructors

  //- Construct from dictionary name
  SampleFunction
  (
      const fvMesh& mesh,
      const word& phaseName,
      const word& dictName,
      const word& probingDictName
  );


  //- Destructor
  virtual ~SampleFunction();
  inline void writeOutput();
  inline void updateProbing();


  //- list of the points to probe per sub dictionary in the probing dictionary
  inline PtrList<List<vector> >& probingPoints();

  //- list of the scalar fields to probe
  inline PtrList<word>& probingFields_scalar();

  //- list of the vector fields to probe
  inline PtrList<word>& probingFields_vector();

  //- list of the tensor fields to probe
  inline PtrList<word>& probingFields_tensor();

  //- Path to the dictionary
  inline fileName& dictPath();

  //- Reference to the scalar fields list
  inline autoPtr<List<List<Field<scalar> > > >& scalarMasterFields();

  //- Change the string to lowercase
  inline void lower_string(string& str);

  //- Change the name of the element list if the uppercase/lowercase name exists already in the list
  inline void changeUpperLowercaseName(PtrList<word>& list);

  //- Verify if the point is in the mesh
  inline bool point_in_mesh(const vector& point);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SampleFunctionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

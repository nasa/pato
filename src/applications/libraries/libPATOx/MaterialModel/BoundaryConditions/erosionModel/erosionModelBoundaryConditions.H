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
    along with OpenFOAM.  If Bprimet, see <http://www.gnu.org/licenses/>.

Class
    Foam::erosionModelBoundaryConditions

Description
    Physic based erosion model for cellMotionU (volume ablation)

SourceFiles
    erosionModelBoundaryConditions.H
    erosionModelBoundaryConditions.C

\*---------------------------------------------------------------------------*/

#ifndef erosionModelBoundaryConditions_H
#define erosionModelBoundaryConditions_H


#include "fvCFD.H"
#include "mathFunctions.H"
#include "IOmanip.H"
#include "IOFunctions.H"
#include "PATOx.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class erosionModelBoundaryConditions Declarations
\*---------------------------------------------------------------------------*/

class erosionModelBoundaryConditions
{
 protected:

  //- Reference to the mesh database
  const fvMesh& mesh_;

  //- Phase name
  const word phaseName_;

  //- Dictionary name
  const word regionName_;

  //- Current patch ID
  const label currentPatchID_;

  //- Dictionary name
  const dictionary dict_;

  //- True if mesh is a dynamic mesh
  Switch dynamicMesh_;

  //- Initialized flag
  Switch initialized_;

  //- volScalarFields
  PMap<volScalarField> scalarFields_;

  //- volVectorFields
  PMap<volVectorField> vectorFields_;

  //- volTensorFields
  PMap<volTensorField> tensorFields_;

  //- dimScalars
  PMap<dimensionedScalar> dimScalars_;

  // Component of dynamic mesh
  word component_;

  //- Physics based erosion model
  const int physicsBasedErosionModel_;

  //- Debug switch
  Switch debug_;

  //- Cell motion velocity of the dynamic mesh [m/s]
  volVectorField * cellMotionU_ptr;

 public:

  // Constructors

  //- Construct
  erosionModelBoundaryConditions(const fvMesh& mesh, const word& phaseName, const word& regionName, const label& currentPatchID, const dictionary dict);

  //- Destructor
  virtual ~erosionModelBoundaryConditions();

  // Member Functions

  //- Initialize
  void initialize();

  //- Update
  virtual void update();

  //- write in os
  void write(Ostream& os) const;

  //- update the motion boundary conditions
  void updateMotionBC();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

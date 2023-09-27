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
    Foam::volumeAverageMaterialFailureMassRemovalModel

Description

SourceFiles
    volumeAverageMaterialFailureMassRemovalModel.H
    volumeAverageMaterialFailureMassRemovalModel.C

\*---------------------------------------------------------------------------*/

#ifndef volumeAverageMaterialFailureMassRemovalModel_H
#define volumeAverageMaterialFailureMassRemovalModel_H

#include "simpleMaterialFailureMassRemovalModel.H"
#include "meshSearch.H"
#include "wedgePolyPatch.H"
#include "IOmanip.H"
#include "simpleSolidMechanicsModel.H"
#include "simpleEnergyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class volumeAverageMaterialFailureMassRemovalModel Declaration
\*---------------------------------------------------------------------------*/

class volumeAverageMaterialFailureMassRemovalModel
:
public simpleMaterialFailureMassRemovalModel
{
 public:

  //- Runtime type information
  TypeName("volumeAverage"); // Disable material failure

  // Constructors

  //- Construct from mesh and phase name
  volumeAverageMaterialFailureMassRemovalModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~volumeAverageMaterialFailureMassRemovalModel();

  //- Member Functions

  //- Update material failure mass removal
  virtual void updateMassRemoval();

  //- Distance of the failure criteria from the surface
  scalar distanceFromSurface(label patchI, label faceI);

  //- Create the list of boundary types for the failure criteria
  wordList failureCriteriaBC();

  //- Clean the fields
  virtual void cleanFields();

 private:

  //- List of boundary types for the failure criteria
  wordList failureCriteriaBC_;

  // Solid Mechanics Model
  simpleSolidMechanicsModel& solidMechanicsModel_;

  // Failure criteria field
  volScalarField& failureCriteria_;

  // Failure mass loss [kg]
  volScalarField& failureMassLoss_;

  // Energy Model
  simpleEnergyModel& energyModel_;

  // Solid density [kg/m3]
  volScalarField& rho_s_;

  //- Cell motion velocity [m/s]
  volVectorField * cellMotionU_ptr;

  // Total failure mass loss [kg]
  scalar totalFailureMassLoss_;

  // Initial total mass [kg]
  scalar initialTotalMass_;

  // Failure criteria
  const int FAILS = 1;

  //- Mesh searching enigne
  autoPtr<meshSearch> meshSearchPtr_;

  //- Output name for the recession rate [m/s]
  fileName failureRecessionRateOutputName_;

  //- OFstream width
  int width_;

  //- OFstream precision
  int precision_;

  //- Output file for the recession rate [m/s]
  OFstream * os_failureRecessionRate_ptr;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

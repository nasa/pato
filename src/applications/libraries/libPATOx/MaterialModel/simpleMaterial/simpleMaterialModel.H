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
    Foam::simpleMaterialModel

Description

SourceFiles
    simpleMaterialModel.C

\*---------------------------------------------------------------------------*/

#ifndef simpleMaterialModel_H
#define simpleMaterialModel_H

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "IOFunctions.H"
#include "simpleMaterialPropertiesModel.H"
#include "simpleGasPropertiesModel.H"
#include "simpleMassModel.H"
#include "simpleEnergyModel.H"
#include "simpleSolidMechanicsModel.H"
#include "simpleIOModel.H"
#include "simpleMaterialChemistryModel.H"
#include "simplePyrolysisModel.H"
#include "simpleVolumeAblationModel.H"
#include "simpleTimeControlModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class simpleMaterialModel Declaration
\*---------------------------------------------------------------------------*/

class simpleMaterialModel
:
public IOdictionary
{

 private:

  // Private Member Functions

 protected:

  // Protected data

  //- Mesh database
  const fvMesh& mesh_;

  //- Material dictionary
  IOdictionary materialDict_;

  //- Debug Switch
  const Switch debug_;

  //- Dynamic mesh flag
  const bool dynamicMesh_;

  // Models
  simpleMaterialChemistryModel& materialChemistryModel_;
  simplePyrolysisModel& pyrolysisModel_;
  simpleMaterialPropertiesModel& materialPropertiesModel_;
  simpleGasPropertiesModel& gasPropertiesModel_;
  simpleMassModel& massModel_;
  simpleEnergyModel& energyModel_;
  simpleSolidMechanicsModel& solidMechanicsModel_;
  simpleVolumeAblationModel& volumeAblationModel_;
  simpleIOModel& ioModel_;
  simpleTimeControlModel& timeControlModel_;

 public:

  //- Runtime type information
  TypeName("simpleMaterialModel");

  scalar updateMinDeltaTw_; // NAmrofel 26/11/2019 

  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      simpleMaterialModel,
      fvMesh,
      (const fvMesh& mesh,const word& dictName),
      (mesh, dictName)
  );


  // Constructors

  //- Construct from mesh and phase name
  simpleMaterialModel
  (
      const fvMesh& mesh,
      const word& dictName
  );


  //- Destructor
  virtual ~simpleMaterialModel();

  //- Loop over all cells and update all thermodynamic properties
  void update();
};

#include "simpleMaterialModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

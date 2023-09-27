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
    Foam::massIOModel

Description

SourceFiles
    massIOModel.H
    massIOModel.C

\*---------------------------------------------------------------------------*/

#ifndef massIOModel_H
#define massIOModel_H

#include "noIOModel.H"
#include "simplePyrolysisModel.H"
#include "simpleMassModel.H"
#include "simpleEnergyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class massIOModel Declaration
\*---------------------------------------------------------------------------*/

class massIOModel
:
public noIOModel
{

 protected:

  //- Start model initialization after parent class
  const Switch startModelInit_;

  //- Energy Model
  const simpleEnergyModel& energyModel_;

  //- Recession [m]
  const volScalarField& recession;

  //- Solid density [kg/m3]
  const volScalarField& rho_s;

  //- Char mass flow rate at the walls [kg/m2/s]
  const volScalarField& mDotCw;

  //- Reference to Mass Model
  const simpleMassModel& MassModel;

  //- Pyrolysis gas mass flow rate at the walls [kg/m2/s]
  const surfaceVectorField& mDotGFace;

  //- Reference to Pyrolysis Model
  const simplePyrolysisModel& PyrolysisModel;

  //- Virgin density [kg/m3]
  const volScalarField& rho_v;

  //- Char density [kg/m3]
  const volScalarField& rho_c;

  //- Name of the moving patch
  const word topPatchName;

  //- ID of the moving patch
  const label topPatchID;

  //- Mass output file name
  const fileName massOutputFile;

  //- Mass output stream
  OFstream os_mass;

  //- Mass loss output file name
  const fileName massLossOutputFile;

  //- Mass loss output stream
  OFstream os_massLoss;

  //- Initial density [kg/m3]
  scalar rho_s0;

  //- Initial volume [m3]
  scalar V0;

  //- Initial position of top patch [m]
  vector initialPosition;

 public:

  //- Runtime type information
  TypeName("mass");

  // Constructors

  //- Construct from mesh and phase name
  massIOModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~massIOModel();

  //- Member Functions

  //- Update
  virtual void update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //

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
    Foam::GradPTimeControlModel

Description
    Time control model based on the gradient of pressure

SourceFiles
    GradPTimeControlModel.H
    GradPTimeControlModel.C

\*---------------------------------------------------------------------------*/

#ifndef GradPTimeControlModel_H
#define GradPTimeControlModel_H

#include "simpleTimeControlModel.H"
#include "simpleMaterialsModel.H"
#include "simpleMassModel.H"
#include "simpleEnergyModel.H"
#include "IOFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class GradPTimeControlModel Declaration
\*---------------------------------------------------------------------------*/

class GradPTimeControlModel
:
public simpleTimeControlModel
{

 protected:

  //- Energy model
  const simpleEnergyModel& energyModel_;

  //- Solid temperature [K]
  const volScalarField& T_;

  //- Mass model
  const simpleMassModel& massModel_;

  //- gas velocity [m/s]
  const volVectorField& vG_;

  //- Interpolation on the faces of the gas velocity [m/s]
  surfaceVectorField& vS_;

  //- Adjust time step flag
  Switch adjustTimeStep;

  //- Interface name
  word stitchingSolidPatchName;

  //- Maximum current number
  scalar maxCo;

  //- adjust start time
  scalar adjustStartTime;

  //- Maximum time step
  scalar maxDeltaT;

  //- Minimum time step
  scalar minDeltaT;

  //- REV length
  scalar REVlength;

  //- Fluid Stitch criteria  //- NAmrofel 26/11/19
  //scalar& MaxDeltaTw_;

  //- Solid Stitch criteria  //- NAmrofel 26/11/19
  scalar minDeltaTw_;

  scalarField Tws_n;
  scalarField Tws_nless1;
  scalarField deltaTws;

 public:

  //- Runtime type information
  TypeName("GradP");

  // Constructors

  //- Construct from mesh and phase name
  GradPTimeControlModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~GradPTimeControlModel();

  //- Member Functions

  //- Update thermodynamic and transport properties in the field
  virtual void update();

  //- Update stitch criteria
  virtual scalar updateMinDeltaTw();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::GradP_ChemYEqnTimeControlModel

Description
    Time control model based on the gradient of pressure and on the gradient of species mass fraction

SourceFiles
    GradP_ChemYEqnTimeControlModel.H
    GradP_ChemYEqnTimeControlModel.C

\*---------------------------------------------------------------------------*/

#ifndef GradP_ChemYEqnTimeControlModel_H
#define GradP_ChemYEqnTimeControlModel_H

#include "simpleTimeControlModel.H"
#include "simpleMassModel.H"
#include "simpleMaterialChemistryModel.H"
#include "IOFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class GradP_ChemYEqnTimeControlModel Declaration
\*---------------------------------------------------------------------------*/

class GradP_ChemYEqnTimeControlModel
:
public simpleTimeControlModel
{

 protected:

  //- Reference mesh
  const fvMesh& mesh_;

  //- Dictionary name (region name)
  const word& dictName_;

  //- Mass model
  simpleMassModel& massModel_;

  //- gas velocity [m/s]
  volVectorField& vG_;

  //- Interpolation on the faces of the gas velocity [m/s]
  surfaceVectorField vS_;

  //- MaterialChemistry model
  simpleMaterialChemistryModel& MaterialChemistryModel_;

  //- Species mass fractions
  const PtrList<volScalarField>& Y_;

  //- Species mass fractions
  const PtrList<volScalarField>& Yold_;

  //- Species names
  const wordList speciesNames_;

  //- Latest estimation of MaterialChemistry integration
  const scalar& dtChem_;

  //- Solid carbon index
  int iCs_;

  //- Adjust time step flag
  Switch adjustTimeStep;

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

  //- Time step: threshold on mass fraction variation rate
  scalar Ythreshold;

  //- Time step: minimum tolerance on mass fraction variation rate
  scalar dYtolMin;

  //- Time step: maximum tolerance on mass fraction variation rate
  scalar dYtolMax;

  //- New time step
  scalar dtChemYEqn;

  //- Minimum MaterialChemistry time step
  scalar minChemDeltaT;

  //- Flag to retain the min chemical time step bewteen MaterialChemistry (dtChem) and transport (dtChemYEqn)
  Switch chemTransEulerStepLimiter;

 public:

  //- Runtime type information
  TypeName("GradP_ChemYEqn");

  // Constructors

  //- Construct from mesh and phase name
  GradP_ChemYEqnTimeControlModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~GradP_ChemYEqnTimeControlModel();

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

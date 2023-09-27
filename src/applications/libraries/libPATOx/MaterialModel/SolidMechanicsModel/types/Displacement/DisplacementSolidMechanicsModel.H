/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2021 PATO team
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
    Foam::DisplacementSolidMechanicsModel

Description
    Finite-volume formulation for small strain elastic thermal pyrolisable
    solid bodies with mesh motion. Based on the work of Philip Cardiff for
    foam-extend-4.1.

SourceFiles
    DisplacementSolidMechanicsModel.H
    DisplacementSolidMechanicsModel.C

\*---------------------------------------------------------------------------*/

#ifndef DisplacementSolidMechanicsModel_H
#define DisplacementSolidMechanicsModel_H

#include "simpleMaterialPropertiesModel.H"
#include "simplePyrolysisModel.H"
#include "simpleSolidMechanicsModel.H"
#include "simpleEnergyModel.H"
#include "dynamicFvMesh.H"
#include "volPointInterpolation.H"
#include "symmetryPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class DisplacementSolidMechanicsModel Declaration
\*---------------------------------------------------------------------------*/

class DisplacementSolidMechanicsModel
:
public simpleSolidMechanicsModel
{

 protected:

  //- Displacement [m]
  volVectorField& D;

  //- Displacement increment [m]
  volVectorField& deltaD;

  //- Gradient of displacement
  volTensorField& gradD;

  //- Poisson's ratio [-]
  volScalarField& nu;

  //- Young's modulus of elasticity [Pa]
  volScalarField& E;

  //- Thermal expansion coefficient [1/K]
  volScalarField& alpha;

  //- Jomaa Function (pyrolysis shrinkage) [-]
  volScalarField& xi;

  //- Lame coefficient [Pa]
  volScalarField& mu_sM;

  //- Lame coefficient [Pa]
  volScalarField& lambda_sM;

  //- threeK [Pa]
  volScalarField& threeK;

  //- threeKalpha [Pa/K]
  volScalarField& threeKalpha;

  //- threeKxi [Pa]
  volScalarField& threeKxi;

  //- sigma [Pa]
  volSymmTensorField& sigma;

  //- divSigmaExp []
  volVectorField& divSigmaExp;

  //- Von Mises stress [Pa]
  volScalarField& sigmaEq;

  //- Strain []
  volSymmTensorField& epsilon;

  //- Von Mises strain []
  volScalarField& epsilonEq;

  //- Position of the cell centers [m]
  volVectorField& cellCenter;

  //- Position of the cell centers in x direction [m]
  volScalarField& cellCenterX;

  //- Position of the cell centers in y direction [m]
  volScalarField& cellCenterY;

  //- Position of the cell centers in z direction [m]
  volScalarField& cellCenterZ;

  //- Displacement of the cell centers in x direction [m]
  volScalarField& DisplacementX;

  //- Displacement of the cell centers in y direction [m]
  volScalarField& DisplacementY;

  //- Displacement of the cell centers in z direction [m]
  volScalarField& DisplacementZ;

  //- Solid displacement control dictionary
  const dictionary& solidControl_;

  //- Number of corrector
  int nCorr_;

  //- Tolerance for convergence criterion
  scalar convergenceTolerance_;

  //- Plane stress
  const Switch& planeStress_;

  //- PyrolysisModel
  simplePyrolysisModel& pyrolysisModel_;

  //- Advancement of the pyrolysis reactions: 1 to 0 [-]
  volScalarField& tau_;

  //- Energy Model
  simpleEnergyModel& energyModel_;

  //- Solid phase density [kg/m3]
  volScalarField& rho;

  //- Solid temperature [K]
  volScalarField& T;

  //- Old (t=t-1) solid temperature [K]
  volScalarField T0;


 public:

  //- Runtime type information
  TypeName("Displacement");

  // Constructors

  //- Construct from mesh and phase name
  DisplacementSolidMechanicsModel
  (
      const fvMesh& mesh,
      const word& dictName
  );

  //- Destructor
  virtual ~DisplacementSolidMechanicsModel();

  //- Member Functions

  //- Update the mass conservation fields
  virtual void update();

  //- Update fields before to solve the mass conservation equation
  void beforeSolve();

  //- Update fields after to solve the mass conservation equation
  void afterSolve();

  //- Apply mesh motion
  void moveMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

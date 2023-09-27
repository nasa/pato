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
    Foam::ElasticThermalSolidMechanicsModel

Description
    Finite-volume formulation for small strain elastic thermal pyrolisable
    solid bodies. Based on the work of Philip Cardiff for foam-extend-4.1

SourceFiles
    ElasticThermalSolidMechanicsModel.H
    ElasticThermalSolidMechanicsModel.C

\*---------------------------------------------------------------------------*/

#ifndef ElasticThermalSolidMechanicsModel_H
#define ElasticThermalSolidMechanicsModel_H

#include "simpleMaterialPropertiesModel.H"
#include "simplePyrolysisModel.H"
#include "simpleSolidMechanicsModel.H"
#include "simpleEnergyModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ElasticThermalSolidMechanicsModel Declaration
\*---------------------------------------------------------------------------*/

class ElasticThermalSolidMechanicsModel
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

  //- ThreeKalpha [Pa/K]
  volScalarField& threeKalpha;

  //- ThreeKxi [Pa]
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

  //- Solid displacement control dictionary
  const dictionary& solidControl_;

  //- Number of corrector
  int nCorr_;

  //- Tolerance for convergence criterion
  scalar convergenceTolerance_;

  //- Plane stress
  const Switch& planeStress_;

  //- PyrolysisModel
  const simplePyrolysisModel& pyrolysisModel_;

  //- Advancement of the pyrolysis reactions: 1 to 0 [-]
  const volScalarField& tau_;

  //- Energy Model
  const simpleEnergyModel& energyModel_;

  //- Solid phase density [kg/m3]
  const volScalarField& rho;

  //- Solid temperature [K]
  const volScalarField& T;

  //- Initial (t=0) Solid temperature [K]
  const volScalarField& T0;

  //- Initialize switch
  Switch initialized_;

 public:

  //- Runtime type information
  TypeName("ElasticThermal");

  // Constructors

  //- Construct from mesh and phase name
  ElasticThermalSolidMechanicsModel
  (
      const fvMesh& mesh,
      const word& regionName
  );

  //- Destructor
  virtual ~ElasticThermalSolidMechanicsModel();

  //- Member Functions

  //- Initialize
  void initialize();

  //- Update the mass conservation fields
  virtual void update();

  //- Update fields before to solve the mass conservation equation
  void beforeSolve();

  //- Update fields after to solve the mass conservation equation
  void afterSolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

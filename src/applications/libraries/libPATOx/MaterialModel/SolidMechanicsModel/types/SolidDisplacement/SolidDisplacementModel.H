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
    Foam::SolidDisplacementModel

Description
    Semi-implicit formulation of mass conservation (with Darcy's law -ie. ,
    momentum conservation- substituted inside the mass conservation equation)

SourceFiles
    SolidDisplacementModel.H
    SolidDisplacementModel.C

\*---------------------------------------------------------------------------*/

#ifndef SolidDisplacementModel_H
#define SolidDisplacementModel_H

#include "simpleMassModel.H"
#include "simpleGasPropertiesModel.H"
#include "simpleMaterialPropertiesModel.H"
#include "simplePyrolysisModel.H"
#include "simpleSolidMechanicsModel.H"

#include "dynamicFvMesh.H"
#include "volPointInterpolation.H"
#include "symmetryPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class SolidDisplacementModel Declaration
\*---------------------------------------------------------------------------*/

class SolidDisplacementModel
:
public simpleSolidMechanicsModel
{

 protected:

  //- Mesh reference
  const fvMesh& mesh_;

  //- Dictionary name
  const word& dictName_;

  //- GasPropertiesModel
  simpleGasPropertiesModel& gasPropertiesModel_;

  //- MaterialPropertiesModel
  simpleMaterialPropertiesModel& materialPropertiesModel_;

  //- PyrolysisModel
  simplePyrolysisModel& pyrolysisModel_;

  //- Gaseous pressure [Pa]
  volScalarField& p;

  //- Gaseous temperature [K]
  volScalarField& T;

  //- Gaseous molar mass [kg/kmol]
  volScalarField& M;

  //- Gaseous viscosity [kg/m/s]
  volScalarField& mu;

  //- Gaseous volume fraction [-]
  volScalarField& eps_g;

  //- Gaseous mass density [-]
  volScalarField& rho_g;

  //- Universal gas constant [J/mol/K]
  const dimensionedScalar R;

  //- Solid permeability [m2]
  volTensorField& K;

  //- Total gas pyrolysis production [kg/m3/s]
  volScalarField& piTotal;

  //- Part of the convective mass flux [s] : q_conv = dx . ( Gamma dx p )
  volTensorField Gamma;

  //- Symmetrical of the part of the convective mass flux [s] : q_conv = dx . ( Gamma dx p )
  volSymmTensorField Gamma_symm;

  //- Part of the storage mass flux [s] : q_sto = dt( Beta p )
  volScalarField Beta;

  //- Gaseous velocity times gas volume fraction [m/s]
  volVectorField& U;

  /*****/

  //- Displacement [m]
  volVectorField& D;

  //- Poisson's ratio [-]
  volScalarField& nu;

  //- Young's modulus of elasticity [Pa]
  volScalarField& E;

  //- Lamé's coefficient [Pa]
  volScalarField mu_sM;

  //- Lamé's coefficient [Pa]
  volScalarField lambda_sM;

  //- threeK [Pa]
  volScalarField threeK;

  //- Thermal expansion coefficient [1/K]
  volScalarField& alpha;

  //- ThreeKalpha [Pa/K]
  volScalarField threeKalpha;

  //- Advancement of the pyrolysis reactions: 1 to 0 [-]
  volScalarField& tau_;

  //- Solid phase density [kg/m3]
  volScalarField& rho;

  //- sigmaD []
  volSymmTensorField sigmaD;

  //- sigma [Pa]
  volSymmTensorField sigma;

  //- divSigmaExp []
  volVectorField divSigmaExp;

  //- Von Mises stress [Pa]
  volScalarField sigmaEq;

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


 public:

  //- Runtime type information
  TypeName("SolidDisplacement");
  // TypeName("DarcyLaw");

  // Constructors

  //- Construct from mesh and phase name
  SolidDisplacementModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~SolidDisplacementModel();

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

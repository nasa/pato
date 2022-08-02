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
    Foam::basicUserThermo

Description

SourceFiles
    basicUserThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicUserThermo_H
#define basicUserThermo_H

#include "fvCFD.H"
#undef Log
#include "mutation++/mutation++.h"
#include "compressibleTransportModel.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     class basicUserThermo Declaration
\*---------------------------------------------------------------------------*/

class basicUserThermo
:
public compressibleTransportModel
{

 private:

  // Private Member Functions

 protected:

  // Protected data

  //- Mesh database
  const fvMesh& mesh_;

  //- Reference variables for the current state -- updated whenever needed
  scalar currentP_;
  scalar currentE_;
  scalar currentT_;
  scalar currentRho_;
  scalar currentPsi_;
  scalar currentMu_;
  scalar currentNu_;
  scalar currentSigma_;
  scalar currentGamma_;
  scalar currentKappa_;
  scalar currentAlpha_;
  scalarList currentYe_;
  scalarList currentYs_;
  scalarList currentXe_;
  scalarList currentXs_;
  scalarList currentRhoe_;
  scalarList currentRhos_;
  scalarList currentFluxesP_;
  scalarList currentFluxesT_;
  scalarList currentFluxesX_;
  scalarList currentOmegas_;

  //- Name of mixture (default is air13_AHF)
  word mixtureName_;

  //- Name of the state model (Equil or ChemNonEq1T)
  word stateModelName_;

  //- Pointer to Mutation::Mixture, used for Mutation++ functions
  std::unique_ptr<Mutation::Mixture> pMutationMix_;

  //- List of species
  wordList speciesList_;

  //- List of elements
  wordList elementsList_;

  //- Number of species
  int nSpecies_;

  //- Number of elements
  int nElements_;

  //- Number of mass equations
  int nMassEqns_;

  //- Number of energy equations
  int nEnergyEqns_;

  //- List of species molar masses [kg/mol]
  scalarList speciesMolarMassList_;

  //- Universal gas constant RU
  scalar RU_ = Mutation::RU;

  //- Presure [Pa]
  volScalarField& p_;

  //- Temperature [K]
  volScalarField& T_;

  //- Mixture total internal energy [m2/s2]
  volScalarField e_;

  //- Mixture mass density [kg/m3]
  volScalarField& rho_;

  //- Compressibility (rho/p) [s2/m2]
  volScalarField& psi_;

  //- Electrical conductivity [?]
  volScalarField& sigma_;

  //- Isentropic exponent (Cp/Cv ratio) [-]
  volScalarField& gamma_;

  //- Thermal diffusivity [kg/m/s]
  volScalarField& alpha_;

  //- Dynamic viscosity [kg/m/s]
  volScalarField& mu_;

  //- Kinematic viscosity [m2/s]
  volScalarField& nu_;

  //- Thermal conductivity [kg*m/s3/K]
  volScalarField& kappa_;

  //- Effective thermal conductivity (including turbulent) [kg*m/s3/K]
  volScalarField& kappaEff_;

  //- Effective thermal diffusivity (including turbulent) [kg/m/s]
  volScalarField& alphaEff_; // see heThermo.H

  //- Thermal diffusivity for energy of mixture [kg/m/s]
  volScalarField& alphahe_;

  //- Element mass fractions [-]
  PtrList<volScalarField> Ye_;

  //- Element mole fractions [-]
  PtrList<volScalarField> Xe_;

  //- Element mass densities [kg/m3]
  PtrList<volScalarField> rhoe_;

  //- Species mass fractions [-]
  PtrList<volScalarField> Ys_;

  //- Species mole fractions [-]
  PtrList<volScalarField> Xs_;

  //- Species mass densities [kg/m3]
  PtrList<volScalarField> rhos_;

  //- Element mass flux per pressure gradient [s]
  PtrList<volScalarField> elementsMassDiffFluxesGradP_;

  //- Element mass flux per temperature gradient [kg/m.s.K]
  PtrList<volScalarField> elementsMassDiffFluxesGradT_;

  //- Element mass flux per mole fractions gradient [kg/m.s]
  PtrList<volScalarField> elementsMassDiffFluxesGradX_;

  //- Element energy flux per pressure gradient [m2/s]
  autoPtr<volScalarField> elementsEnergyDiffFluxesGradP_;

  //- Element energy flux per temperature gradient [kg.m/s3.K]
  autoPtr<volScalarField> elementsEnergyDiffFluxesGradT_;

  //- Element energy flux per mole fractions gradient [kg.m/s3]
  PtrList<volScalarField> elementsEnergyDiffFluxesGradX_;

  //- Elements diffusion velocities [m/s]
  PtrList<volVectorField> diffMassFluxes_;

  //- Species diffusion velocities [m/s]
  PtrList<volVectorField> diffVs_;

  //- Species net mass production rates [kg/m3/s]
  PtrList<volScalarField> omegas_;

  //- Diffusive heat flux vector [kg/s3]
  autoPtr<volVectorField> qDiff_;



  //- Mixture pointer configuration
  //- ATM: would be better to find a way to add it in the constructor
  virtual void configureMixture();

  //- Creates composition fields in PtrList<volScalarField>
  virtual void createCompositionFields(const fvMesh& mesh);

  //- Creates elemental flux fields in PtrList<volScalarField>
  virtual void createElementFluxFields(const fvMesh& mesh);

  //- Creates finite rate fields in PtrList<volScalarField>
  virtual void createFiniteRateFields(const fvMesh& mesh);

 public:

  //- Boolean to define when the code is in the constructor (differentiates setState)
  bool isInConstructor_;

  //- Runtime type information
  TypeName("basicUserThermo");


  //- Declare run-time constructor selection tables
  declareRunTimeSelectionTable
  (
      autoPtr,
      basicUserThermo,
      fvMesh,
      (const fvMesh& mesh,const word& dictName),
      (mesh, dictName)
  );


  // Constructors

  //- Construct from mesh and phase name
  basicUserThermo
  (
      const fvMesh& mesh,
      const word& dictName
  );


  //- Selector
  static autoPtr<basicUserThermo> New
  (
      const fvMesh& mesh,
      const word& dictName
  );

  //- Destructor
  virtual ~basicUserThermo();

  //- Loop over all cells and update all thermodynamic properties
  virtual void update() = 0;


  // Access to thermodynamic/transport/chemistry state variables

  //- Return the mixture name
  virtual const word& getMixtureName() const;

  //- Return the state model name
  virtual const word& getStateModelName() const;

  //- Return the number of species
  virtual const int& nSpecies() const;

  //- Return the number of elements
  virtual const int& nElements() const;

  //- Return the number of mass equations
  virtual const int& nMassEqns() const;

  //- Return the number of energy equations
  virtual const int& nEnergyEqns() const;

  //- Return the elements name list
  virtual const wordList& getElementsList() const;

  //- Returne an element name
  virtual word& getElementName(const int ie);

  //- Return the species name list
  virtual const wordList& getSpeciesList() const;

  //- Return a species name
  virtual word& getSpeciesName(const int is);

  //- Return the species molar masses list
  virtual const scalarList& speciesMolarMassList() const;

  //- Return a species molar mass
  virtual scalar& speciesMolarMass(const int is);

  //- Return the universal gas constant RU
  virtual scalar& RU();

  //- Pressure [Pa]
  virtual volScalarField& p();

  //- Internal energy [K]
  virtual volScalarField& e();

  //- Temperature [K]
  virtual volScalarField& T();

  //- Density [kg/m3]
  virtual const volScalarField& rho() const ;

  //- Density [kg/m3]
  virtual volScalarField& rho();

  //- Compressibility (rho/p) [s2/m2]
  virtual const volScalarField& psi() const ;

  //- Electrical conductivity [S/m]
  virtual volScalarField& sigma();

  //- Heat capacity ratio (cp/cv) [-]
  virtual const volScalarField& gamma() const ;

  //- Elements mass fractions [-]
  virtual PtrList<volScalarField>& elementsMassFractions();

  //- Elements mole fractions [-]
  virtual PtrList<volScalarField>& elementsMoleFractions();

  //- Elements density [kg/m3]
  virtual PtrList<volScalarField>& elementsDensities();

  //- Species mass fractions [-]
  virtual PtrList<volScalarField>& speciesMassFractions();

  //- Species mole fractions [-]
  virtual PtrList<volScalarField>& speciesMoleFractions();

  //- Species density [kg/m3]
  virtual PtrList<volScalarField>& speciesDensities();

  //- Elements mass fluxes per pressure gradient [s]
  virtual PtrList<volScalarField>& elementsMassDiffFluxesGradP();

  //- Elements mass fluxes per temperature gradient [kg/m.s.K]
  virtual PtrList<volScalarField>& elementsMassDiffFluxesGradT();

  //- Elements mass fluxes per mole fraction gradient [kg/m.s]
  virtual PtrList<volScalarField>& elementsMassDiffFluxesGradX();

  //- Elements energy fluxes per pressure gradient [m2/s]
  virtual autoPtr<volScalarField>& elementsEnergyDiffFluxesGradP();

  //- Elements energy fluxes per temperature gradient [kg.m/s3.K]
  virtual autoPtr<volScalarField>& elementsEnergyDiffFluxesGradT();

  //- Elements energy fluxes per mole fraction gradient [kg.m/s3]
  virtual PtrList<volScalarField>& elementsEnergyDiffFluxesGradX();

  //- Elements diffusion velocities [m/s]
  virtual PtrList<volVectorField>& elementsDiffusionMassFluxes();

  //- Species diffusion velocities [m/s]
  virtual PtrList<volVectorField>& speciesDiffusionVelocities();

  //- Species net mass production rates [kg/m3/s]
  virtual PtrList<volScalarField>& speciesNetMassProductionRates();

  //- Diffusive heat flux [kg/s3]
  virtual autoPtr<volVectorField>& diffHeatFlux();


  // Needed by turbulence models

  //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
  virtual const volScalarField& alpha() const ;

  //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
  virtual const scalarField& alpha
  (
      const label patchi
  ) const;

  //- Dynamic viscosity of mixture [kg/m/s]
  virtual tmp<volScalarField> mu() const;

  //- Dynamic viscosity of mixture for patch [kg/m/s]
  virtual tmp<scalarField> mu(const label patchi) const;

  //- Kinematic viscosity of mixture [m^2/s]
  virtual tmp<volScalarField> nu() const;

  //- Kinematic viscosity of mixture for patch [m^2/s]
  virtual tmp<scalarField> nu(const label patchi) const ;

  //- Thermal diffusivity for temperature of mixture [J/m/s/K]
  virtual tmp<volScalarField> kappa() const ;

  //- Thermal diffusivity for temperature
  //  of mixture for patch [J/m/s/K]
  virtual tmp<scalarField> kappa
  (
      const label patchi
  ) const ;

  //- Thermal diffusivity for energy of mixture [kg/m/s]
  virtual tmp<volScalarField> alphahe() const;

  //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
  virtual tmp<scalarField> alphahe(const label patchi) const;

  //- Effective thermal diffusivity for temperature
  //  of mixture [J/m/s/K]
  virtual tmp<volScalarField> kappaEff
  (
      const volScalarField& kappat
  ) const ;

  //- Effective thermal diffusivity for temperature
  //  of mixture for patch [J/m/s/K]
  virtual tmp<scalarField> kappaEff
  (
      const scalarField& kappat,
      const label patchi
  ) const ;

  //- Effective thermal diffusivity of mixture [kg/m/s]
  virtual tmp<volScalarField> alphaEff
  (
      const volScalarField& alphat
  ) const ;

  //- Effective thermal diffusivity of mixture for patch [kg/m/s]
  virtual tmp<scalarField> alphaEff
  (
      const scalarField& alphat,
      const label patchi
  ) const ;


  //- Used to initialize/read the thermodynamic variables

  //- Creates a scalar field with MUST_READ if not found in the mesh database
  inline volScalarField& meshLookupOrConstructScalar
  (
      const fvMesh& mesh,
      const char* name
  );

  //- Creates a scalar field with READ_IF_PRESENT if not found in the mesh database
  inline volScalarField& meshLookupOrConstructScalar
  (
      const fvMesh& mesh,
      const char* name,
      const dimensionedScalar& dimScalar_
  );

  //- Creates a scalar field with READ_IF_PRESENT and a given boundary condition for all patches if not found in the mesh database
  inline volScalarField& meshLookupOrConstructScalar
  (
      const fvMesh& mesh,
      const char* name,
      const dimensionedScalar& dimScalar_,
      const wordList& BC_
  );

  //- Creates a vector field with MUST_READ if not found in the mesh database
  inline volVectorField& meshLookupOrConstructVector
  (
      const fvMesh& mesh,
      const char* name
  );

  //- Creates a vector field  with READ_IF_PRESENT for all patches if not found in the mesh database
  inline volVectorField& meshLookupOrConstructVector
  (
      const fvMesh& mesh,
      const char* name,
      const dimensionedVector& dimVec_
  );

  //- Creates a vector field with READ_IF_PRESENT and a given boundary condition for all patches if not found in the mesh database
  inline volVectorField& meshLookupOrConstructVector
  (
      const fvMesh& mesh,
      const char* name,
      const dimensionedVector& dimVector_,
      const wordList& BC_
  );

  //- Creates a tensor field if not found in the mesh database
  inline volTensorField& meshLookupOrConstructTensor
  (
      const fvMesh& mesh,
      const char* name
  );

  //- Creates a tensor field with given a boundary condition for all patches if not found in the mesh database
  inline volTensorField& meshLookupOrConstructTensor
  (
      const fvMesh& mesh,
      const char* name,
      const dimensionedTensor& dimTens_
  );

  //- Creates a tensor field with READ_IF_PRESENT and a given boundary condition for all patches if not found in the mesh database
  inline volTensorField& meshLookupOrConstructTensor
  (
      const fvMesh& mesh,
      const char* name,
      const dimensionedTensor& dimTensor_,
      const wordList& BC_
  );

  //- Return the internal energy field boundary types
  //  by interrogating the temperature field boundary types
  inline wordList eBoundaryTypes();
};

#include "basicUserThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 PATO-community
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    Based on tilasoldo, Yuusha and chriss85 contribution to OpenFOAM, this new
    thermophysical model has been modified and checked by PATO-community.

    The Interpolation function used in this updated file is that of OpenFoam
    called "interpolation2DTable".

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
    Foam::tabularThermo

Description
    Basic thermodynamic properties based on compressibility

SourceFiles
    tabularThermo.C

\*---------------------------------------------------------------------------*/

#ifndef tabularThermo_H
#define tabularThermo_H

#include "fluidThermo.H"
#include "runTimeSelectionTables.H"
#include "interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class tabularThermo Declaration
\*---------------------------------------------------------------------------*/

class tabularThermo
:
public fluidThermo
{

 protected:

  // Protected data

  //- Compressibility [s^2/m^2]
  volScalarField psi_;

  //- Dynamic viscosity [kg/m/s]
  volScalarField mu_;

  //- Density [kg/m^3]
  interpolation2DTable<scalar> densityTable;

  // Protected Member Functions

  //- Construct as copy (not implemented)
  tabularThermo(const tabularThermo&);


 public:
  //- Runtime type information
  TypeName("tabularThermo");


  //- Declare run-time constructor selection table
  declareRunTimeSelectionTable
  (
      autoPtr,
      tabularThermo,
      fvMesh,
      (const fvMesh& mesh, const word& phaseName),
      (mesh, phaseName)
  );


  // Constructors
  //- Construct from mesh and phase name
  tabularThermo
  (
      const fvMesh&,
      const word& phaseName
  );


  //- Selector
  static autoPtr<tabularThermo> New
  (
      const fvMesh& mesh,
      const word& phaseName=word::null
  );


  //- Destructor
  virtual ~tabularThermo();


  // Member functions

  volScalarField& lookupOrConstruct2
  (
      const fvMesh& mesh,
      const char* name,
      dimensionSet units
  ) const;

  // Fields derived from thermodynamic state variables

  //- Add the given density correction to the density field.
  //  Used to update the density field following pressure solution.
  //  For psiThermo does nothing.
  virtual void correctRho(const volScalarField& deltaRho);

  //- Density [kg/m^3] - uses current value of pressure
  virtual tmp<volScalarField> rho() const;

  //- Compressibility [s^2/m^2]
  virtual const volScalarField& psi() const;


  // Access to transport state variables

  //- Dynamic viscosity of mixture [kg/m/s]
  virtual tmp<volScalarField> mu() const;

  //- Dynamic viscosity of mixture for patch [kg/m/s]
  virtual tmp<scalarField> mu(const label patchi) const; //= 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

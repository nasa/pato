/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2011-2019 OpenFOAM Foundation
                            | Copyright (C) 2020      PATO
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::darcyFlowRatePressureFvPatchScalarField

Description
    Pressure condition for Darcy flow calculated for a uniform velocity field
    normal to the patch adjusted to match the specified flow rate

    For a mass-based flux:
    - the flow rate should be provided in kg/s
    - if \c rho is "none" the flow rate is in m3/s
    - otherwise \c rho should correspond to the name of the density field
    - if the density field cannot be found in the database, the user must
      specify the inlet density using the \c rhoInlet entry

    For a volumetric-based flux:
    - the flow rate is in m3/s


Usage
    \table
        Property           | Description                 | Required | Default value
        massFlowRate       | mass flow rate [kg/s]       | no       |
        volumetricFlowRate | volumetric flow rate [m3/s] | no       |
        rho                | density field name          | no       | rho_g
        rhoInlet           | inlet density               | no       |
    \endtable

    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    <patchName>
    {
        type                darcyFlowRatePressure;
        volumetricFlowRate  0.2;
        value               uniform 1e5;
    }
    \endverbatim

    Example of the boundary condition specification for a mass flow rate:
     \verbatim
    <patchName>
    {
        type                darcyFlowRatePressure;
        massFlowRate        0.2;
        rho                 rho_g;
        rhoInlet            1.0;
        value               uniform 1e5;
    }
    \endverbatim

    The \c flowRate entry is a \c Function1 of time, see Foam::Function1Types.

Note
    - \c rhoInlet is required for the case of a mass flow rate, where the
      density field is not available at start-up
    - The value is positive into the domain (as an inlet)

See also
    Foam::darcyVelocitypressureFvPatchScalarField
    Foam::Function1Types

SourceFiles
    darcyFlowRatePressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef darcyFlowRatePressureFvPatchScalarFields_H
#define darcyFlowRatePressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class darcyFlowRatePressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class darcyFlowRatePressureFvPatchScalarField
:
public fixedGradientFvPatchScalarField
{
  // Private data

  //- Darcy's velocity field
  vectorField U_;

  //- Inlet integral flow rate
  autoPtr<Function1<scalar>> flowRate_;

  //- Is volumetric?
  bool volumetric_;

  //- Name of the density field used to normalize the mass flux
  word rhoName_;

  //- Rho initialisation value (for start; if value not supplied)
  scalar rhoInlet_;

  // Private Member Functions

  //- Update the patch values given the appropriate density type and value
  template<class RhoType>
  vectorField updateValues(const RhoType& rho);

 public:

  //- Runtime type information
  TypeName("darcyFlowRatePressure");


  // Constructors

  //- Construct from patch and internal field
  darcyFlowRatePressureFvPatchScalarField
  (
      const fvPatch&,
      const DimensionedField<scalar, volMesh>&
  );

  //- Construct from patch, internal field and dictionary
  darcyFlowRatePressureFvPatchScalarField
  (
      const fvPatch&,
      const DimensionedField<scalar, volMesh>&,
      const dictionary&
  );

  //- Construct by mapping given darcyFlowRatePressureFvPatchScalarField onto
  //  a new patch
  darcyFlowRatePressureFvPatchScalarField
  (
      const darcyFlowRatePressureFvPatchScalarField&,
      const fvPatch&,
      const DimensionedField<scalar, volMesh>&,
      const fvPatchFieldMapper&
  );

  //- Copy constructor
  darcyFlowRatePressureFvPatchScalarField
  (
      const darcyFlowRatePressureFvPatchScalarField&
  );

  //- Construct and return a clone
  virtual tmp<fvPatchScalarField> clone() const {
    return tmp<fvPatchScalarField>
           (
               new darcyFlowRatePressureFvPatchScalarField(*this)
           );
  }

  //- Copy constructor setting internal field reference
  darcyFlowRatePressureFvPatchScalarField
  (
      const darcyFlowRatePressureFvPatchScalarField&,
      const DimensionedField<scalar, volMesh>&
  );

  //- Construct and return a clone setting internal field reference
  virtual tmp<fvPatchScalarField> clone
  (
      const DimensionedField<scalar, volMesh>& iF
  ) const {
    return tmp<fvPatchScalarField>
           (
               new darcyFlowRatePressureFvPatchScalarField(*this, iF)
           );
  }


  // Member functions

  //- Update the coefficients associated with the patch field
  virtual void updateCoeffs();

  //- Write
  virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

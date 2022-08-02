/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::
    compressible::
    vectorCoupledMixedFvPatchVectorField

Description
    Mixed boundary condition for pressure, to be used for heat-transfer
    on back-to-back baffles.Optional thin thermal layer
    resistances can be specified through thicknessLayers and kappaLayers
    entries.

    The thermal conductivity, \c kappa, can either be retrieved from the
    mesh database using the \c lookup option, or from a \c solidThermo
    or \c fluidThermo thermophysical package.


    Specifies gradient and pressure such that the equations are the same
    on both sides:
    - refGradient = zero gradient
    - refValue = neighbour value
    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

    where KDelta is heat-transfer coefficient K * deltaCoeffs

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        kappa        | thermal conductivity option | yes     |
        kappaName    | name of thermal conductivity field | no  |
        nbr         | name of the field    | no | T
        thicknessLayers | list of thicknesses per layer [m] | no |
        kappaLayers  | list of thermal conductivites per layer [W/m/K] | no |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            compressible::coupledMixed;
        nbr            T;
        kappa           lookup;
        KappaName       kappa;
        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4)
        value           uniform 300;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

    Note: kappa : heat conduction at patch. Gets supplied how to lookup
        calculate kappa:

    - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
    - 'fluidThermo' : use fluidThermo and compressible::RASmodel to calculate
        kappa
    - 'solidThermo' : use solidThermo kappa()
    - 'directionalSolidThermo' directionalKappa()

SourceFiles
    vectorCoupledMixedFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValueToNbrValueFvPatchField_H
#define fixedValueToNbrValueFvPatchField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class vectorCoupledMixedFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedValueToNbrValueFvPatchField
:
public mixedFvPatchField<Type>
//    public pressureCoupledBase
{
  // Private data

  //- Name of field on the neighbour region
  const word nbrName_;

 public:

  //- Runtime type information
  TypeName("fixedValueToNbrValue");


  // Constructors

  //- Construct from patch and internal field
  fixedValueToNbrValueFvPatchField
  (
      const fvPatch&,
      const DimensionedField<Type, volMesh>&
  );

  //- Construct from patch, internal field and dictionary
  fixedValueToNbrValueFvPatchField
  (
      const fvPatch&,
      const DimensionedField<Type, volMesh>&,
      const dictionary&
  );

  //- Construct by mapping given
  //  fixedValueToNbrValuevectorCoupledMixedFvPatchVectorField onto a
  //  new patch
  fixedValueToNbrValueFvPatchField
  (
      const fixedValueToNbrValueFvPatchField<Type>&,
      const fvPatch&,
      const DimensionedField<Type, volMesh>&,
      const fvPatchFieldMapper&
  );

  //- Construct and return a clone
  virtual tmp<fvPatchField<Type> > clone() const {
    return tmp<fvPatchField<Type> >
           (
               new fixedValueToNbrValueFvPatchField<Type>
               (
                   *this
               )
           );
  }

  //- Construct as copy setting internal field reference
  fixedValueToNbrValueFvPatchField
  (
      const fixedValueToNbrValueFvPatchField<Type>&,
      const DimensionedField<Type, volMesh>&
  );

  //- Construct and return a clone setting internal field reference
  virtual tmp<fvPatchField<Type> > clone
  (
      const DimensionedField<Type, volMesh>& iF
  ) const {
    return tmp<fvPatchField<Type> >
           (
               new fixedValueToNbrValueFvPatchField<Type>
               (
                   *this,
                   iF
               )
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

#ifdef NoRepository
#include "fixedValueToNbrValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    along with OpenFOAM.  If Porous_t, see <http://www.gnu.org/licenses/>.

Class
    Foam::GradedPorousMaterialPropertiesModel

Description

SourceFiles
    GradedPorousMaterialPropertiesModel.H
    GradedPorousMaterialPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef GradedPorousMaterialPropertiesModel_H
#define GradedPorousMaterialPropertiesModel_H

#include "PorousMaterialPropertiesModel.H"
#include "simplePyrolysisModel.H"
#include "LinearInterpolationMaterialPropertiesObject.H"
#include "dimensionedTensor.H"
#include "IOFunctions.H"
#include "mathFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class GradedPorousMaterialPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class GradedPorousMaterialPropertiesModel
:
public PorousMaterialPropertiesModel
{

 protected:

  //- Start model initialization after parent class
  const Switch startModelInit_;

  //- Number of solid phases
  const int nSolidPhases_;

  //- Pyrolysis
  simplePyrolysisModel& pyrolysisModel_;

  //- Initial solid phase volume fractions [-]
  PList<volScalarField> solidEpsI_;

  //- Solid phase volume fractions [-]
  PList<volScalarField> solidEps_;

  //- Initial solid phase densities [kg/m3]
  PList<volScalarField> solidRhoI_;

  //- Solid phase densities [kg/m3]
  PList<volScalarField> solidRho_;

  //- File names of the grading properties (distant, volume fraction)
  PList<fileName> gradingFiles_;

  //- Patch names for the grading
  PList<word> gradingPatchNames_;

  //- Patch ID for the grading
  PList<label> gradingPatchID_;

 public:

  //- Runtime type information
  TypeName("GradedPorous");

  // Constructors

  //- Construct from mesh and phase name
  GradedPorousMaterialPropertiesModel
  (
      const fvMesh& mesh,
      const word& dictName
  );

  //- Destructor
  virtual ~GradedPorousMaterialPropertiesModel();

  // Member Functions

  //- Initialize solidEpsI and solidEps
  virtual void initSolidEps();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    along with OpenFOAM.  If Fourier_Radiationt, see <http://www.gnu.org/licenses/>.

Class
    Foam::Fourier_RadiationMaterialPropertiesModel

Description

SourceFiles
    Fourier_RadiationMaterialPropertiesModel.H
    Fourier_RadiationMaterialPropertiesModel.C

\*---------------------------------------------------------------------------*/

#ifndef Fourier_RadiationMaterialPropertiesModel_H
#define Fourier_RadiationMaterialPropertiesModel_H

#include "FourierMaterialPropertiesModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class Fourier_RadiationMaterialPropertiesModel Declaration
\*---------------------------------------------------------------------------*/

class Fourier_RadiationMaterialPropertiesModel
:
public FourierMaterialPropertiesModel
{

 protected:

  //- Start model initialization after parent class
  const Switch startModelInit_;

  //- Radiative heat flux [W/m2]
  volScalarField& Qr_;

  //- Emissivity [-]
  volScalarField& emissivity_;

  //- Boltzmann constant [m2 kg/s2/K]
  const dimensionedScalar& sigmaSB;

  //- Temperature in the background [K]
  volScalarField& Tbackground_;

 public:

  //- Runtime type information
  TypeName("FourierRadiation");

  // Constructors

  //- Construct from mesh and phase name
  Fourier_RadiationMaterialPropertiesModel(const fvMesh& mesh,const word& regionName);

  //- Destructor
  virtual ~Fourier_RadiationMaterialPropertiesModel();

  //- Member Functions

  //- Update
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

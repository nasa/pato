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
    along with OpenFOAM.  If FixedTemperaturet, see <http://www.gnu.org/licenses/>.

Class
    Foam::TabulatedTemperatureEnergyModel

Description

SourceFiles
    TabulatedTemperatureEnergyModel.H
    TabulatedTemperatureEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef TabulatedTemperatureEnergyModel_H
#define TabulatedTemperatureEnergyModel_H

#include "simpleEnergyModel.H"
#include "mathFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class TabulatedTemperatureEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class TabulatedTemperatureEnergyModel
:
public simpleEnergyModel
{

 protected:
  //- origin patch name
  word originPatchName_;

  //- Internal temperature table file name
  fileName intTable_fileName_;

  //- Internal temperature table data
  List<scalarList> intTable_data_;

  //- BC values for the temperature ("name BC", value_BC)
  List<Tuple2<word,scalar>> bcT_;

  //- Temperature [K]
  volScalarField& T_;

 public:

  //- Runtime type information
  TypeName("TabulatedTemperature");

  // Constructors

  //- Construct from mesh and phase name
  TabulatedTemperatureEnergyModel(const fvMesh& mesh,const word& dictName);

  //- Destructor
  virtual ~TabulatedTemperatureEnergyModel();

  //- Member Functions

  //- Update
  virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

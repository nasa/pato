/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    Unit conversion imperial/SI

Details
    BTU <=> Joule
    Rankine <=> Kelvin
    Feet <=> Meter
    Pound <=> Kilogram
    Inch <=> Meter
\*---------------------------------------------------------------------------*/

#ifndef unitConversionBritish_H
#define unitConversionBritish_H

#include "mathematicalConstants.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Conversion from BTU to Joule
inline scalar btuToJoule(const scalar btu)
{
  return (btu*1055.06);
}

//- Conversion from Joule to BTU
inline scalar jouleToBtu(const scalar joule)
{
  return (joule/1055.06);
}

//- Conversion from Rankine to Kelvin
inline scalar rankineToKelvin(const scalar rankine)
{
  return (rankine*5/9);
}

//- Conversion from Kelvin to Rankine
inline scalar kelvinToRankine(const scalar kelvin)
{
  return (kelvin*9/5);
}

//- Conversion from feet to meter
inline scalar feetToMeter(const scalar feet)
{
  return (feet/3.2808399);
}

//- Conversion from meter to feet
inline scalar meterToFeet(const scalar meter)
{
  return (meter*3.2808399);
}

//- Conversion from pound to kg
inline scalar poundToKg(const scalar pound)
{
  return (pound/2.2046226218488);
}

//- Conversion from kg to pound
inline scalar kgToPound(const scalar kg)
{
  return (kg*2.2046226218488);
}

//- Conversion from inches to meter
inline scalar inchToMeter(const scalar inch)
{
  return (inch*0.0254);
}

//- Conversion from meter to inches
inline scalar meterToInch(const scalar meter)
{
  return (meter/0.0254);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

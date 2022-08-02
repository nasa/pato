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
\*---------------------------------------------------------------------------*/

#include "rhoThermo.H"
#include "makeThermo.H"
#include "psiThermo.H"
#include "specie.H"
#include "Rhotabular.H"
#include "hTabularThermo.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"
#include "tabularTransport.H"
#include "heRhoThermo.H"
#include "pureMixture.H"
#include "hePsiThermo.H"
// added jl+cl, 28/02/2019
#include "hPolynomialThermo.H"
#include "polynomialTransport.H"
#include "perfectGas.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * Internal-energy-based * * * * * * * * * * * * * */
// NOT CHECKED
makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    tabularTransport,
    sensibleInternalEnergy,
    hTabularThermo,
    Rhotabular,
    specie
);
// THIS PACKAGE IS USED
makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    tabularTransport,
    sensibleEnthalpy,
    hTabularThermo,
    Rhotabular,
    specie
);
// NOT CHECKED
makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    tabularTransport,
    sensibleEnthalpy,
    hTabularThermo,
    Rhotabular,
    specie
);


// New thermo for PATO : polynomial + Perfect Gas, jl+cl 28/02/2019
makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    perfectGas,
    specie
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

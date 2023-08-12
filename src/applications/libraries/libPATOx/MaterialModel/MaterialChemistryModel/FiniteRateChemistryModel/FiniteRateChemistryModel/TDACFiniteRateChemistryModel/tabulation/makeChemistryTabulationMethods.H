/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#ifndef makeChemistryTabulationMethods_H
#define makeChemistryTabulationMethods_H

#include "FiniteRateChemistryTabulationMethod.H"
#include "noFiniteRateChemistryTabulation.H"
#include "ISAT_FiniteRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChemistryTabulationMethod(SS, Comp, Thermo)                        \
                                                                               \
    typedef FiniteRateChemistryTabulationMethods::SS<Comp, Thermo>                       \
        FiniteRateChemistryTabulationMethod##SS##Comp##Thermo;                           \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        FiniteRateChemistryTabulationMethod##SS##Comp##Thermo,                           \
        (#SS"<" + word(Comp::typeName_())                                      \
      + "," + Thermo::typeName() + ">").c_str(),                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    FiniteRateChemistryTabulationMethod<Comp, Thermo>::                                  \
        adddictionaryConstructorToTable                                        \
        <FiniteRateChemistryTabulationMethod##SS##Comp##Thermo>                          \
        add##FiniteRateChemistryTabulationMethods##SS##Comp##Thermo##ConstructorToTable_;


#define makeChemistryTabulationMethods(CompChemModel, Thermo)                  \
                                                                               \
    typedef FiniteRateChemistryTabulationMethod<CompChemModel, Thermo>                   \
        FiniteRateChemistryTabulationMethod##CompChemModel##Thermo;                      \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        FiniteRateChemistryTabulationMethod##CompChemModel##Thermo,                      \
        (word(FiniteRateChemistryTabulationMethod##CompChemModel##Thermo::typeName_()) + \
        '<' + word(CompChemModel::typeName_()) + "," + Thermo::typeName() + '>'\
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        FiniteRateChemistryTabulationMethod##CompChemModel##Thermo,                      \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    makeChemistryTabulationMethod                                              \
    (                                                                          \
        none,                                                                  \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryTabulationMethod                                              \
    (                                                                          \
        ISAT_FiniteRate,                                                                  \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

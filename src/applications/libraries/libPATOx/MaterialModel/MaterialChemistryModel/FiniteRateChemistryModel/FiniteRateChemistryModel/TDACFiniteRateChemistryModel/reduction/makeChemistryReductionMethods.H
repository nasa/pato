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

#ifndef makeChemistryReductionMethods_H
#define makeChemistryReductionMethods_H

#include "FiniteRateChemistryReductionMethod.H"

#include "noFiniteRateChemistryReduction.H"
#include "DAC_FiniteRate.H"
#include "DRG_FiniteRate.H"
#include "DRGEP_FiniteRate.H"
#include "EFA_FiniteRate.H"
#include "PFA_FiniteRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChemistryReductionMethod(SS, Comp, Thermo)                         \
                                                                               \
    typedef FiniteRateChemistryReductionMethods::SS<Comp, Thermo>                        \
        FiniteRateChemistryReductionMethod##SS##Comp##Thermo;                            \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        FiniteRateChemistryReductionMethod##SS##Comp##Thermo,                            \
        (#SS"<" + word(Comp::typeName_())                                      \
      + "," + Thermo::typeName() + ">").c_str(),                               \
        0                                                                      \
    );                                                                         \
                                                                               \
    FiniteRateChemistryReductionMethod<Comp, Thermo>::                                   \
        adddictionaryConstructorToTable                                        \
        <FiniteRateChemistryReductionMethod##SS##Comp##Thermo>                           \
        add##FiniteRateChemistryReductionMethods##SS##Comp##Thermo##ConstructorToTable_;


#define makeChemistryReductionMethods(CompChemModel, Thermo)                   \
                                                                               \
    typedef FiniteRateChemistryReductionMethod<CompChemModel, Thermo>                    \
        FiniteRateChemistryReductionMethod##CompChemModel##Thermo;                       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        FiniteRateChemistryReductionMethod##CompChemModel##Thermo,                       \
        (word(FiniteRateChemistryReductionMethod##CompChemModel##Thermo::typeName_()) +  \
        '<' + word(CompChemModel::typeName_()) + "," + Thermo::typeName() + '>'\
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        FiniteRateChemistryReductionMethod##CompChemModel##Thermo,                       \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethod                                               \
    (                                                                          \
        none,                                                                  \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethod                                               \
    (                                                                          \
        DAC_FiniteRate,                                                                   \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethod                                               \
    (                                                                          \
        DRG_FiniteRate,                                                                   \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethod                                               \
    (                                                                          \
        DRGEP_FiniteRate,                                                                 \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethod                                               \
    (                                                                          \
        EFA_FiniteRate,                                                                   \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeChemistryReductionMethod                                               \
    (                                                                          \
        PFA_FiniteRate,                                                                   \
        CompChemModel,                                                         \
        Thermo                                                                 \
    );                                                                         \
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

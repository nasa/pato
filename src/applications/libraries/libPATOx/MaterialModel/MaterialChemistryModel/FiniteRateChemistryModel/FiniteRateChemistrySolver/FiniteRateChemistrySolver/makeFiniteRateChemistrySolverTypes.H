/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#ifndef makeFiniteRateChemistrySolverTypes_H
#define makeFiniteRateChemistrySolverTypes_H

#include "FiniteRateChemistrySolver.H"

#include "StandardFiniteRateChemistryModel.H"
#include "TDACFiniteRateChemistryModel.H"

#include "noFiniteRateChemistrySolver.H"
#include "EulerImplicitFiniteRate.H"
#include "odeFiniteRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFiniteRateChemistrySolverType(SS, Comp, Thermo)                              \
                                                                               \
    typedef SS<StandardFiniteRateChemistryModel<Comp, Thermo>> SS##Comp##Thermo;         \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        SS##Comp##Thermo,                                                      \
        (#SS"<" + word(StandardFiniteRateChemistryModel<Comp, Thermo>::typeName_()) + "<"\
        + word(Comp::typeName_()) + "," + Thermo::typeName() + ">>").c_str(),  \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicFiniteRateChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<SS##Comp##Thermo>                      \
        add##SS##Comp##Thermo##thermo##ConstructorTo##BasicFiniteRateChemistryModel##Comp\
##Table_; \
                                                                               \
    typedef SS<TDACFiniteRateChemistryModel<Comp, Thermo>> TDAC##SS##Comp##Thermo;       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        TDAC##SS##Comp##Thermo,                                                \
        (#SS"<" + word(TDACFiniteRateChemistryModel<Comp, Thermo>::typeName_()) + "<"    \
        + word(Comp::typeName_()) + "," + Thermo::typeName() + ">>").c_str(),  \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicFiniteRateChemistryModel<Comp>::                                                \
        add##thermo##ConstructorToTable<TDAC##SS##Comp##Thermo>                \
        add##TDAC##SS##Comp##Thermo##thermo##ConstructorTo##BasicFiniteRateChemistryModel\
##Comp##Table_;


#define makeFiniteRateChemistrySolverTypes(Comp, Thermo)                                 \
                                                                               \
    makeFiniteRateChemistrySolverType                                                    \
    (                                                                          \
        noFiniteRateChemistrySolver,                                                     \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeFiniteRateChemistrySolverType                                                    \
    (                                                                          \
        EulerImplicitFiniteRate,                                                         \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
                                                                               \
    makeFiniteRateChemistrySolverType                                                    \
    (                                                                          \
        odeFiniteRate,                                                                   \
        Comp,                                                                  \
        Thermo                                                                 \
    );                                                                         \
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along with
    OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

\mainpage Free, Open Source Software based on <a href="http://openfoam.org">OpenFOAM</a> released by <a href="https://software.nasa.gov/software/ARC-16680-1A">NASA</a>


\section about About PATO

     The Porous material Analysis Toolbox (PATO) is a modular analysis platform for multiphase porous reactive materials.\n
     It can be run as a simple Fourier heat transfer code or include more advanced features as internal decomposition (pyrolysis, vaporization), \n
     gas-gas and gas-solid chemical interactions (combustion, cracking, coking), gas species transport (convection, diffusion), \n
     and solid morphology evolutions (internal density changes, surface ablation). \n
     PATO is implemented as a C++ top level module of the open source (GNU GPL) computational fluid dynamics software program OpenFOAM. \n
     PATO also uses the open source (GNU LGPL) thermodynamics, transport, and chemistry library <a href="https://www.mutationpp.org/">Mutation++</a> produced by the von Karman Institute for Fluid Dynamics.\n
     Tutorials are provided for different types of geometries (1D, 2D, 2D-axi, 3D), boundary conditions (imposed temperature and pressure, convective boundary layer,\n
     coupled to external flow) and materials (metals, composites, wood).

\section release Open Source Release

    PATO is released Open Source to the general public by <a href="https://software.nasa.gov/software/ARC-16680-1A">NASA</a>.
\section structure PATO Directory Structure

    OpenFOAM comprises of four main directories:
    - \b src: libraries, solvers and utilities
    - \b tutorials: test-cases for material response in 1D, 2D & 3D
    - \b documentation: User's Manual, Developer's Manual & Theroy Manual

\section team Team & contacts

    To contact us about PATO : admin@pato.ac \n

    PATO was originally developed by a team of researchers at the NASA Ames Research Center (California) in collaboration with colleagues from the von Karman Institute for Fluid Dynamics (Belgium).\n


    - \b Modeling \b & \b code \b development: Dr. Jean Lachaud (jean.lachaud@u-bordeaux.fr), Jeremie Meurisse (jeremie.b.meurisse@nasa.gov)
    - \b Modeling \b & \b project \b management: Dr. Nagi N. Mansour (nagi.n.mansour@nasa.gov)
    - \b Modeling: Pr. Thierry E. Magin (thierry.magin@vki.ac.be)
    - \b Mutation++ \b development: J. B. Scoggins (james.scoggins@vki.ac.be)
    - \b Contributors: Dr. Vincent Leroy (vincent.leroy@vki.ac.be), Max Döpke (VKI)

    We would like to thank several colleagues and students for their past or current contributions to verification and/or testing of the code: \n
    Tom van Eekelen (Siemens), Ali Omidy (University of Kentucky), Dr. Francesco Panerai (NASA Ames), Dr. Ioana Cozmuta (NASA Ames), Pr. Alexandre Martin (University of Kentucky),\n
    Kyle Hanquist (University of Michigan), Julien de Muelenaere (VKI, Stanford), Dr. Victor Pozzobon (Mines d’Albi), Fransciso Torres (VKI), Florent Anstett (C la Vie), Dr. Joshua Monk (NASA Ames),\n
    Dr. John Lawless (Redwood Scientific Inc.), Jérémy Mora-Monteros (EPFL), Nikhil Banerji (EPFL), Dr. Eric Stern (NASA Ames), Pr. Pietro Congedo (INRIA Bordeaux), Michaël Rivier (INRIA Bordeaux),\n
    Xiaowen Qin (ENSMA, P’, Poitiers), Pr. Franck Richard (ENSMA, P’, Poitiers), Alessandro Turchi (VKI), Joseph Ferguson (NASA Ames, Stanford).\n

    Finally, we would like to thank the people that patiently and professionally developed OpenFOAM over the years to make it the amazing CFD platform that it became.\n
    In particular, Pr. Hrvoje Jasak and Dr. Henrik Rusche are gratefully acknowledged for their interest in PATO from the beginning and for their guidance during the second\n
    OpenFoam summer school (Numerical Modeling of Coupled Problems in Applied Physics with OpenFoam, Faculty of Mechanical Engineering and Naval Architecture, Zagreb, Croatia; 02–15 September 2009).

\*---------------------------------------------------------------------------*/

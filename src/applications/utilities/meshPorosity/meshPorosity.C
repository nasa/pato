/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    meshPorosity is a post-postprocessor to compute and print the total, voids, and fibres volumes, and porosity.
    It works considering the total domain composed of two buffer domains and a porous domain, like sketched: 
                    ------------------------------------------
                   |            |                  |          |
           inlet-->|   buffer   | porous material  |  buffer  |<--outlet
                   |   domain   |     domain       |  domain  |
                   |            |                  |          |
                    ------------------------------------------
    In order to identify the porous domain, the script looks for the patch named "fibres".
    It is then necessary to identify the walls of the porous domain with this patch name.

    The script works also in the case in which the porous domain corresponds to the whole domain 
    (no buffer domains). In this case there is no need to have a patch called "fibres.

\*---------------------------------------------------------------------------*/

// Linking librairies (class & function containers) needed by the program
#include "fvCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

    // executing case initializations
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    scalar volFibres = 0.0;         // volume of the fibres
    scalar volVoids  = 0.0;         // volume of the voids
    scalar volTotal  = 0.0;         // total volume of the domain
    scalar porosity  = 0.0;         // porosity

    label patchID;
    
    scalar xMaxDomain = mesh.bounds().max().component(0);        // to characterize the dimensions of the whole domain 
    scalar xMinDomain = mesh.bounds().min().component(0);
    scalar yMaxDomain = mesh.bounds().max().component(1);
    scalar yMinDomain = mesh.bounds().min().component(1);
    scalar zMaxDomain = mesh.bounds().max().component(2);
    scalar zMinDomain = mesh.bounds().min().component(2);

    forAll(mesh.V(), index)                 
    {
        volVoids += mesh.V()[index];   
    }

    volTotal  = (xMaxDomain-xMinDomain) * (yMaxDomain-yMinDomain) * (zMaxDomain-zMinDomain);

    if (mesh.boundaryMesh().findPatchID("fibres")>0)         // case in which the porous domain does not coincide with the whole domain
    {
        patchID = mesh.boundaryMesh().findPatchID("fibres");
        scalar xMaxFibres = max(mesh.Cf().boundaryField()[patchID].component(0));     // to characterize the dimensions of the fibres domain 
        scalar xMinFibres = min(mesh.Cf().boundaryField()[patchID].component(0));
        scalar yMaxFibres = max(mesh.Cf().boundaryField()[patchID].component(1));
        scalar yMinFibres = min(mesh.Cf().boundaryField()[patchID].component(1));
        scalar zMaxFibres = max(mesh.Cf().boundaryField()[patchID].component(2));
        scalar zMinFibres = min(mesh.Cf().boundaryField()[patchID].component(2));

        // to get the voids volume of only the fibres domain, it is needed to remove the buffer domain

        scalar correction1 = (xMaxDomain-xMinDomain)*(yMaxDomain-yMinDomain)*(zMinFibres-zMinDomain); 
        scalar correction2 = (xMaxDomain-xMinDomain)*(yMaxDomain-yMinDomain)*(zMaxDomain-zMaxFibres); 
        scalar correction3 = (xMinFibres-xMinDomain)*(yMaxDomain-yMinDomain)*(zMaxFibres-zMinFibres); 
        scalar correction4 = (xMaxDomain-xMaxFibres)*(yMaxDomain-yMinDomain)*(zMaxFibres-zMinFibres); 
        scalar correction5 = (xMaxFibres-xMinFibres)*(yMinFibres-yMinDomain)*(zMaxFibres-zMinFibres);
        scalar correction6 = (xMaxFibres-xMinFibres)*(yMaxDomain-yMaxFibres)*(zMaxFibres-zMinFibres);

        volTotal  = (xMaxFibres-xMinFibres) * (yMaxFibres-yMinFibres) * (zMaxFibres-zMinFibres);            
        volVoids  = volVoids - (correction1+correction2+correction3+correction4+correction5+correction6);   
    }

    volFibres = volTotal - volVoids;  
    porosity  = volVoids / volTotal; 

    Info<< "Total  Volume = "  << volTotal << "  [m^3] \n"
        << "Voids  Volume = " << volVoids <<  "  [m^3]  \n"
        << "Fibres Volume = " << volFibres << "  [m^3] \n"
        << "Porosity = " << porosity << "\n" << endl;


    return 0;
}


// ************************************************************************* //

// ************************************************************************* //

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
    processSets is a post-postprocessor

Usage
    Run processSets after the OpenFoam postprocessor "sample" to grab 'cloud'
    data in the 'sets' directories and re-write them in a single file.
    Currently writen in a restrictive manner to extract data from 'list_Ta_rho_s'.
    Change read and output file names to read and print other fields, and recompile.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOsampledSets.H"
#include "IOsampledSurfaces.H"
#include "OFstream.H"
#include "IFstream.H"
#include "IOstream.H"
#include "ISstream.H"
#include "RectangularMatrix.H"
#include "coordSet.H"
#include <sys/stat.h>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    // creating output directory
    if (!isDir("output")){mkdir("output", 0777);}

    std::stringstream sstrmSets;
    OFstream list_T ("output/list_T");  // opens an output file called list_T
    OFstream list_rho ("output/list_rho");  // opens an output file called list_rho
    OFstream list_p ("output/list_p");  // opens an output file called list_p

    // count number of rows (ie. probed values) at initial time
    word time0 = timeDirs[0].name();
    word fileName = "postProcessing/sampleDict/porousMat/" + time0 + "/plot_Ta_p_rho_s.xy";
    IFstream readTemp (fileName);  // opens a temporary input file                                                                                                                                                                                                                                                         
    IFstream readTemp2 (fileName);  // opens a temporary input file                                    
    Info << "The 'Sets' file must have 6 columns: x, y, z, Ta, p, rho_s." << nl;
    int columnTable = 6; int n_probes = 0;  scalar rawTableFrac, rawTableInt;   int i_raw = 0;
    scalar temp;
    while (true)
    {
        readTemp >> temp;
        if (readTemp.eof() == 1)
            break;
        i_raw++;
    }
    rawTableFrac = modf(static_cast<scalar>(i_raw) / static_cast<scalar>(columnTable), &rawTableInt);
    if (rawTableFrac != 0)
    {
        Info << "The 'Sets' file does not have a number of value that is a multiple of 6." << nl;
        Info << rawTableInt << " sets of 6 values and " << rawTableFrac * columnTable << " values have been read." << nl;
        exit (1);
    }
    else
    {
        n_probes = i_raw / columnTable;
        Info << n_probes << " probes have been found at initial time" << nl;
    }

    RectangularMatrix<scalar> Table0(n_probes, columnTable);
    RectangularMatrix<vector> probeList(n_probes, 1);
    for (int x = 0; x < n_probes; x++)
    {
        for (int i = 0; i < columnTable; i++)
        {
            readTemp2 >> Table0[x][i];
            if (i < 3)
                probeList[x][0][i] = Table0[x][i];
        }
    }
    // write header of file 'list_T'
    list_T << "// time (s)";
    for (int i = 0; i < n_probes; i++)
        list_T << " TC" << i + 1 << " (K)";
    list_T << nl;

    list_p << "// time (s)";
    for (int i = 0; i < n_probes; i++)
        list_p << " p" << i + 1 << " (Pa)";
    list_p << nl;

    list_rho << "// Solid density (ie. excluding gas) in kg/m³" << nl;
    list_rho << "// time (s)";
    for (int i = 0; i < n_probes; i++)
        list_rho << "   rho" << i + 1 ;
    list_rho << nl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        list_T << runTime.timeName() << "       ";
        list_p << runTime.timeName() << "       ";
        list_rho << runTime.timeName() << "     ";

        sstrmSets << "postProcessing/sampleDict/porousMat/" << runTime.timeName() <<  "/plot_Ta_p_rho_s.xy";
        IFstream readFile (sstrmSets.str().c_str());  // opens an input file
        IFstream readFile2 (sstrmSets.str().c_str());  // opens an input file

        if (readFile.good() == false) // checks that the input file is opened
        {
            Info << "readFile not found" << nl;
            exit (1); // exits otherwise
        }
        int n_probes_i = 0;
        scalar tempo;
        while (true)
        {
            readFile >> tempo;
            if (readFile.eof() == 1)
                break;
            n_probes_i++;
        }
        n_probes_i = n_probes_i / columnTable;

        RectangularMatrix<scalar> Table_i(n_probes_i, columnTable);
        for (int x = 0; x < n_probes_i; x++)
        {
            for (int i = 0; i < columnTable; i++)
            {
                readFile2 >> Table_i[x][i];
            }
        }

        RectangularMatrix<vector> probeList_i(n_probes_i, 1);
        for (int x = 0; x < n_probes_i; x++)
        {
            for (int i = 0; i < 3; i++)
            {
                probeList_i[x][0][i] = Table_i[x][i];
            }
        }

        int i = 0; int j = 0;
        while (i < n_probes_i && j < n_probes)
        {
            if (probeList_i[i][0] == probeList[j][0])
            {   list_T << Table_i[i][3] << "    ";
                list_p << Table_i[i][4] << "    ";
                list_rho << Table_i[i][5] << "  ";
                i++; j++;
            }
            else
            {
                list_T << 300 << "  ";
                list_rho << 0 << "  ";
                list_p << 0 << "    ";
                j++;
            } // note: no i++ to wait for next try
        }

        list_T << nl;
        list_p << nl;
        list_rho << nl;
        sstrmSets.str(""); // clear stream
    }

    Info << "File list_T, list_p, and list_rho created and saved in directory 'output'.\n" << endl;

    return 0;
}

// ************************************************************************* //

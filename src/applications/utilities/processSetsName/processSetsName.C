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
#include "fileName.H"
#include <sys/stat.h>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("inputFile");
    argList::validArgs.append("outputFile");
    argList::addOption // number of probed fields
    (
        "npf",
        "npfvalue"
    );
    argList::addOption // probe is allowed to move (surface)
    (
        "movingProbe",
        "movingProbeValue"
    );
    argList::addOption // number of probed fields
    (
        "region",
        "regionName"
    );

    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    Info << "Input File = " << argv[1] << nl;
    Info << "Base name of output files = " << argv[2] << nl;

    int npfvalue = 0;
    args.optionReadIfPresent("npf", npfvalue);
    int movingProbeValue = 0;
    args.optionReadIfPresent("movingProbe", movingProbeValue);
    word regionName = "porousMat";
    args.optionReadIfPresent("region", regionName);

    if (npfvalue < 1)
    {
        Info << "Use option -npf and provide number of probed fields in "
             "sampleDict (e.g. -npf 3)"
             << nl;
        exit(1);
    }

    std::stringstream sstrmSets, sstrmDict, sstrmFile;

    // creating a list of ofstream and naming them - with gcc 5.4.0
    /*std::vector<std::ofstream> list;
    for (int i = 0; i < npfvalue; i++)
    {
        std::stringstream fileName;
        fileName << "output/" << argv[2] << "_field_" << i;
        list.emplace_back(std::ofstream{ fileName.str().c_str() });
        fileName.str("");
    }*/

    // creating output directory
    if (!isDir("output")){mkdir("output", 0777);}

    // creating a list of ofstream and naming them - with gcc 5.0.0
    PtrList<OFstream> list(npfvalue);
    for (int i = 0; i < npfvalue; i++)
    {
        std::stringstream fileName;
        fileName << "output/" << argv[2] << "_field_" << i;
        list.set(i, new OFstream(fileName.str()));
        fileName.str("");
    }


    // verifying number of raws (ie. probed values) at initial time
    std::stringstream inputFile;
    word time0 = timeDirs[0].name();
    word folderName = "postProcessing/sampleDict/" + regionName + "/"  + time0 +"/";
    inputFile << folderName << argv[1] ;
    IFstream readTemp (inputFile.str().c_str());  // opens a temporary input file (to count)
    IFstream readTemp2 (inputFile.str().c_str());  // opens a second temporary input file (to use)
    inputFile.str("");
    Info << "The 'Sets' file must have " << npfvalue + 3 << " columns : x, y, z + nfp probed fields." << nl;
    int columnTable = npfvalue + 3;
    int n_probes = 0;
    scalar rawTableFrac, rawTableInt;
    int i_raw = 0;
    scalar temp;

    while (true)
    {
        readTemp >> temp;
        if (readTemp.eof() == 1) {break;}
        i_raw++;
    }

    rawTableFrac = modf(static_cast<scalar>(i_raw) / static_cast<scalar>(columnTable), &rawTableInt);
    if (rawTableFrac != 0)
    {
        Info << "The 'Sets' file does not have a number of value that is a multiple of " << columnTable  << nl;
        Info << "Please check your file and number of probed fields (npf)." << nl;
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
            {
                probeList[x][0][i] = Table0[x][i];
            }
        }
    }
    // write header of output files
    for ( int j = 0 ; j < npfvalue ; j++)
    {
        list[j] << "// time (s)";
        for (int i = 0; i < n_probes; i++)
        {
            list[j] << "	probe" << i + 1;
        }

        list[j] << nl;
    }

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        for ( int j = 0 ; j < npfvalue ; j++)
        {
            list[j] << runTime.timeName() << "		";
        }

        sstrmSets << "postProcessing/sampleDict/"
                  << regionName << "/"
                  << runTime.timeName() << "/" << argv[1];
        IFstream readFile(sstrmSets.str().c_str());  // opens an input file (to count)
        IFstream readFile2(sstrmSets.str().c_str()); // opens an input file (to use)

        if (readFile.good() == false) // checks that the input file is opened
        {
            Info << "Input file not found in time directory: " << runTime.timeName() << nl;
            exit(1); // exits otherwise
        }

        int n_probes_i = 0;
        scalar tempo;

        while (true)
        {
            readFile >> tempo;
            if (readFile.eof() == 1)
            {
                break;
            }

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
            if (probeList_i[i][0] == probeList[j][0] || (movingProbeValue == 1))
            {
                for (int k = 0 ; k < npfvalue ; k++)
                {
                    list[k] << Table_i[i][k + 3] << "	";
                }

                i++;
                j++;
            }
            else
            {
                for (int k = 0 ; k < npfvalue ; k++)
                {
                    list[k] << 0 << "	";
                }

                j++;
            } // note: no i++ to wait for next try
        }

        for (int k = 0 ; k < npfvalue ; k++)
        {
            list[k] << nl;
        }
        sstrmSets.str(""); // clear stream
    }

    Info << "Output files created and saved in directory 'output'.\n" << endl;

    return 0;
}

// ************************************************************************* //

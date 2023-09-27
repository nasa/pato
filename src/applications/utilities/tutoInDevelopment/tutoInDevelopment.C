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
along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
tutoInDevelopment

Description
Send an error message saying that the tutorial is in development

\*---------------------------------------------------------------------------*/

// include headers
#include "fvCFD.H"

// namespace
using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  // Options
  argList::noParallel();
  argList::removeOption("case");
  argList::removeOption("fileHandler");
  argList::removeOption("libs");
  argList::removeOption("noFunctionObjects");
  argList::addOption("errorMsg","errorMsg","Error message.");
  argList::addOption("debug","debug","Debug Switch flag to run trough this application.");
  argList args(argc, argv);
  if (!args.check()) {
    FatalError.exit();
  }
  // Debug
  Switch debug = false;
  args.optionReadIfPresent("debug",debug);
  if (debug) {
    Info << "Run through this application (debug=" << debug << ")." << endl;
  } else {
    // Error message
    word errorMsg="This tutorial is still in development.";
    args.optionReadIfPresent("errorMsg", errorMsg);
    FatalError << errorMsg << exit(FatalError);
  }
  return 0;
}

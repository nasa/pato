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

Class
    Foam::TecplotReader

Description
    Read tecplot format

SourceFiles
    TecplotReader.H
    TecplotReader.C

\*---------------------------------------------------------------------------*/

#ifndef TecplotReader_H
#define TecplotReader_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvCFD.H"
#include "IOFunctions.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class TecplotReader Declaration
\*---------------------------------------------------------------------------*/

class TecplotReader
{

 protected:

  //- Data
  List<scalarList> data_;

  //- file name
  const fileName fileName_;

  //- Line stream
  autoPtr<IStringStream> lineStream_;

  //- title
  word title_;

  //- variables
  wordList variables_;

  //- Number of zones
  int nZones_;

  //- Number of nodes in the I direction per zone
  scalarList I_zones_;

  //- Number of nodes in the J direction per zone
  scalarList J_zones_;

  //- Number of nodes in the K direction per zone
  scalarList K_zones_;

  //- Zone index
  int zoneI;

  //- Zone type
  wordList zoneType;

  //- I,J,K vector per region and per node
  List<List<vector> > ijk_nodes_;

  //- X,Y,Z vector per region and per node
  List<List<vector> > xyz_nodes_;



 public:

  //- Runtime type information
//  TypeName("TecplotReader");

  // Constructors

  //- Construct from mesh and phase name
  TecplotReader(const fileName& file);

  //- Destructor
  virtual ~TecplotReader();

  //- Member Functions

  //- Read tecplot file
  void readFile();

  //- Read the header of the tecplot file
  void readHeader(IFstream& dataFile_);

  //- Read the zone headers of the tecplot file
  void readZoneHeader(IFstream& dataFile_);

  //- Read the zone data of the tecplot file
  void readZoneData(IFstream& dataFile_);

  //- Write the plot 3D ASCII file
  void writePlot3D(word name, List<scalarList>& data);

  //- Split the data into xyz_nodes and ijk_nodes
  void splitData(const List<scalarList>& data);

  //- Return the number of zones
  int nZones();

  //- Return the number of nodes in the I direction per zone
  scalarList I_zones();

  //- Return the number of nodes in the J direction per zone
  scalarList J_zones();

  //- Return the number of nodes in the K direction per zone
  scalarList K_zones();

  //- Return the data
  const List<scalarList>& data();

  //- Return ijk_nodes
  const List<List<vector>>& ijk_nodes();

  //- Return xyz_nodes
  const List<List<vector>>& xyz_nodes();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

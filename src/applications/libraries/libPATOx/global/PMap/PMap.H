/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::PMap

Description
    Pointers Map

SourceFiles
    PMap.C

\*---------------------------------------------------------------------------*/

#ifndef PMap_H
#define PMap_H

#include "fvCFD.H"
#include <map>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class T> class PMap;

/*---------------------------------------------------------------------------*\
                     class PMap Declaration
\*---------------------------------------------------------------------------*/

template <class T>
class PMap
{

 private:

  // Private Member Functions

  //- Map of pointers
  std::map<word, T*> ptrs_;

 public:

  // Types
  typedef typename std::map<word, T*>::iterator iterator;

  // Constructors

  //- Null Constructor
  PMap() {}

  //- Destructor
  ~PMap() {}

  //- Insert ptr element
  inline void insert(const word w, T* t);

  //- Insert ref element
  inline void insert(const word w, T& t);

  //- Insert st element
  inline void insert(const word w, ITstream& st);

  //- Get size of the pointers map
  inline int size() const;

  //- Clear the pointers map
  inline void clear();

  //- Find a key in the pointers map
  inline iterator find(const word w);

  //- end of the pointers map
  inline iterator end();

  //- Found a key in the pointers map
  inline bool found(const word w);

  //- Found a key in the pointers map (const function)
  inline bool found(const word w) const;

  //- Return element reference
  inline const T& operator[](const word w) const;

  //- Return element reference
  inline T& operator[](const word w);

  //- Initialize the pointers list
  inline void operator=(const PMap<T>& pmap);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PMapI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

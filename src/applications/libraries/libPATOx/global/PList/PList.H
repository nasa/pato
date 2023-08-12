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
    Foam::PList

Description
    Pointers List

SourceFiles
    PList.C

\*---------------------------------------------------------------------------*/

#ifndef PList_H
#define PList_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class T> class PList;

/*---------------------------------------------------------------------------*\
                     class PList Declaration
\*---------------------------------------------------------------------------*/

template <class T>
class PList
{

 private:

  // Private Member Functions

  //- List of pointers
  List<T*> ptrs_;

 public:

  // Constructors

  //- Null Constructor
  PList() {}

  //- Construct from int
  PList(const int n) {
    ptrs_.resize(n);
  }

  //  List<MapField*>
  //- Construct from List of pointers
  PList(const List<T*>& p) {
    ptrs_.resize(p.size());
    forAll(p,i) {
      this->set(i, p[i]);
    }
  }
//  //- Construct from other PList
//  PList(const PList<T>& p) {
//    ptrs_.resize(p.size());
//    forAll(p,i){
//        this->ptrs_[i] = (p[i]).clone().ptr();
//    }
//  }

//  //- Construct from other List
//  PList(const List<T>& p) {
//    ptrs_.resize(p.size());
//    forAll(p,i){
//        this->ptrs_[i] = new T(p[i]);
//    }
//  }

  //- Destructor
  ~PList() {}

  //- Set ptr element
  inline void set(const label l, T* t);

  //- Set ref element
  inline void set(const label l, T& t);

  //- Set st element
  inline void set(const label l, ITstream& st);

  //- Append ptr element
  inline void append(T* t);

  //- Append ref element
  inline void append(T& t);

  //- Append st element
  inline void append(ITstream& st);

  //- Resize the pointers list
  inline void resize(const int n);

  //- Get size of the pointers list
  inline int size() const;

  //- Clear the pointers list
  inline void clear();

  //- Return const element reference
  inline const T& operator[](const label l) const;

  //- Return element reference
  inline T& operator[](const label l);

  //- Initialize the pointers list
  inline void operator=(const std::initializer_list<T*>& list);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

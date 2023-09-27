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

\*---------------------------------------------------------------------------*/

#include "PMap.H"

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class T>
inline void Foam::PMap<T>::insert(const word w, T* ptr)
{
  if (found(w)) {
    FatalError << "PMap: key " << (string) w << " already exists." << exit(FatalError);
  }
  ptrs_[w]=ptr;
}

template<class T>
inline void Foam::PMap<T>::insert(const word w, T& ref)
{
  if (found(w)) {
    FatalError << "PMap: key " << (string) w << " already exists." << exit(FatalError);
  }
  ptrs_[w]=&ref;
}

template<class T>
inline void Foam::PMap<T>::insert(const word w, Foam::ITstream& st)
{
  if (found(w)) {
    FatalError << "PMap: key " << (string) w << " already exists." << exit(FatalError);
  }
  ptrs_[w]=new T(st);
}

template<class T>
inline int Foam::PMap<T>::size() const
{
  return ptrs_.size();
}

template<class T>
inline void Foam::PMap<T>::clear()
{
  ptrs_.clear();
}

template<class T>
inline typename Foam::PMap<T>::iterator Foam::PMap<T>::find(const word w)
{
  return ptrs_.find(w);
}

template<class T>
inline typename Foam::PMap<T>::iterator Foam::PMap<T>::end()
{
  return ptrs_.end();
}

template<class T>
inline bool Foam::PMap<T>::found(const word w)
{
  return find(w) != end();
}

template<class T>
inline bool Foam::PMap<T>::found(const word w) const
{
  return find(w) != end();
}

template<class T>
inline const T& Foam::PMap<T>::operator[](const word w) const
{
  if (!found(w)) {
    FatalError << "PMap: key " << (string) w << " does not exist." << exit(FatalError);
  }
  return *ptrs_[w];
}

template<class T>
inline T& Foam::PMap<T>::operator[](const word w)
{
  if (!found(w)) {
    FatalError << "PMap: key " << (string) w << " does not exist." << exit(FatalError);
  }
  return *ptrs_[w];
}

template<class T>
inline void Foam::PMap<T>::operator=(const PMap<T>& pmap)
{
  ptrs_=pmap.ptrs_;
}


// ************************************************************************* //


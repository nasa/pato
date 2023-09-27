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

#include "PList.H"

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class T>
inline void Foam::PList<T>::set(const label l, T* ptr)
{
  ptrs_[l]=ptr;
}

template<class T>
inline void Foam::PList<T>::set(const label l, T& ref)
{
  ptrs_[l]=&ref;
}

template<class T>
inline void Foam::PList<T>::set(const label l, Foam::ITstream& st)
{
  ptrs_[l]=new T(st);
}

template<class T>
inline void Foam::PList<T>::append(T* ptr)
{
  ptrs_.append(ptr);
}

template<class T>
inline void Foam::PList<T>::append(T& ref)
{
  ptrs_.append(&ref);
}

template<class T>
inline void Foam::PList<T>::append(Foam::ITstream& st)
{
  ptrs_.append(new T(st));
}

template<class T>
inline void Foam::PList<T>::resize(const int n)
{
  ptrs_.resize(n);
}

template<class T>
inline int Foam::PList<T>::size() const
{
  return ptrs_.size();
}

template<class T>
inline void Foam::PList<T>::clear()
{
  ptrs_.clear();
}

template<class T>
inline const T& Foam::PList<T>::operator[](const label l) const
{
  return *ptrs_[l];
}


template<class T>
inline T& Foam::PList<T>::operator[](const label l)
{
  return *ptrs_[l];
}

template<class T>
inline void Foam::PList<T>::operator=(const std::initializer_list<T*>& init_list)
{
  ptrs_=init_list;
}


// ************************************************************************* //


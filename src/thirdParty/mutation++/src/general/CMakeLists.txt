#
# Copyright 2014-2017 von Karman Institute for Fluid Dynamics (VKI)
#
# This file is part of MUlticomponent Thermodynamic And Transport
# properties for IONized gases in C++ (Mutation++) software package.
#
# Mutation++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Mutation++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with Mutation++.  If not, see
# <http://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 2.6)

add_sources(mutation++
    Mixture.cpp
    MixtureOptions.cpp
)

# Build the checkmix code
add_executable(checkmix checkmix.cpp)
target_link_libraries(checkmix mutation++)
install(TARGETS checkmix DESTINATION bin)

# Build the mppequil code
add_executable(mppequil mppequil.cpp)
target_link_libraries(mppequil mutation++)
install(TARGETS mppequil DESTINATION bin)

add_executable(bprime bprime.cpp)
target_link_libraries(bprime mutation++)
install(TARGETS bprime DESTINATION bin)

# Install the header files
install(FILES GlobalOptions.h DESTINATION include/mutation++)
install(FILES mutation++.h DESTINATION include/mutation++)
install(FILES Mixture.h DESTINATION include/mutation++)
install(FILES MixtureOptions.h DESTINATION include/mutation++)
install(FILES Constants.h DESTINATION include/mutation++)
install(FILES Errors.h DESTINATION include/mutation++)




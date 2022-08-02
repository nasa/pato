#!/usr/bin/python3
# -*- coding: utf-8 -*-

################################################################################
'''
    Based on Copyright (C) 2018 Yuusha
    Modified by PATO-community
'''
################################################################################
'''
License
    This file is part of Yuusha contribution to OpenFOAM.

    This contribution is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This contribution is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this contribution.  If not, see <http://www.gnu.org/licenses/>.
'''
################################################################################

from re import sub

class thermophysicalTable() :
    """Base class for manipulate thermophysical table
    needed by OpenFOAM tabulated thermophysical properties model."""

    def __init__(self) :
        self.table = []
        self.fileName = ''

    def read(self, fileName='') :
        """ Read a tabulated thermophysical file and store it into table"""
        
        with open(fileName, 'r') as f :
            for line in f :
                tmpList = sub('\(+|\)+', ' ', line).strip().split()
                try :
                    count = len(tmpList) - 1
                    tupleList=[]
                    for i in range(1,2):
                     tupleList.append((tmpList[1], tmpList[j]))
                    self.table.append([tmpList[0], tupleList])
                except IndexError :
                    continue
                    
  
    def readInver(self, fileName='') :
        """ Read a tabulated thermophysical file and store it into table"""

        with open(fileName, 'r') as f :
            for line in f :
                tmpList = sub('\(+|\)+', ' ', line).strip().split()
                try :
                    count = len(tmpList) - 1
                    tupleList=[]
                    for i in range(1):
                     tupleList.append((tmpList[0], tmpList[j]))
                    self.table.append([tmpList[1], tupleList])
                except IndexError :
                    continue

    def write(self, fileName='') :
        """ Write a thermopysical file """
        del self.table[0]
        with open(fileName, 'w') as f :
            f.write('(\n')
            for element in self.table :
                stringOut=""
                for elem in element[1][:] :
                    stringOut += '({} {}) '.format(elem[0], elem[1])
                f.write('( {} ({}))\n'.format(element[0], stringOut))
            f.write(')')

    def importation(self, fileName='', ext = '') :
        """ Import thermophysical data from a file.
        Only CSV is implemented yet """

    def transpose(self) :
        """ Invert lines and columns in a tabulated thermophysical list """

        tableTmp = list(self.table)
        self.table = []
        for element in tableTmp :
            for elem in element[1] :
                find = False
                for pres in self.table :
                    if elem[0] == pres[0] :
                        pres[1].append((element[0], elem[1]))
                        find = True
                        break
                if not find :
                    self.table.append([elem[0], [(element[0], elem[1])]])



thermo = thermophysicalTable()
j=2
thermo.read('flowEquilProp_NASA9.txt')
thermo.transpose()
thermo.write('hTable')

thermo = thermophysicalTable()
j=3
thermo.read('flowEquilProp_NASA9.txt')
thermo.transpose()
thermo.write('cpTable')

thermo = thermophysicalTable()
j=4
thermo.readInver('flowEquilProp_NASA9.txt')
thermo.transpose()
thermo.write('muTable') 

thermo = thermophysicalTable()
j=5
thermo.readInver('flowEquilProp_NASA9.txt')
thermo.transpose()
thermo.write('kappaTable')

thermo = thermophysicalTable()
j=6
thermo.readInver('flowEquilProp_NASA9.txt')
thermo.transpose()
thermo.write('densityTable')

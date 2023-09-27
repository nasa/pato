#!/usr/bin/env python3

# ---------------------------------------------------------------------------

# Notices:

#     Copyright © 2010 United States Government as represented by the
#     Administrator of the National Aeronautics and Space Administration.  All
#     Rights Reserved.

# Disclaimers

#     No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY
#     OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
#     LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
#     SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE
#     SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION,
#     IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES
#     NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY
#     PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE
#     PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT
#     SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND
#     LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL
#     SOFTWARE, AND DISTRIBUTES IT "AS IS."

#     Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
#     AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS,
#     AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT
#     SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES
#     ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR
#     RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL
#     INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS
#     AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT
#     PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE
#     THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.

#   ---------------------------------------------------------------------------

import numpy as np
import os

path = os.getcwd()

tol = 1e-6           # define criteria of the minimum error
error = 10
h = 1000             # define the initial hv
dh = 1              # define the step size to approximate  derviative
k = 0.1
dk = 0.04
i = 0               # number of iterations
Te = np.array([[293., 308.15716783, 322.13049385, 333.50629854, 342.91414445,
                350.78196863]])  # experiment data of temperature
s1 = 'Hv0             Hv0             [1 -1 -3 -1 0 0 0]              '
s3 = ';\n'

s4 = 'kiCoef          kiCoef          [0 0 0 0 0 0 0]         '
s6 = 'kjCoef          kjCoef          [0 0 0 0 0 0 0]         '
s7 = 'kkCoef          kkCoef          [0 0 0 0 0 0 0]         '

n = 6    # define the length of temperature data

while (error > tol):
    os.chdir(path)
    os.system(path+"/Allclean  >> /dev/null 2>&1")
    # ..........................#
    f = open(path+"/constant/CarbonFiberPreform/constantProperties", 'r+')
    flist = f.readlines()
    s2 = str(h)
    s5 = str(k)
    # s2 = s2[1:-1]
    flist[38] = s1+s2+s3
    flist[58] = s4+s5+s3
    flist[59] = s6+s5+s3
    flist[60] = s7+s5+s3
    f = open(path+"/constant/CarbonFiberPreform/constantProperties", 'w+')
    f.writelines(flist)
    f.close()

    # change  Hv0 in constantProperties
    # ...........................................#

    # Allrun_parallel, solve the energy equations
    os.system(path+"/Allrun_basic  >> /dev/null 2>&1")
    # os.system('/home/sliu/OpenFOAM/python3/python3  script.py')

    T1 = np.zeros(n)
    index = 0
    for time in range(0, n, 1):
        os.chdir(path+"/postProcessing/singleGraph/porousMat/"+str(time)+"/")
        T1[index] = np.loadtxt(open('line_p_Tg_Ta.xy', 'r'))[100, 3]
        index = index + 1
    print("Temperature at iteration n: {}".format(T1))

    # os.system('PATOx/Allclean')
    # ............................#
    # Grab temperature data at each time

    hh = h + dh               # approximate derviative by (s(h + dh)-s(h))/(dh)
    f = open(path+"/constant/CarbonFiberPreform/constantProperties", 'r+')
    flist = f.readlines()
    s2 = str(hh)
    flist[38] = s1+s2+s3
    f = open(path+"/constant/CarbonFiberPreform/constantProperties", 'w+')
    f.writelines(flist)
    f.close()
    os.system(path+"/Allrun_basic  >> /dev/null 2>&1")

    T2 = np.zeros(n)
    index = 0
    for time in range(0, n, 1):
        os.chdir(path+"/postProcessing/singleGraph/porousMat/"+str(time)+"/")
        T2[index] = np.loadtxt(open('line_p_Tg_Ta.xy', 'r'))[100, 3]
        index = index + 1
    print("Temperature at time n+1: {}".format(T2))

    # the difference between experiment data (Te)and predicted data(T1)
    dT = Te - T1
    Z = (T2 - T1)/dh          # approximate derviative by (s(h + dh)-s(h))/(dh)
    A1 = np.dot(dT, Z.T)
    A2 = np.dot(Z, Z.T)
    A3 = A1/A2
    h1 = h+A3   # solve the equation (12) in the article.  estimate of next hv

    errorh = abs((h1-h)/h1)     # solve the error between two hv
    h1 = str(h1)
    h1 = h1[1:-1]
    h1 = float(h1)
    # h = h1

    # ...........................................#
    # change  Hv0 in constantProperties

    kk = k + dk               # approximate derviative by (s(h + dh)-s(h))/(dh)
    f = open(path+"/constant/CarbonFiberPreform/constantProperties", 'r+')
    flist = f.readlines()
    s2 = str(h)
    flist[38] = s1+s2+s3
    s5 = str(kk)
    flist[58] = s4+s5+s3
    flist[59] = s6+s5+s3
    flist[60] = s7+s5+s3
    f = open(path+"/constant/CarbonFiberPreform/constantProperties", 'w+')
    f.writelines(flist)
    f.close()
    os.system(path+"/Allrun_basic  >> /dev/null 2>&1")

    T3 = np.zeros(n)
    index = 0
    for time in range(0, n, 1):
        os.chdir(path+"/postProcessing/singleGraph/porousMat/"+str(time)+"/")
        T3[index] = np.loadtxt(open('line_p_Tg_Ta.xy', 'r'))[100, 3]
        index = index + 1
    print("Temperature at time n+2: {}".format(T3))

    # the difference between experiment data (Te) and predicted data(T1)
    dT = Te - T1
    Y = (T3 - T1)/dk          # approximate derviative by (s(h + dh)-s(h))/(dh)
    A1 = np.dot(dT, Y.T)
    A2 = np.dot(Y, Y.T)
    A3 = A1/A2
    k1 = k+1*A3  # solve the equation (12) in the article.  estimate of next hv

    errork = abs((k1-k)/k1)     # solve the error between two hv
    k1 = str(k1)
    k1 = k1[1:-1]
    k1 = float(k1)

    k1 = max(0.000001, k1)
    k1 = min(k1, 4)
    h = h1
    k = k1
# ...........................................#
# change  Hv0 in constantProperties

    error = max(errorh, errork)
    S = np.dot(dT, dT.T)
    i = i+1

    print("{:<12} {:<18} {:<21} {:<18} {:<21} {:<10}"
          .format("iteration", "h_v", "errorHv", "k", "errorK", "S(i-1)"))
    print("{:<12} {:<18} {:<21} {:<18} {:<21} {:<10}\n"
          .format(i, h, errorh[0], k, errork[0], S[0][0]))

    if i > 200:
        break

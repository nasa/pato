/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    Math functions (derivative, integration, interpolation)

Details
    - \b derivative: second order, finite difference
    - \b integration: first order, trapezoidal
    - \b interpolation: first order, linear

\*---------------------------------------------------------------------------*/

#ifndef mathFunctionsFoam_H
#define mathFunctionsFoam_H

#include "fvCFD.H"
#include "IOFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Compute the distance between two vectors
inline scalar dist(vector v1, vector v2)
{
  return sqrt(pow(v1.x()-v2.x(),2) + pow(v1.y()-v2.y(),2) + pow(v1.z()-v2.z(),2));
}

//- Comparison between two scalar values
inline bool assertEquals(scalar val1, scalar val2)
{
  bool result_ = false;
  if(val1==0||val2==0) {
    if(fabs(val2-val1)<1e-5) {
      result_= true;
    } else {
      result_ = false;
    }
  } else {
    if(fabs(1-val2/val1)<1e-5) {
      result_ = true;
    } else {
      result_ = false;
    }
  }
  return result_;
}

//- derivative of a discrete function (divided differences)
inline scalarList secondOrderDerivative(scalarList x, scalarList f_x)
{
  if (x.size() != f_x.size()) {
    FatalErrorInFunction << "Derivative failed. x and f_x have a different size." << exit(FatalError);
  }
  labelList visitOrder;
  sortedOrder(x, visitOrder);
  x = scalarList(x, visitOrder);
  f_x = scalarList(f_x, visitOrder);
  scalarList df_dx(f_x.size());

  forAll(f_x, i) {
    //1st order: forward difference
    if (i==0) {

      df_dx[i] = (f_x[i+1] - f_x[i])/(x[i+1]-x[i]);
      continue;
    }

    //1st order: backward difference
    if (i == f_x.size()-1) {
      df_dx[i] = (f_x[i] - f_x[i-1])/(x[i]-x[i-1]);
      continue;
    }

    //2d order: central difference
    if (i!=0 && i!= f_x.size()-1) {
      df_dx[i]=(f_x[i+1] - f_x[i-1])/(x[i+1]-x[i-1]);
      continue;
    }

  }
  return df_dx;
}

//- Linear interpolation from a table
inline scalar linearInterpolation(scalarList x, scalarList f_x, scalar x0)
{
  if (x.size() != f_x.size()) {
    FatalErrorInFunction << "Linear interpolation failed." << exit(FatalError);
  }
  scalar f_x0;
  labelList visitOrder;
  sortedOrder(x, visitOrder);
  x = scalarList(x, visitOrder);
  f_x = scalarList(f_x, visitOrder);

  if (x0 < x[0]) {
//    WarningInFunction << "Interpolation lower bound." << endl;
    f_x0 = f_x[0];
    return  f_x0;
  }
  if (x0 > x[x.size()-1]) {
//      WarningInFunction << "Interpolation upper bound." << endl;
    f_x0 = f_x[x.size()-1];
    return  f_x0;
  }
  forAll(x, i) {
    if (x[i] ==x0) {
      f_x0 = f_x[i];
      return  f_x0;
    }
  }
  forAll(x, i) {
    if (i>0) {
      if (x0 > x[i-1] && x0 < x[i]) {
        scalar f_x0 = (f_x[i] - f_x[i-1])/(x[i] - x[i-1])*(x0 - x[i-1]) + f_x[i-1];
        return f_x0;
      }
    }
  }

  FatalErrorInFunction << "Linear interpolation failed." << exit(FatalError);
  return 0;
}

//- Bilinear interpolation from a table (x,y,f(x,y))
inline scalar bilinearInterpolation(scalarList x, scalarList y, scalarList f_xy, scalar x0, scalar y0)
{
  // Verify x,y,f
  if (x.size() != f_xy.size()) {
    FatalErrorInFunction << "Bilinear interpolation failed. x.size() != f_xy.size()" << exit(FatalError);
  }
  if (y.size() != f_xy.size()) {
    FatalErrorInFunction << "Bilinear interpolation failed. y.size() != f_xy.size()" << exit(FatalError);
  }
  if (x.size() != y.size()) {
    FatalErrorInFunction << "Bilinear interpolation failed. x.size() != y.size()" << exit(FatalError);
  }
  if (x.size() == 0) {
    FatalErrorInFunction << "Bilinear interpolation failed. Table size == 0" << exit(FatalError);
  }

  // sort x
  labelList visitOrder;
  sortedOrder(x, visitOrder);
  x = scalarList(x, visitOrder);
  y = scalarList(y, visitOrder);
  f_xy = scalarList(f_xy, visitOrder);

  // construct table x (remove the same values)
  scalarList table_x;
  List<scalarList> table_y;
  List<scalarList> table_f;
  forAll(x, i) {
    if (!foundInList<scalar>(x[i],table_x)) {
      table_x.append(x[i]);
    }
  }

  // construct table y and table f
  table_y.resize(table_x.size());
  table_f.resize(table_x.size());
  forAll(table_x, i) {
    forAll(y, j) {
      if(x[j]==table_x[i]) {
        table_y[i].append(y[j]);
        table_f[i].append(f_xy[j]);
      }
    }
  }

  // sort table y and table f
  forAll(table_y, xI) {
    sortedOrder(table_y[xI], visitOrder);
    table_y[xI] = scalarList(table_y[xI], visitOrder);
    table_f[xI] = scalarList(table_f[xI], visitOrder);
  }

  // verify not same values for table y
  forAll(table_y, xI) {
    forAll(table_y[xI],yI) {
      if(yI > 0) {
        if(table_y[xI][yI] == table_y[xI][yI-1]) {
          FatalErrorInFunction << "Bilinear interpolation failed. Two y values are the same." << exit(FatalError);
        }
      }
    }
  }
  scalar f_x0y0=0; // interpolated solution
  scalar f_x1y1=0;
  scalar f_x2y1=0;
  scalar f_x1y2=0;
  scalar f_x2y2=0;
  scalar x1=0;
  scalar x2=0;
  scalar y11=0;
  scalar y12=0;
  scalar y21=0;
  scalar y22=0;


  if (x0 < min(table_x)) {
    //    WarningInFunction << "Interpolation lower bound." << endl;
    x0 = min(table_x);
  }
  if (x0 > max(table_x)) {
    //    WarningInFunction << "Interpolation upper bound." << endl;
    x0 = max(table_x);
  }

  forAll(table_x, i) {
    if (i>0) {
      if (x0 >= table_x[i-1] && x0 <= table_x[i]) {
        x1 = table_x[i-1];
        x2 = table_x[i];

        scalar y0_1 = y0;
        if (y0_1 < min(table_y[i-1])) {
          //    WarningInFunction << "Interpolation lower bound." << endl;
          y0_1 = min(table_y[i-1]);
        }
        if (y0_1 > max(table_y[i-1])) {
          //    WarningInFunction << "Interpolation upper bound." << endl;
          y0_1 = max(table_y[i-1]);
        }

        scalar y0_2 = y0;
        if (y0_2 < min(table_y[i])) {
          //    WarningInFunction << "Interpolation lower bound." << endl;
          y0_2 = min(table_y[i]);
        }
        if (y0_2 > max(table_y[i])) {
          //    WarningInFunction << "Interpolation upper bound." << endl;
          y0_2 = max(table_y[i]);
        }
        forAll(table_y[i-1], j) {
          if (j>0) {
            if (y0_1 >= table_y[i-1][j-1] && y0_1 <= table_y[i-1][j]) {
              y11 = table_y[i-1][j-1];
              y12 = table_y[i-1][j];
              f_x1y1 = table_f[i-1][j-1];
              f_x1y2 = table_f[i-1][j];
            }
          }
        }
        forAll(table_y[i], j) {
          if (j>0) {
            if (y0_2 >= table_y[i][j-1] && y0_2 <= table_y[i][j]) {
              y21 = table_y[i][j-1];
              y22 = table_y[i][j];
              f_x2y1 = table_f[i][j-1];
              f_x2y2 = table_f[i][j];
            }
          }
        }
        break;
      }
    }
  }

  scalar yR1= (y21 - y11) / (x2 - x1) * ( x0 - x1) + y11;
  scalar yR2= (y22 - y12) / (x2 - x1) * ( x0 - x1) + y12;
  if (((x2-x1)==0 ) || ((yR2 - yR1)==0)) {
    Info << "x1=" <<  x1 << endl;
    Info << "x2=" <<  x2 << endl;
    Info << "yR1=" << yR1  << endl;
    Info << "yR2=" << yR2 << endl;
    FatalErrorInFunction << "Bilinear interpolation failed. ((x2-x1)==0 ) ||  ((yR2 - yR1)==0)" << exit(FatalError);
  }
  // Bilinear interpolation
  scalar R1 = (x2- x0)/(x2-x1)* f_x1y1 + (x0 - x1)/(x2-x1)* f_x2y1;
  scalar R2 = (x2- x0)/(x2-x1)* f_x1y2 + (x0 - x1)/(x2-x1)* f_x2y2;
  f_x0y0 = (yR2 - y0)/(yR2 - yR1)*R1 + (y0 - yR1)/(yR2 - yR1)*R2;

  return f_x0y0;
}


//- Integrate f(x) along the x axis using the composite trapezoidal rule. I = \int^b_a f(x) dx
inline scalar trapezoidalIntegration(scalarList x, scalarList f_x, scalar a, scalar b)
{
  if (x.size() != f_x.size()) {
    FatalErrorInFunction << "Trapezoidal integration failed." << exit(FatalError);
  }
  labelList visitOrder;
  sortedOrder(x, visitOrder);
  x = scalarList(x, visitOrder);
  f_x = scalarList(f_x, visitOrder);

  scalar a_ = a;
  scalar b_ = b;
  if (a_ > b_) {
    a = b_;
    b = a_;
  }
  scalar integral = 0;
  forAll(f_x, i) {
    if ((x[i]>=a) && (x[i]<=b)) {
      if(i>0) {
        scalar delta_x = x[i] - x[i-1];
        integral += (f_x[i-1]+f_x[i])/2*delta_x;
      }
    }
  }

  if (a_ > b_) {
    return -integral;
  }

  return integral;
}

//- Integrate f(x) along the x axis using the composite trapezoidal rule. I = \int_dom(x) f(x) dx
inline scalar trapezoidalIntegrationFullDomain(scalarList x, scalarList f_x)
{

  scalar a = min(x);
  scalar b = max(x);

  return trapezoidalIntegration(x, f_x, a, b);
}

//- Calculation enthalpy from heat capacity
inline scalar calculationEnthalpy(scalar T, scalar h_ref, scalar T_ref, scalarList T_list, scalarList cp_list)
{
  scalar h =  h_ref + trapezoidalIntegration(T_list, cp_list, T_ref, T);
  return h;
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

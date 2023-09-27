/*---------------------------------------------------------------------------* \
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
  createBprimeTable

  Description
  Create B' table using multi elemental composition at the surface

  \*---------------------------------------------------------------------------*/

// include headers
#include "fvCFD.H"
#include "PATOx.H"
#include "IOFunctions.H"
#include "mathFunctions.H"
#undef Log // conflict between OpenFoam and Mutation++ on the alias "Log" -> undefined here from OpenFoam.
#include <mutation++/mutation++.h>
#include "IOmanip.H"

// namespace
using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Compute the equilibrium
Eigen::VectorXd equilibrium(const Mutation::Mixture& mix, const double& p, const double& T, const Eigen::VectorXd& Z_wsg, Eigen::VectorXd& X_sg )
{
  int ne = mix.nElements(); // number of elements
  Eigen::VectorXd Zx_wsg(Eigen::VectorXd::Zero(ne)); // Gaseous mass fraction (solid + gas)
  Eigen::VectorXd Zx_g(Eigen::VectorXd::Zero(ne)); // Gaseous mole fraction (gas)
  Eigen::VectorXd Z_g(Eigen::VectorXd::Zero(ne)); // Gaseous mass fraction (gas)
  const int ng = mix.nGas(); // Number of gases

  mix.convert<Mutation::Thermodynamics::YE_TO_XE>(Z_wsg.data(), Zx_wsg.data()); // Convert to mass fraction

  mix.equilibriumComposition(T, p, Zx_wsg.data(), X_sg.data(), Mutation::Thermodynamics::IN_PHASE); //IN_PHASE); //GLOBAL); // Compute the equilibrium in Mutation++
  Zx_g = mix.elementMatrix().topRows(ng).transpose()*X_sg.head(ng); // Gaseous mass fraction from the solid, liquid and gaseous mole fraction
  double m_z = Zx_g.sum(); // sum
  Zx_g /= m_z; // sum = 1
  mix.convert<Mutation::Thermodynamics::XE_TO_YE>(Zx_g.data(), Z_g.data()); // Convert to mass fraction
  return Z_g;
}

// Compute the surface mass balance
void surfaceMassBalance(const Mutation::Mixture& mix, const List<Tuple2<word,scalar>>& elemOnSurface, const Eigen::VectorXd& z_e, const Eigen::VectorXd& z_pg, const double T,
                        const double p, const double Bg, double &Bc, double &hw, Eigen::VectorXd& X_sg, const double& Bc_max)
{
  const int ns = mix.nSpecies(); // number of species
  const int ne = mix.nElements(); // number of elements
  const int ng = mix.nGas(); // number of gaseous species

  Eigen::VectorXd z_wsg(Eigen::VectorXd::Zero(ne)); // Solid, liquid and gaseous mass fraction at the wall
  Eigen::VectorXd z_wg(Eigen::VectorXd::Zero(ne)); // Gaseous mass fraction at the wall

  // Initialize the wall element fractions to be the pyrolysis gas fractions
  double sum = 0.0;
  z_wsg=z_e+Bg*z_pg;

  // Use "large" amount to simulate infinite char
  double large = max(100.0*Bg, 200.0);
  for (int eI=0; eI<elemOnSurface.size(); eI++) {
    int i = mix.elementIndex(elemOnSurface[eI].first());
    double large_elem = (elemOnSurface[eI].second() * large);
    z_wsg[i] += large_elem;
  }
  sum=z_wsg.sum();
  z_wsg/=sum;
  z_wg = equilibrium(mix,  p, T, z_wsg, X_sg); // Gaseous mass fraction at the wall

  double z_w_surf = 0;
  double z_e_surf = 0;
  double z_pg_surf = 0;

  for (int eI=0; eI<elemOnSurface.size(); eI++) {
    int i = mix.elementIndex(elemOnSurface[eI].first());
    z_w_surf += z_wg[i];
    z_e_surf += z_e[i];
    z_pg_surf += z_pg[i];
  }

  if (z_w_surf == 1) {
    Bc = Bc_max;
  } else {
    Bc = (z_e_surf + Bg*z_pg_surf - z_w_surf*(1.0 + Bg)) / (z_w_surf - 1.0); // Dimensionless char blowing rate
    Bc = max(Bc, 0.0);
  }

  double * hwi = new double[ns]; // Wall species enthalpies
  for(int i=0; i<ns; i++) {
    hwi[i]=0;
  }
  const double RU = Mutation::RU;

  // Compute the gas enthalpy
  mix.speciesHOverRT(T, hwi);
  double mwg = 0.0;
  for (int j = 0; j < ng; ++j) {
    mwg += mix.speciesMw(j) * X_sg[j];
  }
  hw = 0.0;
  for (int i = 0; i < ng; ++i) {
    hw += X_sg[i] * hwi[i];
  }
  hw *= RU * T / mwg;
}

// Compute SiC-SiO2
List<Tuple2<word,scalar>> computeSiCSiO2SurfComp(
    const Mutation::Mixture& mix, const double& zC)
{
  Eigen::ArrayXd mm (mix.nSpecies());
  forAll(mm, i) {
    mm[i] = mix.speciesMw(i);
  }

  // Arranging Species
  const int ns_e = 3;
  const int eSi = 0;
  string s_Si = "Si";
  const int eO = 1;
  string s_O = "O";
  const int eC = 2;
  string s_C = "C";

  // Computing molar masses
  const double mmSi = mm(mix.speciesIndex(s_Si));
  const double mmO  = mm(mix.speciesIndex(s_O));
  const double mmC  = mm(mix.speciesIndex(s_C));
  const double mmSiC = mmSi + mmC;
  const double mmSiO2 = mmSi + 2.*mmO;

  // Mass Fraction of Surface Species
  double zSiC = zC + mmSi/mmC*zC;
  double zSiO2 = 1.-zSiC;

  if (zSiO2 < 0.) {
    std::cerr << "Such a composition is not possible. yC = " << zC << std::endl;
    exit(EXIT_FAILURE);
  }

  List<Tuple2<word,double>> zs_e(ns_e);
  zs_e[eSi].first() = s_Si;
  zs_e[eSi].second() = mmSi/mmSiC*zSiC + mmSi/mmSiO2*zSiO2;
  zs_e[eO].first() = s_O;
  zs_e[eO].second() = 2.*mmO/mmSiO2*zSiO2;
  zs_e[eC].first() = s_C;
  zs_e[eC].second() = zSiC*mmC/mmSiC;

  return zs_e;
}

int main(int argc, char *argv[])
{
  // Read the input file
  argList::noParallel();
  argList::validArgs.append("inputFile");
  argList args(argc, argv);

  if (!args.check()) {
    FatalError.exit();
  }

  fileName file_ = args[1];

  Foam::Time runTime(file_.path(),file_.name());

  IOdictionary dict(
      IOobject
      (
          "",
          "",
          runTime.db(),
          IOobject::MUST_READ,
          IOobject::NO_WRITE
      )
  );

  // Delta
  scalar delta_Bc = dict.lookupOrDefault<scalar>("delta_Bc",0.0);
  scalar delta_Hw = dict.lookupOrDefault<scalar>("delta_Hw",0.0);

  // Switch
  Switch printPressure = dict.lookupOrDefault<Switch>("printPressure","yes");

  // Mutation++ mixture name
  word mixtureMutationPPName = dict.lookup("mixtureName");

  // Mutation++ mixture
  Mutation::Mixture mix= Mutation::Mixture(mixtureMutationPPName);

  int ne = mix.nElements();// Number of elements
  int ns = mix.nSpecies();  // Number of species

  int NT = readScalar(dict.lookup("nTemp"));  // Number of temperatures
  scalarList temperatures(NT);  // Temperatures
  scalar minT = readScalar(dict.lookup("minTemp")); // Minimum temperature
  scalar maxT = readScalar(dict.lookup("maxTemp")); // Maximum temperature


  // Computes the list of temperatures
  if (NT == 1) {
    temperatures[0]=minT;
  } else {
    for(int i = 0; i<NT; i++) {
      temperatures[i] = minT + i * (maxT-minT)/(NT-1);
    }
  }
  Info << "Temperatures: nTemp=" << NT << ", minTemp=" << minT << ", maxTemp=" << maxT << endl;

  scalarList BgList = dict.lookup("BgList");  // B'g list
  Info << "B'g = " << BgList << endl;

  scalar maxBc = readScalar(dict.lookup("maxBc"));
  scalar minBc = readScalar(dict.lookup("minBc"));

  if ( minBc < 0 ) {
    FatalError << "minBc must be more than 0" << exit(FatalError);
  }

  // Surface Composition based on different modes.
  List<Tuple2<word,scalar>> elemOnSurface;
  int mode_surf = dict.lookupOrDefault<int>("modeSurfComp",0);
  if (mode_surf == 0) {
    List<Tuple2<word,scalar>> tmp = dict.lookup("elemOnSurface");
    elemOnSurface = tmp; // list of the elements on the surface + mass fractions
    if ( elemOnSurface.size()==0) {
      FatalError << "elemOnSurface size must be >= 1." << exit(FatalError);
    }
    forAll(elemOnSurface, eI) {
      bool isElem = false;
      for (int i = 0; i < ne; i++) {
        if (elemOnSurface[eI].first()==mix.elementName(i)) {
          isElem = true;
          break;
        }
      }
      if (!isElem) {
        FatalError << elemOnSurface[eI].first() << " element not found in the mixture." << exit(FatalError);
      }
    }
  } else if (mode_surf == 1) {
    const scalar zC = readScalar(dict.lookup("massCSurface"));
    elemOnSurface = computeSiCSiO2SurfComp(mix, zC);
  } else {
    FatalError << "Unrecognised mode for imposing surface composition: " << mode_surf << exit(FatalError);
  }

  scalarList pressures = dict.lookup("pressures"); // Pressure [Pa]
  Info << "Pressures = " << pressures << endl;

  // Initialize the pyrolysis gas elemental composition
  Info << "Pyrolysis gas elemental composition:" << endl;
  Eigen::VectorXd z_pg(Eigen::VectorXd::Zero(ne));

  for (int i = 0; i < ne ; i++) {
    word elemName = mix.elementName(i);
    z_pg[i] = dict.lookupOrDefault<scalar>("z_pg["+elemName+"]",0.0);
    Info << "z_pg[" << elemName << "]=" << z_pg[i] << endl;
  }

  // Initialize the boundary layer edge elemental composition
  Info << "Boundary layer edge elemental composition:" << endl;
  Eigen::VectorXd z_e(Eigen::VectorXd::Zero(ne));

  for (int i = 0; i < ne ; i++) {
    word elemName = mix.elementName(i);
    z_e[i] = dict.lookupOrDefault<scalar>("z_e["+elemName+"]",0.0);
    Info << "z_e[" << elemName << "]=" << z_e[i] << endl;
  }


  double * Xw = new double[ns]; // Wall elemental mole composition
  for(int i=0; i<ns; i++) {
    Xw[i]=0;
  }
  double Bc = 0; // Dimensionless char ablation rate B'ca [-]
  double hw = 0; // Wall enthalpy
  Eigen::VectorXd X_sg(Eigen::VectorXd::Zero(ns)); // Solid, liquid and gaseous mole fraction
  List<List<scalarList>> dataBc(pressures.size()); // All the B'ca [pI][bI][tI]

  const word tableName(dict.lookup("tableName")); // name of the B' table
  OFstream os_table(tableName);
  const int numberBc= readScalar(dict.lookup("numberBc")); // number of B'c in the table
  scalar slopeMin = dict.lookupOrDefault<scalar>("slopeMin",1e-20); // B'c minimum slope


  // Loop over the B'g
  forAll(pressures, pI) {
    scalar p = pressures[pI]; // Temperature [K]

    word dirName = "p_" + name(p);
    if (printPressure) {
      if (!isDir(dirName)) {
        system("mkdir "+dirName);
        Info << "Creates \"" << dirName << "\" folder"<< endl;
      } else {
        FatalError << dirName << " directory exists already." << exit(FatalError);
      }
    } else {
      Info << "Compute " << dirName << endl;
    }

    PtrList<OFstream> os_out;
    PtrList<OFstream> os_out2;
    PtrList<OFstream> os_out3;
    if (printPressure) {
      os_out.resize(BgList.size());
      os_out2.resize(BgList.size());
      os_out3.resize(BgList.size());
    }
    dataBc[pI].resize(BgList.size());

    // Loop over the pressures
    forAll(BgList, bI) {
      scalar Bg = BgList[bI]; // DimLensionless pyrolysis gas blowing rate, B'g [-]

      if (printPressure) {
        // Output files for B'ca and species equilibrium composition
        os_out.set(bI, new OFstream(dirName+"/Bc_Bg"+name(Bg)));
        os_out[bI] << "//T[K] B'ca[-]" << endl;
        os_out2.set(bI, new OFstream(dirName+"/hw_Bg"+name(Bg)));
        os_out2[bI] << "//T[K] hw[J/kg]" << endl;
        os_out3.set(bI, new OFstream(dirName+"/species_Bg"+name(Bg)));
        os_out3[bI] << "//T[K] ";
        for (int i = 0; i < ns; i++) {
          os_out3[bI] << "X_" << mix.speciesName(i) << "[-] ";
        }
        os_out3[bI] << endl;
      }
      dataBc[pI][bI].resize(temperatures.size());

      // Loop over the temperatures
      forAll(temperatures,tI) {
        scalar T = temperatures[tI]; // Temperature [K]
        surfaceMassBalance(mix, elemOnSurface, z_e, z_pg, T, p, Bg, Bc, hw, X_sg, maxBc); // compute the surface mass balance

        if (Bc>=maxBc) {
          Bc=maxBc;
        }
        if (Bc < minBc || Bc < 0) {
          Bc=minBc;
        }

        if (printPressure) {
          // Write in the output files
          os_out[bI] << T << " " << Bc << endl;
          os_out2[bI] << T << " " << hw << endl;
          os_out3[bI] << T << " ";
          for (int i = 0; i < ns; i++) {
            os_out3[bI] << X_sg[i] << " ";
          }
          os_out3[bI] << endl;
        }

        // Store the B'c values
        dataBc[pI][bI][tI]=Bc;
      }
    }
  }


  // Write B'c table
  Info << "Writing " << tableName << endl;
  os_table.setf(ios_base::scientific, ios_base::floatfield);
  os_table.precision(5);
  int width=15;
  os_table << setf(ios_base::left) << setw(width) << "// p[Pa]";
  os_table << setf(ios_base::left) << setw(width) << "Bg[-]";
  os_table << setf(ios_base::left) << setw(width) << "Bc[-]";
  os_table << setf(ios_base::left) << setw(width) << "T[K]";
  os_table << setf(ios_base::left) << setw(width) << "hw[J/kg]";
  os_table << endl;

  // Loop over the pressures
  forAll(pressures, pI) {
    const double p = pressures[pI];
    // Loop over the B'g
    forAll(BgList, bI) {
      const double Bg = BgList[bI];
      // Loop over the temperatures
      forAll(temperatures, tI) {
        double& Bc = dataBc[pI][bI][tI];
        const double T = temperatures[tI];

        // Force the slope to be minimum slopeMin
        if (tI>0) {
          const double T_old = temperatures[tI-1];
          const double Bc_old = dataBc[pI][bI][tI-1];
          const double slope = (Bc-Bc_old)/(T-T_old);
          if (slope<slopeMin) {
            Bc = Bc_old + slopeMin*(T-T_old);
          }
        }

        // Compute the table points based on log(B'c)
        if ( Bc >= maxBc || tI == NT - 1) { // Last temperature or maximum B'c
          List<double> inverseTw_crop(tI+1);
          List<double> inverseBc_crop(tI+1);
          for(int i = 0 ; i <= tI ; i++) {
            inverseTw_crop[i]=temperatures[i];
            inverseBc_crop[i]=dataBc[pI][bI][i];
          }

          // Compute the first log(B'c) and the delta(log(B'c)) depending on the number of B'c (numberBc)
          double log_Bc0 = Foam::log(dataBc[pI][bI][0]);
          double deltaBc_log = (Foam::log(Bc) - log_Bc0 )/ (numberBc-1);

          for (int bJ=0; bJ<numberBc; bJ++) {
            double Bc_crop = Foam::exp(log_Bc0+deltaBc_log*bJ);
            double T_crop = linearInterpolation(Foam::log(inverseBc_crop), inverseTw_crop,Foam::log(Bc_crop));
            double hw_crop=0;

            // Compute the surface mass balance at the cropped temperature
            surfaceMassBalance(mix, elemOnSurface, z_e, z_pg, T_crop, p, Bg, Bc_crop, hw_crop, X_sg, maxBc);

            // Modifications
            Bc_crop+=delta_Bc;
            hw_crop+=delta_Hw;

            if (Bc_crop>=maxBc) {
              Bc_crop=maxBc;
            }
            if (Bc_crop < minBc || Bc_crop < 0) {
              Bc_crop =  minBc;
            }

            // Write the B' table
            os_table << setw(width) << p;
            os_table << setw(width) << Bg;
            os_table << setw(width) << Bc_crop;
            os_table << setw(width) << T_crop;
            os_table << setw(width) << hw_crop;
            os_table << endl;

            if (Bc_crop==maxBc) {
              break;
            }
          }
          break;
        }
      }
    }
  }

  return 0;
}

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
    along with OpenFOAM.  If Profilet, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ProfileIOModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ProfileIOModel::ProfileIOModel
(
    const fvMesh& mesh,
    const word& dictName
)
  :
noIOModel(mesh, dictName),
plot1DProfileList_(simpleIOModel::materialDict_.subDict("IO").lookup("plot1DProfileList")),
plot1DMassLoss_(simpleIOModel::materialDict_.subDict("IO").lookup("plot1DMassLoss")),
plot1DFolder_(simpleIOModel::materialDict_.subDict("IO").template lookupOrDefault<word>("plot1DFolder","output/"+dictName_+"/profile")),
rho_s_(meshLookupOrConstructScalar(mesh, "rho_s")),
V0_(0.0),
rho_s0_(0.0),
initOutput_(simpleIOModel::initOutput()),
createFolder_(createFolder()),
os_massLoss_(plot1DFolder_+"/massLoss"),
topPatch_(simpleIOModel::materialDict_.subDict("IO").lookup("topPatch")),
bottomPatch_(simpleIOModel::materialDict_.subDict("IO").lookup("bottomPatch"))
{
  // Verify fields exist
  foundFieldsInMesh(mesh_, plot1DProfileList_);

  // Mass loss init
  os_massLoss_ << "// t(s) "
               << "        "
               << "Mean density(kg/m³)"
               << "        "
               << "V/V0(-)"
               << "        "
               << "m/m0(-)"
               << "        "
               << endl;

  forAll(rho_s_, i) {
    V0_+=mesh.V()[i];
    rho_s0_+=mesh.V()[i]*rho_s_[i];
  }
  rho_s0_=rho_s0_/V0_;
  update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ProfileIOModel::~ProfileIOModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ProfileIOModel::update()
{
  noIOModel::writeOutput();

  if (mesh_.time().outputTime() || (mesh_.time().value() == mesh_.time().startTime().value())) {
    forAll(plot1DProfileList_, fieldI) {
      OFstream os_profile_(plot1DFolder_+"/"+plot1DProfileList_[fieldI]+"_"+mesh_.time().timeName());
      if (mesh_.objectRegistry::foundObject<volScalarField>(plot1DProfileList_[fieldI])) {
        volScalarField field_ = mesh_.objectRegistry::lookupObject<volScalarField>(plot1DProfileList_[fieldI]);
        os_profile_ << "// location(m) " << " " <<  plot1DProfileList_[fieldI] << "(" <<  field_.dimensions()  << ")"  << nl;
        scalarList x_;
        scalarList y_;
        scalarList z_;
        scalarList f_;

        label topPatchID = mesh_.boundaryMesh().findPatchID(topPatch_);
        label bottomPatchID = mesh_.boundaryMesh().findPatchID(bottomPatch_);
        x_.append(mesh_.Cf().boundaryField()[bottomPatchID][0].x());
        y_.append(mesh_.Cf().boundaryField()[bottomPatchID][0].y());
        z_.append(mesh_.Cf().boundaryField()[bottomPatchID][0].z());
        f_.append(field_.boundaryField()[bottomPatchID][0]);
        forAll(field_, fieldI) {
          x_.append(mesh_.C()[fieldI].x());
          y_.append(mesh_.C()[fieldI].y());
          z_.append(mesh_.C()[fieldI].z());
          f_.append(field_[fieldI]);
        }
        x_.append(mesh_.Cf().boundaryField()[topPatchID][0].x());
        y_.append(mesh_.Cf().boundaryField()[topPatchID][0].y());
        z_.append(mesh_.Cf().boundaryField()[topPatchID][0].z());
        f_.append(field_.boundaryField()[topPatchID][0]);

        word axis1D_ = "none";
        forAll(x_, fieldI) {
          if (fieldI > 0) {
            if (!assertEquals(x_[fieldI-1],x_[fieldI])) {
              axis1D_="x";
              break;
            }
          }
        }
        forAll(y_, fieldI) {
          if (fieldI > 0) {
            if (!assertEquals(y_[fieldI-1],y_[fieldI])) {
              if (axis1D_ == "x") {
                FatalErrorInFunction << "You can not use plot1DProfileList with a 3D mesh." << exit(FatalError);
              }
              axis1D_="y";
              break;
            }
          }
        }

        forAll(z_, fieldI) {
          if (fieldI > 0) {
            if (!assertEquals(z_[fieldI-1],z_[fieldI])) {
              if (axis1D_ == "x" || axis1D_ == "y_" ) {
                FatalErrorInFunction << "You can not use plot1DProfileList with a 3D mesh." << exit(FatalError);
              }
              axis1D_="z";
              break;
            }
          }
        }
        if (axis1D_ == "none") {
          FatalErrorInFunction << "The 1D axis is not found in plot1DProfile output." << exit(FatalError);
        }
        forAll(x_, fieldI) {
          if (axis1D_=="x") {
            os_profile_ <<  x_[fieldI] << " ";

          }
          if (axis1D_=="y") {
            os_profile_ <<  y_[fieldI] << " ";

          }
          if (axis1D_=="z") {
            os_profile_ <<  z_[fieldI] << " ";
          }
          os_profile_ << f_[fieldI] << endl;
        }
      }
      if (mesh_.objectRegistry::foundObject<volVectorField>(plot1DProfileList_[fieldI])) {
        FatalErrorInFunction << "volVectorField " << plot1DProfileList_[fieldI] << " not implemented" << exit(FatalError);
        volVectorField field_ = mesh_.objectRegistry::lookupObject<volVectorField>(plot1DProfileList_[fieldI]);
      }
      if (mesh_.objectRegistry::foundObject<volTensorField>(plot1DProfileList_[fieldI])) {
        FatalErrorInFunction << "volTensorField " << plot1DProfileList_[fieldI] << " not implemented" << exit(FatalError);
        volTensorField field_ = mesh_.objectRegistry::lookupObject<volTensorField>(plot1DProfileList_[fieldI]);
      }
      if (mesh_.objectRegistry::foundObject<surfaceScalarField>(plot1DProfileList_[fieldI])) {
        FatalErrorInFunction << "surfaceScalarField " << plot1DProfileList_[fieldI] << " not implemented" << exit(FatalError);
        surfaceScalarField field_ = mesh_.objectRegistry::lookupObject<surfaceScalarField>(plot1DProfileList_[fieldI]);
      }
      if (mesh_.objectRegistry::foundObject<surfaceVectorField>(plot1DProfileList_[fieldI])) {
        FatalErrorInFunction << "surfaceVectorField " << plot1DProfileList_[fieldI] << " not implemented" << exit(FatalError);
        surfaceVectorField field_ = mesh_.objectRegistry::lookupObject<surfaceVectorField>(plot1DProfileList_[fieldI]);
      }
      if (mesh_.objectRegistry::foundObject<surfaceTensorField>(plot1DProfileList_[fieldI])) {
        FatalErrorInFunction << "surfaceTensorField " << plot1DProfileList_[fieldI] << " not implemented" << exit(FatalError);
        surfaceTensorField field_ = mesh_.objectRegistry::lookupObject<surfaceTensorField>(plot1DProfileList_[fieldI]);
      }
    }

    if(plot1DMassLoss_) {
      scalar VtimeName = 0;
      scalar rho_s_mean = 0;
      forAll(rho_s_, i) {
        VtimeName+=mesh_.V()[i];
        rho_s_mean+=mesh_.V()[i]*rho_s_[i];
      }
      rho_s_mean = rho_s_mean / VtimeName;
      os_massLoss_ <<  mesh_.time().timeName() << "	"
                   << rho_s_mean << "	" << VtimeName/V0_ << "	" << rho_s_mean/rho_s0_ << endl;
    }
  }
}

Switch Foam::ProfileIOModel::createFolder()
{
  system("mkdir -p " + plot1DFolder_);
  return true;
}

// ************************************************************************* //

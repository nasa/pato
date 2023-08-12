/*---------------------------------------------------------------------------*\

Notices:

    Copyright © 2010 United States Government as represented by the
    Administrator of the National Aeronautics and Space Administration.  All
    Rights Reserved.

Disclaimers

    No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY
    OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
    LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
    SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
    PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE
    SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION,
    IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES
    NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY
    PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE
    PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT
    SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND
    LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL
    SOFTWARE, AND DISTRIBUTES IT "AS IS."

    Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
    AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS,
    AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT
    SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES
    ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR
    RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL
    INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS
    AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT
    PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE
    THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.

  ---------------------------------------------------------------------------

Application

    inletPower

Description

    Post-processor to calculate inlet power

Usage
    Run inletPower after calculation in order to obtain the power on specified
    field. Reference temperature, patch and temperature field name should be
    provided. You should provide the region name. You can select which time
    step you want to post process.
    inletPower [options] <patch>

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "fvCFD.H"
#include "PATOx.H"
#include "IOFunctions.H"
#include "OFstream.H"
#include "IOobject.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  argList::noParallel();
  //- Patch name
  argList::validArgs.append("patch");
  //- Initial temperature
  argList::addOption
  (
      "Tref",
      "T0",
      "Reference temperature for power calculation (default : 0 K)"

  );
  //- Temperature field name (mainly Tg, Ta or T in PATO)
  argList::addOption
  (
      "Tfield",
      "Tname",
      "Temperature field name (default : Tg)"

  );
  //- Field name for velocity (default : "U")
  argList::addOption
  (
      "field",
      "fieldName",
      "Velocity field name (default : U)"
  );
  timeSelector::addOptions();
#include "addRegionOption.H"
#include "setRootCase.H"
  scalar T0_ = 0;
  args.optionReadIfPresent("Tref", T0_);
  word Tfield = "Tg";
  args.optionReadIfPresent("Tfield", Tfield);
  word field = "U";
  args.optionReadIfPresent("field", field);
  word region = "";
  args.optionReadIfPresent("region", region);
  if (!args.check()) {
    FatalError.exit();
  }
  if (region.empty()) {
    FatalErrorInFunction<< "You need to specify a region name."
                        << exit(FatalError);
  }

  // creating output directory
  if (!isDir("output/"+region)) {
    mkDir("output/"+region);
  }
  OFstream output("output/"+region+"/inputpower");

#include "createTime.H"
  instantList timeDirs = timeSelector::select0(runTime, args);
  word patchName_ = args[1];
  fvMesh mesh
  (
      IOobject
      (
          region,
          runTime.timeName(),
          runTime,
          IOobject::MUST_READ
      )
  );

  // Create models
  simpleGasPropertiesModel& gasPropertiesModel_
  (
      meshLookupOrConstructModel<simpleGasPropertiesModel>
      (
          mesh,
          region,
          simpleGasPropertiesModel::modelName
      )
  );
  simpleEnergyModel& energyModel_
  (
      meshLookupOrConstructModel<simpleEnergyModel>
      (
          mesh,
          region,
          simpleEnergyModel::modelName
      )
  );
  simpleMassModel& massModel_
  (
      meshLookupOrConstructModel<simpleMassModel>
      (
          mesh,
          region,
          simpleMassModel::modelName
      )
  );

  // References to fields
  volScalarField& T_(energyModel_.refVolField<scalar>(Tfield));
  volScalarField& p_(massModel_.refVolField<scalar>("p"));
  volVectorField& U_(massModel_.refVolField<vector>("U"));
  volScalarField& cp_(gasPropertiesModel_.refVolField<scalar>("cp_g"));
  volScalarField& rho_(gasPropertiesModel_.refVolField<scalar>("rho_g"));

  const label patchID = mesh.boundaryMesh().findPatchID(patchName_);
  scalar power_ = 0;

  Info<< "Patch : " << patchName_ << endl;
  Info<< "Initial temperature Ti = " << T0_ << endl;
  Info<< "Temperature field name : " << Tfield << endl;
  Info<< "Field U : "<< field << endl << nl;

  // Write Header of output file
  output << "// time (s)" << " P inlet (W)" << endl;

  forAll(timeDirs, timei) {
    // Update time, mesh and properties
    runTime.setTime(timeDirs[timei], timei);
    Info<< "Time = " << runTime.timeName() << endl;

    // Check the existence of 0/U. Skip if it does not exist.
    if (runTime.value() == 0 && !isFile("0/"+region+"/"+field)) {
      WarningInFunction
          << "Velocity field : " << field
          << " does not exist for time 0. " << nl
          << "    Skipping."
          << endl;
      continue;
    }

    volScalarField T_read
    (
        IOobject
        (
            Tfield.c_str(),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    T_.ref()=T_read();
    T_.boundaryFieldRef()==T_read.boundaryField();

    volScalarField p_read
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    p_.ref()=p_read();
    p_.boundaryFieldRef()==p_read.boundaryField();

    volVectorField U_read
    (
        IOobject
        (
            field.c_str(),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    U_.ref()=U_read();
    U_.boundaryFieldRef()==U_read.boundaryField();

    gasPropertiesModel_.update();
    scalarField magSf(mesh.magSf().boundaryField()[patchID]);
    scalarField cpp(cp_.boundaryField()[patchID]);
    scalarField rhop(rho_.boundaryField()[patchID]);
    vectorField Up(U_.boundaryField()[patchID]);
    scalarField magU(mag(U_.boundaryField()[patchID]));

    // Calculate temperature difference
    forAll(T_, celli) {
      T_[celli] -= T0_;
    }
    forAll(mesh.boundaryMesh(), patchi) {
      forAll(T_.boundaryFieldRef()[patchi], facei) {
        T_.boundaryFieldRef()[patchi][facei] -= T0_;
      }
    }

    scalarField dT(T_.boundaryField()[patchID]);

    // Calculate power on given patch
    power_ = sum(magSf*cpp*rhop*magU*dT);

    // Write output
    Info<< "Inlet power = " << power_ << " W" << endl;
    output<< runTime.timeName() << "           " << power_ << endl;
  }
}

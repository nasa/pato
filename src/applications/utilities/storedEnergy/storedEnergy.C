/*---------------------------------------------------------------------------*\
 *
 N otices:
 *

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

 storedEnergy

 Description

 Post-processor to calculate stored energy in the solid phase.

 Usage
 Run storedEnergy after calculation in order to obtain the stored energy on
solid phase. Reference temperature, temperature field name and region name
should be  provided. You can select which time step you want to post process.
 storedEnergy [options]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "fvCFD.H"
#include "PATOx.H"
#include "IOFunctions.H"
#include "OFstream.H"
#include "IOobject.H"


#include <inttypes.h>


using namespace Foam;

int main(int argc, char *argv[])
{
  //- Temperature field name (mainly Ts, Ta or T in PATO)
  argList::addOption
  (
      "Tfield",
      "Tname",
      "Temperature field name (default : Ta)"

  );
  timeSelector::addOptions();
#include "addRegionOption.H"
#include "setRootCase.H"
  word Tfield = "Ta";
  args.optionReadIfPresent("Tfield", Tfield);
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
  if (!isDir("output/"+region) && Pstream::master()) {
    mkDir("output/"+region);
  }
  OFstream output("output/"+region+"/storedEnergy");

#include "createTime.H"
  instantList timeDirs = timeSelector::select0(runTime, args);
  scalar energy_ = 0;
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

  runTime.setTime(runTime.startTime(), runTime.startTimeIndex());

  // Create the models
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
  simplePyrolysisModel& pyrolysisModel_
  (
      meshLookupOrConstructModel<simplePyrolysisModel>
      (
          mesh,
          region,
          simplePyrolysisModel::modelName
      )
  );
  simpleMaterialPropertiesModel& materialPropertiesModel_
  (
      meshLookupOrConstructModel<simpleMaterialPropertiesModel>
      (
          mesh,
          region,
          simpleMaterialPropertiesModel::modelName
      )
  );
  materialPropertiesModel_.initialize();

  // References to fields
  volScalarField T0_(energyModel_.refVolField<scalar>(Tfield));
  volScalarField& T_(energyModel_.refVolField<scalar>(Tfield));
  volScalarField& cp_(energyModel_.refVolField<scalar>("cp"));
  volScalarField& rho_(energyModel_.refVolField<scalar>("rho_s"));
  volScalarField& p_(massModel_.refVolField<scalar>("p"));

  Info<< "Temperature field name : " << Tfield << endl;

  // Write Header of output file
  if (Pstream::master()) {
    output << "// time (s)" << " E stored (J)" << endl;
  }

  forAll(timeDirs, timei) {
    // Update time, mesh and properties
    runTime.setTime(timeDirs[timei], timei);
    Info<< "Time = " << runTime.timeName() << endl;

    // Read and transfer variables
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

    materialPropertiesModel_.update();

    // Calculate temperature difference
    T_ -= T0_;

    energy_ = gSum(rho_*cp_*T_*mesh.V());

    // Write output
    if (Pstream::master()) {
      Info<< "Stored energy = " << energy_ << " J" << endl;
      output<< runTime.timeName() << "           " << energy_ << endl;
    }
  }

}

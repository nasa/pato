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

    compressibleDarcyForchheimer

Description

    Transient solver for compressible flows in porous media using
    Darcy-Forchheimer equation.
    Implict resolution in Pressure, anisotropic permeability.
    Explicit treatment of Forchheimer term.

\*---------------------------------------------------------------------------*/

// Linking librairies (class & function containers) needed by the program
#include "fvCFD.H"

// A C++ program shall contain a global function named main,
// which is the designated start of the program.
int main(int argc, char *argv[])
{

  // executing case initializations
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"

  // Fill Gamma Field
  Gamma.replace
  (
      0,
      K.component(Tensor<double>::XX)
      /(
          mu + Beta.component(Tensor<double>::XX)*p*M/(R*T)*magV
          *K.component(Tensor<double>::XX)*eps
      )
  );
  Gamma.replace
  (
      1,
      K.component(Tensor<double>::XY)
      /(
          mu + Beta.component(Tensor<double>::XY)*p*M/(R*T)*magV
          *K.component(Tensor<double>::XY)*eps
      )
  );
  Gamma.replace
  (
      2,
      K.component(Tensor<double>::XZ)
      /(
          mu + Beta.component(Tensor<double>::XZ)*p*M/(R*T)*magV
          *K.component(Tensor<double>::XZ)*eps
      )
  );
  Gamma.replace
  (
      3,
      K.component(Tensor<double>::YX)
      /(
          mu + Beta.component(Tensor<double>::YX)*p*M/(R*T)*magV
          *K.component(Tensor<double>::YX)*eps
      )
  );
  Gamma.replace
  (
      4,
      K.component(Tensor<double>::YY)
      /(
          mu + Beta.component(Tensor<double>::YY)*p*M/(R*T)*magV
          *K.component(Tensor<double>::YY)*eps
      )
  );
  Gamma.replace
  (
      5,
      K.component(Tensor<double>::YZ)
      /(
          mu + Beta.component(Tensor<double>::YZ)*p*M/(R*T)*magV
          *K.component(Tensor<double>::YZ)*eps
      )
  );
  Gamma.replace
  (
      6,
      K.component(Tensor<double>::ZX)
      /(
          mu + Beta.component(Tensor<double>::ZX)*p*M/(R*T)*magV
          *K.component(Tensor<double>::ZX)*eps
      )
  );
  Gamma.replace
  (
      7,
      K.component(Tensor<double>::ZY)
      /(
          mu + Beta.component(Tensor<double>::ZY)*p*M/(R*T)*magV
          *K.component(Tensor<double>::ZZ)*eps
      )
  );
  Gamma.replace
  (
      8,
      K.component(Tensor<double>::ZZ)
      /(
          mu + Beta.component(Tensor<double>::ZZ)*p*M/(R*T)*magV
          *K.component(Tensor<double>::ZZ)*eps
      )
  );
  Gamma.field() *= p.ref()*M.value()/(R.value()*T.value());
  Gamma.boundaryFieldRef() *=
      p.boundaryFieldRef()*M.value()/(R.value()*T.value());

  // solving over the mesh
  while (runTime.loop()) {
    Info<< "Time = " << runTime.timeName() << nl << endl;

    // semi-implicit resolution of the pressure field
    solve(fvm::ddt(Eta, p) - fvm::laplacian(Gamma, p));

    // Partial update of tensor Gamma
    // This is the velocity part
    Gamma.replace
    (
        0,
        K.component(Tensor<double>::XX)
        /(
            mu + Beta.component(Tensor<double>::XX)*p*M/(R*T)*magV
            *K.component(Tensor<double>::XX)*eps
        )
    );
    Gamma.replace
    (
        1,
        K.component(Tensor<double>::XY)
        /(
            mu + Beta.component(Tensor<double>::XY)*p*M/(R*T)*magV
            *K.component(Tensor<double>::XY)*eps
        )
    );
    Gamma.replace
    (
        2,
        K.component(Tensor<double>::XZ)
        /(
            mu + Beta.component(Tensor<double>::XZ)*p*M/(R*T)*magV
            *K.component(Tensor<double>::XZ)*eps
        )
    );
    Gamma.replace
    (
        3,
        K.component(Tensor<double>::YX)
        /(
            mu + Beta.component(Tensor<double>::YX)*p*M/(R*T)*magV
            *K.component(Tensor<double>::YX)*eps
        )
    );
    Gamma.replace
    (
        4,
        K.component(Tensor<double>::YY)
        /(
            mu + Beta.component(Tensor<double>::YY)*p*M/(R*T)*magV
            *K.component(Tensor<double>::YY)*eps
        )
    );
    Gamma.replace
    (
        5,
        K.component(Tensor<double>::YZ)
        /(
            mu + Beta.component(Tensor<double>::YZ)*p*M/(R*T)*magV
            *K.component(Tensor<double>::YZ)*eps
        )
    );
    Gamma.replace
    (
        6,
        K.component(Tensor<double>::ZX)
        /(
            mu + Beta.component(Tensor<double>::ZX)*p*M/(R*T)*magV
            *K.component(Tensor<double>::ZX)*eps
        )
    );
    Gamma.replace
    (
        7,
        K.component(Tensor<double>::ZY)
        /(
            mu + Beta.component(Tensor<double>::ZY)*p*M/(R*T)*magV
            *K.component(Tensor<double>::ZZ)*eps
        )
    );
    Gamma.replace
    (
        8,
        K.component(Tensor<double>::ZZ)
        /(
            mu + Beta.component(Tensor<double>::ZZ)*p*M/(R*T)*magV
            *K.component(Tensor<double>::ZZ)*eps
        )
    );

    // Update tensor Gamma
    Gamma.field() *= p.ref()*M.value()/(R.value()*T.value());
    Gamma.boundaryFieldRef() *=
        p.boundaryFieldRef()*M.value()/(R.value()*T.value());

    // compute average gas velocity (Darcy's law) and gas flow rate
    magV=mag(v);
    v = - ( (Gamma*R*T/(p*M)) & fvc::grad(p))*1/eps;
    q = v*p*M/(R*T);

    runTime.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
  }

  // exiting
  Info<< "End\n" << endl;

  return 0;
}


// ************************************************************************* //

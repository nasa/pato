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

    buoyantDarcy2T

Description

    Transient compressible, buoyant Darcy solver for porous media.
    Implict resolution in Temperature, anisotropic thermal conductivity.

\*---------------------------------------------------------------------------*/

// Linking libraries (class & function containers) needed by the program
#include "fvCFD.H"
#include "PATOx.H"
#include "fvOptions.H"

// A C++ program shall contain a global function named main,
// which is the designated start of the program.
int main(int argc, char *argv[])
{

  // executing case initializations
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"

  // solving over the mesh
  while (runTime.loop()) {
    Info<< "Time = " << runTime.timeName() << nl << endl;

    // update scalar Beta
    Beta = eps_g*M/(R*Tg);

    // update tensor Gamma
    Gamma = rho_g*K /mu;

    // update tensor Theta
    Theta = Gamma*M/(R*Tg);

    // Gravitational mass flow rate
    phiG = linearInterpolate(Theta & g) & mesh.Sf();

    // Semi-implicit calculation of the pressure field
    fvScalarMatrix pEqn
    (
        fvm::ddt(Beta, p) - fvm::laplacian(Gamma, p) + fvm::div(phiG, p)
    );

    solve(pEqn);

    // print min and max pressures in the terminal during the run
    Info << "p_max = p"
         << max(p).value()
         << " / p_min = "
         << min(p).value()
         << endl;

    // average gas velocity (Darcy's law) from the new pressure gradient (cell centered)
    rho_g = p * M / (R * Tg);
    rho_g.correctBoundaryConditions();

    U = -(K/mu & (fvc::grad(p) - rho_g*g));
    U.correctBoundaryConditions();

    vG = U/eps_g;
    vG.correctBoundaryConditions();

    // Mass flow rate
    mDotG = U*rho_g;
    phi_g = linearInterpolate(mDotG) & mesh.Sf();

    // Calculation of gas temperature
    fvScalarMatrix TgEqn
    (
        eps_g*rho_g*cp_g*fvm::ddt(Tg) 			// Storage
        + eps_g*cp_g*fvm::div(phi_g, Tg)              	// Convection
        - fvm::laplacian(k_g, Tg) 				// Conduction
        + fvm::Sp(Hv0, Tg) - Hv0*Ta				// Heat exchange
        ==
        fvOptions(rho_g*cp_g, Tg)                           // Source term
    );

    fvOptions.constrain(TgEqn);

    TgEqn.solve();

    fvOptions.correct(Tg);

    // Calculation of solid temperature
    fvScalarMatrix TaEqn
    (
        (1. - eps_g)*rho_s*cp_s*fvm::ddt(Ta) 		// Storage
        - fvm::laplacian(k_a, Ta) 				// Conduction
        + Hv0*(Ta - Tg)      				        // Heat exchange
        ==
        fvOptions(rho_s*cp_s, Ta)                           // Source term
    );

    fvOptions.constrain(TaEqn);

    TaEqn.solve();

    fvOptions.correct(Ta);

    // update tensor Theta
    Theta = Gamma*p*M/(R*Tg);

    // Gravitational mass flow rate
    phiG = linearInterpolate(Theta & g) & mesh.Sf();

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

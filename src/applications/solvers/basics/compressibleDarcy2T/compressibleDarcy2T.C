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

    compressibleDarcy2T

Description 

    Transient heat transfer solver for porous media.
    Implict resolution in Temperature, anisotropic thermal conductivity.

\*---------------------------------------------------------------------------*/

// Linking librairies (class & function containers) needed by the program
#include "fvCFD.H"

// A C++ program shall contain a global function named main, which is the designated start of the program.
int main(int argc, char *argv[])
{

    // executing case initializations
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // solving over the mesh
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

// update scalar Beta
	Beta = eps * M / (R * T_g);

// update tensor Gamma
	Gamma = p * M / (mu * R * T_g) * K;

// Semi-implicit calculation of the pressure field
	solve
	(
    		fvm::ddt(Beta, p) - fvm::laplacian(Gamma, p)
	);

	T_g.correctBoundaryConditions();

// Update average gas velocity (Darcy's law), gas density and gas flow rate
	u_g = - (K & fvc::grad(p))/(eps*mu);
	rho_g = p * M / (R * T_g);
	q = u_g * rho_g;

// Update scalar phi_g
	phi_g = linearInterpolate(q) & mesh.Sf();

// Calculation of solid temperature
	solve
	(
		(1. - eps) * rho_s * cp_s * fvm::ddt(T_s) 		// Storage
		- fvm::laplacian(k_s, T_s) 				// Conduction
		+ fvm::Sp(hgs, T_s) - hgs * T_g				// Heat exchange
	);

// Calculation of gas temperature
	solve
	(
                eps * rho_g * cp_g * fvm::ddt(T_g) 			// Storage
		+ eps * cp_g * fvm::div(phi_g, T_g)			// Convection
		- fvm::laplacian(k_g, T_g) 				// Conduction
		+ fvm::Sp(hsg, T_g) - hsg * T_s				// Heat exchange
	);


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

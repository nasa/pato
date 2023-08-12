/* TEST VERSION*/

// Do any mesh changes
mesh.update();
// Make the fluxes relative to the mesh-motion
fvc::makeRelative(phi, rho, U);

/*
// Store momentum to set rhoUf for introduced faces.
            volVectorField rhoU("rhoU", rho*U);
	    surfaceVectorField rhoUf = fvc::interpolate(rho*U);

            // Do any mesh changes
            mesh.update();

            if (mesh.changing() && correctPhi)
            {
                // Calculate absolute flux from the mapped surface velocity
                phi = mesh.Sf() & rhoUf;

                #include "correctPhi.H"

                // Make the fluxes relative to the mesh-motion
                fvc::makeRelative(phi, rho, U);
            }
        }

        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }
*/

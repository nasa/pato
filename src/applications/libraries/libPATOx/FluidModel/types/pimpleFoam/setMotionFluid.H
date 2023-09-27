/* TEST VERSION*/

// Do any mesh changes
if (oCorr == 0 || moveMeshOuterCorrectors)
{
  mesh.update();
  if (mesh.changing()) {
    MRF.update();

    if (correctPhi) {
      // Calculate absolute flux
      // from the mapped surface velocity
      phi = mesh.Sf() & Uf();

#include "correctPhi.H"

#include "continuityErrors.H"

      // Make the flux relative to the mesh motion
      fvc::makeRelative(phi, U);
    }
    /*  TODO : Add checkMeshCourantNo */
  }
}

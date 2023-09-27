{
  // Update multi region coupled criteria

  Tw_n = T.boundaryField()[mesh.boundaryMesh().findPatchID(stitchingFluidPatchName)];
  Tw_nless1 = T.oldTime().boundaryField()[mesh.boundaryMesh().findPatchID(stitchingFluidPatchName)];

  //if (ii == 0) {
  //ii = 1;
  //Tw_nless1 = Tw_n*0.1;
  //}

  deltaTw = mag( (Tw_n - Tw_nless1) / (mesh.time()).deltaT().value() );
  //Tw_nless1 = Tw_n;
  maxDeltaTw_ = max(deltaTw);
}

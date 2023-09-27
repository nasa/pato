word historyPatchName_tmp(mesh.solutionDict().subDict("solidMechanics").lookup("historyPatch"));
historyPatchName = historyPatchName_tmp;
historyPatchID = mesh.boundaryMesh().findPatchID(historyPatchName);
if(historyPatchID == -1)
{
  Warning << "history patch " << historyPatchName
          << " not found. Force-displacement will not be written"
          << endl;
} else if(Pstream::master())
{
  Info << "Force-displacement for patch " << historyPatchName
       << " will be written to forceDisp.dat"
       << endl;
  mkDir("history");
  filePtr = new OFstream("history"/historyPatchName+"_forceDisp.dat");
  OFstream& forceDispFile = *filePtr;
  forceDispFile << "#Disp(mm)\tForce(N)" << endl;
}

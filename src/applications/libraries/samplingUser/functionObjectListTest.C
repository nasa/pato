/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "functionObjectListTest.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "argList.H"
#include "timeControlFunctionObject.H"
#include "IFstream.H"
#include "dictionaryEntry.H"
#include "stringOps.H"
#include "Tuple2.H"
#include "etcFiles.H"
#include "fileOperation.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

Foam::fileName Foam::functionObjectListTest::functionObjectDictPath
(
    "caseDicts/postProcessing"
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::functionObjectTest* Foam::functionObjectListTest::remove
(
    const word& key,
    label& oldIndex
)
{
  functionObjectTest* ptr = 0;

  // Find index of existing functionObject
  HashTable<label>::iterator fnd = indices_.find(key);

  if (fnd != indices_.end()) {
    oldIndex = fnd();

    // Retrieve the pointer and remove it from the old list
    ptr = this->set(oldIndex, 0).ptr();
    indices_.erase(fnd);
  } else {
    oldIndex = -1;
  }

  return ptr;
}


void Foam::functionObjectListTest::listDir
(
    const fileName& dir,
    HashSet<word>& foMap
)
{
  // Search specified directory for functionObject configuration files
  {
    fileNameList foFiles(readDir(dir));
    forAll(foFiles, f) {
      if (foFiles[f].ext().empty()) {
        foMap.insert(foFiles[f]);
      }
    }
  }

  // Recurse into sub-directories
  {
    fileNameList foDirs(readDir(dir, fileType::directory));
    forAll(foDirs, fd) {
      listDir(dir/foDirs[fd], foMap);
    }
  }
}


void Foam::functionObjectListTest::list()
{
  HashSet<word> foMap;

  fileNameList etcDirs(findEtcDirs(functionObjectDictPath));

  forAll(etcDirs, ed) {
    listDir(etcDirs[ed], foMap);
  }

  Info<< nl
      << "Available configured functionObjects:"
      << foMap.sortedToc()
      << nl;
}


Foam::fileName Foam::functionObjectListTest::findDict(const word& funcName)
{
  // First check if there is a functionObject dictionary file in the
  // case system directory
  fileName dictFile = stringOps::expand("$FOAM_CASE")/"system"/funcName;

  if (isFile(dictFile)) {
    return dictFile;
  } else {
    fileNameList etcDirs(findEtcDirs(functionObjectDictPath));

    forAll(etcDirs, i) {
      dictFile = search(funcName, etcDirs[i]);
      if (!dictFile.empty()) {
        return dictFile;
      }
    }
  }

  return fileName::null;
}


bool Foam::functionObjectListTest::readFunctionObject
(
    const string& funcNameArgs,
    dictionary& functionsDict,
    HashSet<word>& requiredFields,
    const word& region
)
{
  // Parse the optional functionObject arguments:
  //     'Q(U)' -> funcName = Q; args = (U); field = U
  //
  // Supports named arguments:
  //     'patchAverage(patch=inlet, p)' -> funcName = patchAverage;
  //         args = (patch=inlet, p); field = p

  word funcName(funcNameArgs);

  int argLevel = 0;
  wordList args;

  List<Tuple2<word, string>> namedArgs;
  bool namedArg = false;
  word argName;

  word::size_type start = 0;
  word::size_type i = 0;

  for
  (
      word::const_iterator iter = funcNameArgs.begin();
      iter != funcNameArgs.end();
      ++iter
  ) {
    char c = *iter;

    if (c == '(') {
      if (argLevel == 0) {
        funcName = funcNameArgs(start, i - start);
        start = i+1;
      }
      ++argLevel;
    } else if (c == ',' || c == ')') {
      if (argLevel == 1) {
        if (namedArg) {
          namedArgs.append
          (
              Tuple2<word, string>
              (
                  argName,
                  funcNameArgs(start, i - start)
              )
          );
          namedArg = false;
        } else {
          args.append
          (
              string::validate<word>(funcNameArgs(start, i - start))
          );
        }
        start = i+1;
      }

      if (c == ')') {
        if (argLevel == 1) {
          break;
        }
        --argLevel;
      }
    } else if (c == '=') {
      argName = string::validate<word>(funcNameArgs(start, i - start));
      start = i+1;
      namedArg = true;
    }

    ++i;
  }

  // Search for the functionObject dictionary
  fileName path = findDict(funcName);

  if (path == fileName::null) {
    WarningInFunction
        << "Cannot find functionObject file " << funcName << endl;
    return false;
  }

  // Read the functionObject dictionary
  IFstream fileStream(path);
  dictionary funcsDict(fileStream);
  dictionary* funcDictPtr = &funcsDict;

  if (funcsDict.found(funcName) && funcsDict.isDict(funcName)) {
    funcDictPtr = &funcsDict.subDict(funcName);
  }

  dictionary& funcDict = *funcDictPtr;

  // Insert the 'field' and/or 'fields' entry corresponding to the optional
  // arguments or read the 'field' or 'fields' entry and add the required
  // fields to requiredFields
  if (args.size() == 1) {
    funcDict.set("field", args[0]);
    funcDict.set("fields", args);
    requiredFields.insert(args[0]);
  } else if (args.size() > 1) {
    funcDict.set("fields", args);
    requiredFields.insert(args);
  } else if (funcDict.found("field")) {
    requiredFields.insert(word(funcDict.lookup("field")));
  } else if (funcDict.found("fields")) {
    requiredFields.insert(wordList(funcDict.lookup("fields")));
  }

  // Insert named arguments
  forAll(namedArgs, i) {
    IStringStream entryStream
    (
        namedArgs[i].first() + ' ' + namedArgs[i].second() + ';'
    );
    funcDict.set(entry::New(entryStream).ptr());
  }

  // Insert the region name if specified
  if (region != word::null) {
    funcDict.set("region", region);
  }

  // Merge this functionObject dictionary into functionsDict
  dictionary funcArgsDict;
  funcArgsDict.add(string::validate<word>(funcNameArgs), funcDict);
  functionsDict.merge(funcArgsDict);

  return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectListTest::functionObjectListTest
(
    const Time& t,
    const bool execution
)
  :
PtrList<functionObjectTest>(),
digests_(),
indices_(),
time_(t),
parentDict_(t.controlDict()),
execution_(execution),
updated_(false)
{
}


Foam::functionObjectListTest::functionObjectListTest
(
    const Time& t,
    const dictionary& parentDict,
    const bool execution
)
  :
PtrList<functionObjectTest>(),
digests_(),
indices_(),
time_(t),
parentDict_(parentDict),
execution_(execution),
updated_(false)
{
}


Foam::autoPtr<Foam::functionObjectListTest> Foam::functionObjectListTest::New
(
    const argList& args,
    const Time& runTime,
    dictionary& controlDict,
    HashSet<word>& requiredFields
)
{

  autoPtr<functionObjectListTest> functionsPtr;

  controlDict.add
  (
      dictionaryEntry("functions", controlDict, dictionary::null)
  );

  dictionary& functionsDict = controlDict.subDict("functions");

  word region = word::null;

  // Set the region name if specified
  if (args.optionFound("region")) {
    region = args["region"];
  }

  if
  (
      args.optionFound("dict")
      || args.optionFound("func")
      || args.optionFound("funcs")
  ) {
    if (args.optionFound("dict")) {
      controlDict.merge
      (
          IOdictionary
          (
              IOobject
              (
                  args["dict"],
                  runTime,
                  IOobject::MUST_READ_IF_MODIFIED
              )
          )
      );
    }

    if (args.optionFound("func")) {
      readFunctionObject
      (
          args["func"],
          functionsDict,
          requiredFields,
          region
      );

    }

    if (args.optionFound("funcs")) {
      wordList funcs(args.optionLookup("funcs")());

      forAll(funcs, i) {
        readFunctionObject
        (
            funcs[i],
            functionsDict,
            requiredFields,
            region
        );
      }
    }


    functionsPtr.reset(new functionObjectListTest(runTime, controlDict));
  } else {
    functionsPtr.reset(new functionObjectListTest(runTime));
  }

  functionsPtr->start();

  return functionsPtr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjectListTest::~functionObjectListTest()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjectListTest::clear()
{
  PtrList<functionObjectTest>::clear();
  digests_.clear();
  indices_.clear();
  updated_ = false;
}


Foam::label Foam::functionObjectListTest::findObjectID(const word& name) const
{
  forAll(*this, objectI) {
    if (operator[](objectI).name() == name) {
      return objectI;
    }
  }

  return -1;
}


void Foam::functionObjectListTest::on()
{
  execution_ = true;
}


void Foam::functionObjectListTest::off()
{
  // For safety, also force a read() when execution is turned back on
  updated_ = execution_ = false;
}


bool Foam::functionObjectListTest::status() const
{
  return execution_;
}


bool Foam::functionObjectListTest::start()
{
  return read();
}


bool Foam::functionObjectListTest::execute()
{
  bool ok = true;


  if (execution_) {
    if (!updated_) {
      read();
    }

    forAll(*this, objectI) {
      ok = operator[](objectI).execute() && ok;
      //ok = operator[](objectI).write() && ok;
    }
  }

  return ok;
}

Foam::List<Foam::List<Foam::Field<Foam::scalar> > > Foam::functionObjectListTest::executeScaling()
{
  bool ok = true;
  if (execution_) {
    if (!updated_) {
      read();
    }

    forAll(*this, objectI) {
      ok = operator[](objectI).execute() && ok;
      return operator[](objectI).scalarMasterFields_;
    }
  }
  List<List<Field<scalar> > > nullMasterFields;
  return nullMasterFields;
}

void Foam::functionObjectListTest::executeScalingAll
(
    autoPtr<List<List<Field<scalar> > > >& scalarMasterFields_ref,
    autoPtr<List<List<Field<vector> > > >& vectorMasterFields_ref,
    autoPtr<List<List<Field<sphericalTensor> > > >& sphericalTensorMasterFields_ref,
    autoPtr<List<List<Field<symmTensor> > > >& symmTensorMasterFields_ref,
    autoPtr<List<List<Field<tensor> > > >& tensorMasterFields_ref
)
{
  bool ok = true;
  if (execution_) {
    if (!updated_) {
      read();
    }

    forAll(*this, objectI) {
      ok = operator[](objectI).execute() && ok;
      scalarMasterFields_ref.reset(new List<List<Field<scalar> > >(operator[](objectI).scalarMasterFields_));
      vectorMasterFields_ref.reset(new List<List<Field<vector> > >(operator[](objectI).vectorMasterFields_));
      sphericalTensorMasterFields_ref.reset(new List<List<Field<sphericalTensor> > >(operator[](objectI).sphericalTensorMasterFields_));
      symmTensorMasterFields_ref.reset(new List<List<Field<symmTensor> > >(operator[](objectI).symmTensorMasterFields_));
      tensorMasterFields_ref.reset(new List<List<Field<tensor> > >(operator[](objectI).tensorMasterFields_));

    }
  }
}


bool Foam::functionObjectListTest::end()
{
  bool ok = true;

  if (execution_) {
    if (!updated_) {
      read();
    }

    forAll(*this, objectI) {
      ok = operator[](objectI).end() && ok;
    }
  }

  return ok;
}


bool Foam::functionObjectListTest::adjustTimeStep()
{
  bool ok = true;

  if (execution_) {
    if (!updated_) {
      read();
    }

    forAll(*this, objectI) {
      ok = operator[](objectI).adjustTimeStep() && ok;
    }
  }

  return ok;
}


bool Foam::functionObjectListTest::read()
{
  bool ok = true;
  updated_ = execution_;

  // Avoid reading/initializing if execution is off
  if (!execution_) {
    return true;
  }

  // Update existing and add new functionObjects
  const entry* entryPtr = parentDict_.lookupEntryPtr
                          (
                              "functions",
                              false,
                              false
                          );

  if (entryPtr) {
    PtrList<functionObjectTest> newPtrs;
    List<SHA1Digest> newDigs;
    HashTable<label> newIndices;

    label nFunc = 0;

    if (!entryPtr->isDict()) {
      FatalIOErrorInFunction(parentDict_)
          << "'functions' entry is not a dictionary"
          << exit(FatalIOError);
    }

    const dictionary& functionsDict = entryPtr->dict();

    const_cast<Time&>(time_).libs().open
    (
        functionsDict,
        "libs",
        functionObjectTest::dictionaryConstructorTablePtr_
    );

    newPtrs.setSize(functionsDict.size());
    newDigs.setSize(functionsDict.size());

    forAllConstIter(dictionary, functionsDict, iter) {
      const word& key = iter().keyword();

      if (!iter().isDict()) {
        if (key != "libs") {
          IOWarningInFunction(parentDict_)
              << "Entry " << key << " is not a dictionary" << endl;
        }

        continue;
      }

      const dictionary& dict = iter().dict();
      bool enabled = dict.lookupOrDefault("enabled", true);

      newDigs[nFunc] = dict.digest();

      label oldIndex;
      functionObjectTest* objPtr = remove(key, oldIndex);

      if (objPtr) {
        if (enabled) {
          // Dictionary changed for an existing functionObject
          if (newDigs[nFunc] != digests_[oldIndex]) {
            ok = objPtr->read(dict) && ok;
          }
        } else {
          // Delete the disabled functionObject
          delete objPtr;
          objPtr = NULL;
          continue;
        }
      } else if (enabled) {
        autoPtr<functionObjectTest> foPtr;

        FatalError.throwExceptions();
        FatalIOError.throwExceptions();
        try {
          if
          (
              dict.found("writeControl")
              || dict.found("outputControl")
          ) {
            /*
                          foPtr.set
                          (
                              new functionObjects::timeControl(key, time_, dict)
                          );
            */
          } else {
            foPtr = functionObjectTest::New(key, time_, dict);
          }
        } catch (Foam::IOerror& ioErr) {
          Info<< ioErr << nl << endl;
          ::exit(1);
        } catch (Foam::error& err) {
          WarningInFunction
              << "Caught FatalError " << err << nl << endl;
        }
        FatalError.dontThrowExceptions();
        FatalIOError.dontThrowExceptions();

        if (foPtr.valid()) {
          objPtr = foPtr.ptr();
        } else {
          ok = false;
        }
      }

      // Insert active functionObjects into the list
      if (objPtr) {
        newPtrs.set(nFunc, objPtr);
        newIndices.insert(key, nFunc);
        nFunc++;
      }
    }

    newPtrs.setSize(nFunc);
    newDigs.setSize(nFunc);

    // Updating the PtrList of functionObjects deletes any
    // existing unused functionObjects
    PtrList<functionObjectTest>::transfer(newPtrs);
    digests_.transfer(newDigs);
    indices_.transfer(newIndices);
  } else {
    PtrList<functionObjectTest>::clear();
    digests_.clear();
    indices_.clear();
  }

  return ok;
}


void Foam::functionObjectListTest::updateMesh(const mapPolyMesh& mpm)
{
  if (execution_) {
    forAll(*this, objectI) {
      operator[](objectI).updateMesh(mpm);
    }
  }
}


void Foam::functionObjectListTest::movePoints(const polyMesh& mesh)
{
  if (execution_) {
    forAll(*this, objectI) {
      operator[](objectI).movePoints(mesh);
    }
  }
}


// ************************************************************************* //

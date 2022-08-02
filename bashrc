export PATO_VERSION=PATO-dev

# Figure out the complete path to the sourced script
if [ -n "$BASH_SOURCE" ]; then
  this_script=$BASH_SOURCE
elif [ -n "$ZSH_VERSION" ]; then
  setopt function_argzero
  this_script=$0
else
  echo 1>&2 "Unsupported shell. Please use bash or zsh."
  exit 2
fi

# If PATO_DIR is undefined, guess it from script location
if [ -z $PATO_DIR ]; then
   export PATO_DIR=`dirname $this_script`
fi

# Configure PATO
export PATH=$PATO_DIR/install/bin:$PATH
export LIB_PATO=$PATO_DIR/src/applications/libraries
export PATO_UNIT_TESTING=$PATO_DIR/src/applications/utilities/tests
export PATO_TUTORIALS=$PATO_DIR/tutorials
if [ "$(uname)" == "Darwin" ]; then
   export DYLD_LIBRARY_PATH=$PATO_DIR/install/lib:$DYLD_LIBRARY_PATH
else
   export LD_LIBRARY_PATH=$PATO_DIR/install/lib:$LD_LIBRARY_PATH
fi

# Configure M++ variables
export MPP_DIRECTORY=$PATO_DIR/src/thirdParty/mutation++
export MPP_DATA_DIRECTORY=$PATO_DIR/data/ThermoTransportChemistry/mutation++
export PATH=$MPP_DIRECTORY/install/bin:$PATH
if [ "$(uname)" == "Darwin" ]; then
   export DYLD_LIBRARY_PATH=$MPP_DIRECTORY/install/lib:$DYLD_LIBRARY_PATH
else
   export LD_LIBRARY_PATH=$MPP_DIRECTORY/install/lib:$LD_LIBRARY_PATH
fi

# Configure Python Plot directory 
export PYTHONPATH=$PATO_DIR/src/thirdParty/patoPlot/:$PYTHONPATH 
export PYTHONDONTWRITEBYTECODE=1 # do not write __pycache__ folder

# Documentation
BUILD_DOCUMENTATION=no

# Configure useful aliases
alias pato='cd $PATO_DIR'
alias solo='cd $PATO_DIR/src/applications/solvers'
alias utio='cd $PATO_DIR/src/applications/utilities'
alias libo='cd $PATO_DIR/src/applications/libraries'
alias tuto='cd $PATO_DIR/tutorials'
alias 1D='cd $PATO_DIR/tutorials/1D'
alias 1='cd $PATO_DIR/tutorials/1D/AblationTestCase_1.0'
alias 2D='cd $PATO_DIR/tutorials/2D'
alias 3D='cd $PATO_DIR/tutorials/3D'
alias muto='cd $MPP_DIRECTORY'

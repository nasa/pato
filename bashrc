export PATO_VERSION=PATO-dev

# Figure out the complete path to the sourced script
if [ -n "$BASH_SOURCE" ]; then
  this_script=$BASH_SOURCE
elif [ -n "$ZSH_VERSION" ]; then
  setopt function_argzero
  this_script=$0
else
  echo 1>&2 "Unsupported shell. Please use bash or zsh."
  return
fi

# If PATO_DIR is undefined, guess it from script location
if [ -z $PATO_DIR ]; then
   export PATO_DIR=`dirname $this_script`
fi

# Configure PATO
export PATH=$PATO_DIR/install/bin:$PATH
export PATH=$PATO_DIR/src/applications/utilities/tutoInDevelopment:$PATH
export LIB_PATO=$PATO_DIR/src/applications/libraries
export PATO_UNIT_TESTING=$PATO_DIR/src/applications/utilities/tests
export PATO_TUTORIALS=$PATO_DIR/tutorials
export UNAME=$(uname)
if [ "$(uname)" = "Darwin" ]; then
   export DYLD_LIBRARY_PATH=$PATO_DIR/install/lib:$DYLD_LIBRARY_PATH
else
   export LD_LIBRARY_PATH=$PATO_DIR/install/lib:$LD_LIBRARY_PATH
fi

# sed
if [ "$(uname)" = "Darwin" ]; then
   if [ ! `command -v gsed` ];  then
      echo 1>&2 Error: gsed command not found.
      return
   fi
   sed_cmd=gsed
else
   if [ ! `command -v sed` ];  then
      echo 1>&2 Error: sed command not found.
      return
   fi
   sed_cmd=sed
fi

# Remove FOAM_EXTEND_SRC from CPATH
REMOVE_PATH=$FOAM_EXTEND_SRC/finiteVolume/lnInclude
OLDPATH="$CPATH"; NEWPATH=""; colon=""
while [ "${OLDPATH#*:}" != "$OLDPATH" ]
do  entry="${OLDPATH%%:*}"; search=":${OLDPATH#*:}:"
    [ "${search#*:$entry:}" = "$search" ] && NEWPATH="$NEWPATH$colon$entry" && colon=:
    OLDPATH="${OLDPATH#*:}"
done
NEWPATH="$NEWPATH:$OLDPATH"
CPATH="$(echo "$NEWPATH" |$sed_cmd -e "s#\(^\|:\)$(echo "$REMOVE_PATH" |$sed_cmd -e 's/[^^]/[&]/g' -e 's/\^/\\^/g')\(:\|/\{0,1\}$\)#\1\2#" -e 's#:\+#:#g' -e 's#^:\|:$##g')"

# Configure foam-extend
if [[ ! -z "${FOAM_EXTEND_SRC}" ]] && [[ ! -z "${FOAM_EXTEND_LIB}" ]]; then
   if [ "$(uname)" = "Darwin" ]; then
      export DYLD_LIBRARY_PATH=$FOAM_EXTEND_LIB:$DYLD_LIBRARY_PATH
   else
      export LD_LIBRARY_PATH=$FOAM_EXTEND_LIB:$LD_LIBRARY_PATH
   fi
   CPATH=$FOAM_EXTEND_SRC/finiteVolume/lnInclude:$CPATH
fi
export CPATH

# PATO debug
if [ -z PATO_DEBUG ]; then
    export PATO_DEBUG="NO"
fi
debug_line="GFLAGS += -DFULLDEBUG -g -O0"
options_file="$PATO_DIR/src/applications/libraries/libPATOx/Make/options"
if [ "$PATO_DEBUG" = "YES" ]; then
   if ! grep -qF "$debug_line" $options_file; then
      $sed_cmd -i "\$a$debug_line" $options_file
   fi
else
   if grep -qF "$debug_line" $options_file; then
      $sed_cmd -i "s/$debug_line//" $options_file
      $sed_cmd -i "\$d" $options_file
   fi
fi

# Configure M++ variables
export MPP_DIRECTORY=$PATO_DIR/src/thirdParty/mutation++
export MPP_DATA_DIRECTORY=$PATO_DIR/data/ThermoTransportChemistry/mutation++
export PATH=$MPP_DIRECTORY/install/bin:$PATH
if [ "$(uname)" = "Darwin" ]; then
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
alias bco='cd $PATO_DIR/src/applications/libraries/libPATOx/MaterialModel/BoundaryConditions'
alias tuto='cd $PATO_DIR/tutorials'
alias 1D='cd $PATO_DIR/tutorials/1D'
alias 1='cd $PATO_DIR/tutorials/1D/AblationTestCase_1.0'
alias 2D='cd $PATO_DIR/tutorials/2D'
alias 3D='cd $PATO_DIR/tutorials/3D'
alias muto='cd $MPP_DIRECTORY'

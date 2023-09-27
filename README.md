<div align="center">
<img src=
"http://pato.ac/wp-content/uploads/2022/06/2022-06-07_PATO_logo_transparent.png"
width="35%">
</div>
<br />

[![Release](https://img.shields.io/badge/Release-github/pato-blue)](https://github.com/nasa/pato)
[![Anaconda-Server Badge](https://anaconda.org/pato.devel/pato/badges/version.svg)](https://anaconda.org/pato.devel/pato)
[![Documentation](https://img.shields.io/badge/Documentation-3.1-brightgreen)](https://github.com/nasa/pato/blob/pato-3.1/documentation/user_guide/PATO_v3.1_User_Guide.pdf)
[![Information](https://img.shields.io/badge/Information-pato.ac-brightgreen)](https://pato.ac)



# Porous material Analysis Toolbox based on OpenFoam (PATO) 3.1

## Notices

Copyright © 2023 United States Government as represented by the Administrator of the National Aeronautics and Space Administration.  All Rights Reserved.

## Disclaimers

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT. 

## Legal Information

* Distributed under NASA Open Source Agreement Version 1.3.
* Please visit https://pato.ac for more information.

## Installation Guideline Summary

### Binaries using conda (UNIX)

To install PATO, a conda distribution must be installed on your machine. 
To test whether conda is installed, run "conda --version" from a terminal to see if the command is recognized. 
If not, conda can be installed by following the instructions 
[here](https://docs.anaconda.com/anaconda/install/index.html).

Once the conda command is working, all the PATO components can be installed by executing 
the following commands in a terminal (note that the "solving environment" step can sometimes take up to 15 mins):

```bash
conda config --add channels conda-forge 
conda config --add channels pato.devel 
conda config --set channel_priority strict 
conda create -y --name pato -c conda-forge -c pato.devel pato
```

On **UNIX** (i.e. Mac or Linux), the openfoam_for_pato and foam-extend_for_pato packages are installed. 

PATO relies on a conda environment in order to manage its software dependencies and environment variables.
It is therefore important to always activate the environment before using any of PATO's functionalities. 
Once the installation is complete, the PATO software can be tested by running (note that the tests can take up to 2 hours):

```bash
conda activate pato
runtests
```

To uninstall PATO and all the installed dependencies, execute the following command to delete the PATO environment:

```bash
conda remove -y --name pato --all
```

For PATO developers that have access to the GitLab repository, a dev version of the PATO conda package is available.
The dev version will install OpenFOAM/foam-extend and clone the latest version of PATO from the GitLab repository.

First, ssh needs to be configured for Gitlab.
* Create the SSH key in .ssh: `ssh-keygen`
* Add the public key to https://gitlab.com/-/profile/keys
* Modify $HOME/.ssh/config as follows

```bash
Host gitlab.com
  Hostname gitlab.com
  PreferredAuthentications publickey
  IdentityFile ~/.ssh/id_rsa
```

To test the ssh connection, try to clone the GitLab repository.
```bash
git clone git@gitlab.com:PATO/PATO-dev.git
```

To install the dev version of the PATO conda package (note that the clone and build can take up to 1 hour):

```bash
conda create -y --name pato-dev -c conda-forge -c pato.devel pato=dev
```

Once the installation is complete, the PATO software can be tested by running (note that the tests can take up to 2 hours):

```bash
conda activate pato-dev
runtests
```

To uninstall PATO dev and all the installed dependencies, execute the following command to delete the PATO environment:

```bash
conda remove -y --name pato-dev --all
```

### Build from Source (UNIX)

0.  Recommended system (for a first installation):
    * OS: Ubuntu Linux 16.04
    * Shell: bash
    * Text editor: gedit
    * Plotting software: python v2.7.5 or newer
    * Optional installation: numpy, matplotlib

1.  Install OpenFOAM 7 to your home directory by compiling it from sources - see for example www.openfoam.org.

2.  Copy your PATO version in the directory of your choice, for example, copy it to `$HOME/PATO` or directly clone it from the git repository.

3.  Source the PATO `bashrc` from your `$HOME/.bashrc` file (depending on where you installed PATO), for example, add this line to your `$HOME/.bashrc`:
    ```bash
    export PATO_DIR=$HOME/PATO/PATO-dev
    source $PATO_DIR/bashrc
    ```

4.  To compile PATO, open a terminal and execute the following commands:
    ```bash
    source $HOME/.bashrc
    pato
    ./AllwcleanAllclean
    ./Allwmake
    ```

5.  If you can use Gedit, let's optimise its layout for PATO files. In the menu, go to edit/preferences:
    * select font : courrier 10 pitch /11
    * uncheck "enable text wrapping"
    * set the indentation size to 4 spaces (do not use tabs)

6.  Your system is ready to run PATO tests and tutorials.
    ```bash
    runtests
    ```

## NASA Registration

NASA requires users to register for statistical purposes on a good will basis.
    Please send us an email at admin@pato.ac with your Name and Affiliation.
    We are also interested in knowing a little about your projects if this is something you can share.

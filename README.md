# Porous material Analysis Toolbox based on OpenFoam (PATO) 2.3

## Notices

Copyright © 2010 United States Government as represented by the Administrator of the National Aeronautics and Space Administration.  All Rights Reserved.

Disclaimers

* No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
* Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.

## Legal Information

* Distributed under NASA Open Source Agreement Version 1.3.
* Please visit https://pato.ac for more information.

## Installation Guideline Summary For Personal Computers

0.  Recommended system (for a first installation):
    * OS: Ubuntu Linux 16.04
    * Shell: bash
    * Text editor: gedit
    * Plotting software: python v2.7.5 or newer
    * Optional installation: numpy, matplotlib

1.  Install OpenFOAM 7 to your home directory by compiling it from sources - see for example www.openfoam.org.

2.  Copy your PATO version in the directory of your choice, for example, copy it to `$HOME/PATO` or directly clone it from the git repository.

3.  Source the PATO `bashrc` from your `bashrc` file (depending on where you installed PATO), for example, add this line to your `.bashrc`:
    ```bash
    source $HOME/PATO/PATO-dev/bashrc
    ```

4.  To compile PATO, open a terminal and execute the following commands:
    ```bash
    source .bashrc
    pato
    ./AllwcleanAllclean
    ./Allwmake
    ```

5.  If you can use Gedit, let's optimise its layout for PATO files. In the menu, go to edit/preferences:
    * select font : courrier 10 pitch /11
    * uncheck "enable text wrapping"
    * set the indentation size to 4 spaces (do not use tabs)

6.  Your system is ready to run PATO tutorials.

7.  NASA requires users to register for statistical purposes on a good will basis.
    Please send us an email at admin@pato.ac with your Name and Affiliation.
    We are also interested in knowing a little about your projects if this is something you can share.

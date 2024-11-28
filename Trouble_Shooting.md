# Installation
   For errors in installation phase, we recommend the following steps to try. Note that the error messages could be different for Mac Users regarding M1/M2/M3 chips.

- In terminal:
1. Ensure Xcode Command Line Tools are installed
  ```
  xcode-select --install
  ```
2. Ensure Homebrew is installed
   ```
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```
3. Install 'gfortran' and other dependencies
  ```
  brew install gfortran
  brew install llvm
  brew install gcc
  ```

   
- Check R configuration
1. create a file named Makevars in your .R directory, put the following lines in your .R/Makevars file
   ```
   FC = $path_gfortran
   F77 = $path_gfortran
   F
   ```
- In R:
1. Install devtools
   ```
   install.packages("devtools")
   library(devtools)
   ```
2. Install causl
   ```
   install_github("rje42/causl")
   ```
   

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

   

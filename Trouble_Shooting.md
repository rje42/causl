# Installation
For errors during the installation phase, we recommend trying the following steps. Please note that error messages may vary for Mac users, particularly those using M1, M2, or M3 chips.
## In terminal:
1. Ensure Xcode Command Line Tools are installed
   ```
   xcode-select --install
   ```
2. Ensure Homebrew is installed
   ```
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```
3. Install 'gfortran' and other dependencies, including gcc
   ```
   brew install gfortran
   brew install llvm
   brew install gcc
   ```
4. Add `gcc` to your path
   - Option 1: Create aliases
      Add the following lines to your shell configuration file (~/.zshrc for Zsh or ~/.bashrc for Bash):
      ```
      alias gcc='/opt/homebrew/bin/gcc-14'
      alias g++='/opt/homebrew/bin/g++-14'
      ```
      (Here we assume your gcc is version 14. you should run
      ```
      gcc-14 --version
      g++-14 --version
      ```
      to verify it.)
      Then apply the changes:
      ```
      source ~/.zshrc
      ```
   - Option 2: Create Symlinks
      ```
      sudo ln -s /opt/homebrew/bin/gcc-14 /usr/local/bin/gcc
      sudo ln -s /opt/homebrew/bin/g++-14 /usr/local/bin/g++
      ```

## Check R configuration:
1. create a file named Makevars in your .R directory, put the following lines in your .R/Makevars file
   ```
   FC = (path_gfortran)
   F77 = (path_gfortran)
   FLIBS=-L(path_gcc) -lgfortran -lquadmath -lm
   ```
   where you should replace the (path_gfortran) with path returned by running the following command in your terminal:
   ```
   brew info gfortran
   ```
   you should be able to get something similar to:
   ```
   FC = /opt/homebrew/bin/gfortran
   F77 = /opt/homebrew/bin/gfortran
   ```
   and (path_gcc) with output by running this command in your terminal:
   ```
   brew info gcc
   ```
   and you will get a line similar to:
   ```
   -L/opt/homebrew/opt/gcc/lib/gcc/14 -lgfortran -lquadmath -lm
   ```
## In R:
1. Install devtools
   ```
   install.packages("devtools")
   library(devtools)
   ```
2. Install causl
   ```
   install_github("rje42/causl")
   ```

If you encounter any issues, please send us the error messages!

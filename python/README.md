# Python Communication

We provide examples of customizing data generation in `causl` and calling the function in Python.

## Requirements

Ensure you have installed the Python package `rpy2`. If it’s not already installed, run:

```bash
pip install rpy2
```

## Necessary files
You’ll need two files:
1. An `.R` file that defines and wraps the data-generation process.
- Make sure the function returns the simulated data. 
- Expose any parameters you’d like to customize as function arguments so you can specify them when calling from Python.
- See `data_causl.R' in this folder for an example. 

2. A `.py` file that defines the Python function for passing those parameters.
- We assume all data formats in Python use NumPy.
- Convert incoming NumPy arrays to the corresponding R types before calling your R function.
- See `data_causl.py' in this folder for an example.


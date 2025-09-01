import numpy as np
import pandas as pd

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri,numpy2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr


import io
import contextlib
import warnings
import os
import sys

here = os.path.dirname(__file__)
r_script = os.path.join(here, 'data_causl.R')

warnings.filterwarnings("ignore", message="R is not initialized by the main thread")

# Suppress R output
# If you would like to see the R output, please comment this out.
@contextlib.contextmanager
def suppress_r_output():
    r_output = io.StringIO()
    with contextlib.redirect_stdout(r_output), contextlib.redirect_stderr(r_output):
        yield


def gen_causl_example1(n=10000):
    pandas2ri.activate()
    # Source the ./data.r script for data.causl dgp function
    with suppress_r_output():
        # get the folder containing utils.py
        robjects.r['source'](r_script)
        gen_causl = robjects.globalenv['gen_causl_example1']
        r_dataframe = gen_causl(n=n)
    # Use the localconverter context manager to convert the R dataframe to a Pandas DataFrame
    with localconverter(robjects.default_converter + pandas2ri.converter):
        df = robjects.conversion.rpy2py(r_dataframe)
    return df

def gen_causl_example2(n=10000, nI = 3, nX= 1, nO = 1, nS = 1, ate = 2, beta_cov = 0, strength_instr = 3, strength_conf = 1, strength_outcome = 1, binary_intervention=True):
    pandas2ri.activate()
    # Source the ./data.r script for data.causl dgp function
    with suppress_r_output():
        # get the folder containing utils.py
        robjects.r['source'](r_script)
        gen_causl = robjects.globalenv['gen_causl_example2']
        r_dataframe = gen_causl(n=n, nI=nI, nX=nX, nO=nO, nS=nS, ate=ate, beta_cov=beta_cov, strength_instr=strength_instr, strength_conf=strength_conf, strength_outcome=strength_outcome, binary_intervention=binary_intervention)
    # Use the localconverter context manager to convert the R dataframe to a Pandas DataFrame
    with localconverter(robjects.default_converter + pandas2ri.converter):
        df = robjects.conversion.rpy2py(r_dataframe)
    return df

## Test example
gen_causl_example1(n=100)
''' External Python packages
* `pickle`
* `requests`
* itertools.product (`product`)
* `re`
* Bio.Seq.Seq (`BioSeq`)
* pandas (`pd`)
* numpy (`np`)
* matplotlib.pyplot (`plt`)
* matplotlib.patches (`patches`)
* `faicons`
* shiny (`App`, `Inputs`, `Outputs`, `Session`, `reactive`, `render`, `ui`)

Other:
* `ROOT`: root directory filepath
'''

__all__ = ["pickle", "requests", "product", "re", "BioSeq", "pd", "np", "plt", "patches", "faicons",
           "App", "Inputs", "Outputs", "Session", "reactive", "render", "ui", "ROOT"]

import pickle
import requests
from itertools import product
import re
from Bio.Seq import Seq as BioSeq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import faicons
from shiny import App, Inputs, Outputs, Session, reactive, render, ui

from pathlib import Path
ROOT = Path(__file__).parent.parent

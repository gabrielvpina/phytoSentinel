import os, glob
from Bio import SeqIO
import pandas as pd
import collections
import pyhmmer

from .data import pfam_metadata
from .data import pfam_complete_desc


ResultHMM = collections.namedtuple("Result", ["query", "subjID", "bitscore", "start", "end"])


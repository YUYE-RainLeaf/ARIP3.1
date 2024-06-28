from pathlib import Path
from typing import *

import numpy as np
from numpy import ndarray
from numpy.typing import NDArray
from pandas import DataFrame, Series

# Numpy dtype
DTYPE = np.float64

FileLike = object()
AtomModel = List[str]

Range  = Tuple[float, float]
Point  = NDArray[DTYPE]    # [D=3]
Points = NDArray[DTYPE]    # [N, D=3]

a_pair = Tuple[str, str]
Atom = Tuple[str, float, Points]        # Atom name, radius, coordinate
AtomGen = Generator[Atom, None, None]
PointAtomGen = Generator[Tuple[Point, Atom], None, None]
Contact = Tuple[str, float,  str]       # radius, surface, volume

Angles   = Dict[str, Tuple[float, float]]
SASA     = Dict[str, float]
Surfaces = Dict[a_pair, List[Union[float, str]]]
Volumes  = Dict[a_pair, float]

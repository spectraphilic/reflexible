# Import the public functions here
from .legacy_structures import Header, Structure, BinaryFile, FDC
from .flexpart_read import (
    read_header, read_trajectories, read_command, read_releases)
from .grid_read import read_grid, get_slabs, fill_grids
from .helpers import get_fpdirs

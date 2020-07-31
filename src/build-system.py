import mbuild as mb
import atools
from mbuild.lib.atoms import H
from atools.lib.chains import Alkylsilane
from atools.recipes import DualSurface, SilicaInterface, SurfaceMonolayer
from atools.lib.chains.alkylsilane_internal import Alkylsilane as AlkylsilaneInternal

def init_system(job):
    c_a = 17
    c_b = 17
    c_c = 17
    c_d = 17
    seed = 12345
    pattern_type = 'random'
    terminal_groups

# -*- coding: utf-8 -*-
import textwrap
import glob
import json
import logging
import os
import pathlib

import flow
import numpy as np
import scipy
import signac
from flow import FlowProject, environments
from foyer import Forcefield
from mbuild.formats.lammpsdata import write_lammpsdata
from mbuild.lib.atoms import H

from atools.fileio import write_monolayer_ndx, read_ndx
from atools.lib.chains import Alkylsilane
from atools.recipes import DualSurface, SilicaInterface, SurfaceMonolayer
from util.index_groups import generate_index_groups
from atools.lib.chains.alkylsilane_internal import Alkylsilane as AlkylsilaneInternal


def _setup_logger(logger_name, log_file, level=logging.INFO):
    mylog = logging.getLogger(logger_name)
    formatter = logging.Formatter("%(asctime)s : %(message)s")
    fileHandler = logging.FileHandler(log_file, mode="a")
    fileHandler.setFormatter(formatter)

    mylog.setLevel(level)
    mylog.addHandler(fileHandler)
    return mylog


def _mdrun_str(op_name):
    """return string formatted based on name of operation"""
    msg = (
        "mdrun_mpi_sp -v -deffnm {op} -s {op}.tpr -cpi {op}.cpt".format(op=op_name)
    )
    return msg


def __trjcat_naive_shear(shear_force, file_list, extension):
    msg = "gmx_sp trjcat -f {} -o shear_{}nN_combined.{}".format(
        " ".join(file_list), shear_force, extension
    )
    return msg


def _wraprun_add_tasks_helper(n_cpus, working_dir, command):
    """For serial jobs, build up a large string of executables, sep by :"""
    msg = "cd {}; {} ".format(working_dir, command)
    return msg

"""Configuration of the project enviroment.

The environments defined in this module can be auto-detected.
This helps to define environment specific behaviour in heterogenous
environments.
"""
class CoriEnvironment(flow.environment.DefaultSlurmEnvironment):

    hostname_pattern = 'cori'
    template='cori.sh'
    cores_per_node=64
    
    @classmethod
    def add_args(cls, parser):
        super(flow.environment.SlurmEnvironment, cls).add_args(parser)
        parser.add_argument(
                '-w', "--walltime", type=float, help="Walltime"
                )
        parser.add_argument(
                "--name", type=str, help="Operation name"
                )
        parser.add_argument(
                "--queue", type=str, help="Queue type"
                )
        parser.add_argument(
                "--mail", type=str, help="Mail option", default="NONE"
                )

class Project(FlowProject):
    pass


@Project.label
def is_initialized(job):
    result = bool(
        job.isfile("init.gro")
        and job.isfile("init.top")
        and job.isfile("init.lammps")
        and job.isfile("init.ndx")
    )
    return result


@Project.label
def overlaps_fixed(job):
    return job.isfile("minimize.xtc")


@Project.label
def converted_lmp(job):
    return job.isfile("minimize.gro")


@Project.label
def made_ndx(job):
    return job.isfile("init2.ndx")


@Project.label
def is_em_grompp(job):
    return job.isfile("em.tpr")


@Project.label
def em_completed(job):
    return job.isfile("em.gro")


@Project.label
def is_nvt_grompp(job):
    return job.isfile("nvt.tpr")


@Project.label
def nvt_completed(job):
    return job.isfile("nvt.gro")


@Project.label
def is_compress_grompp(job):
    return job.isfile("compress.tpr")


@Project.label
def compress_completed(job):
    return job.isfile("compress.gro")


@Project.label
def is_shear_5nN_grompp(job):
    return job.isfile("shear_5nN.tpr")


@Project.label
def shear_5nN_completed(job):
    return bool(glob.glob(job.fn("shear_5nN*.gro")))


@Project.label
def is_shear_15nN_grompp(job):
    return job.isfile("shear_15nN.tpr")


@Project.label
def shear_15nN_completed(job):
    return bool(glob.glob(job.fn("shear_15nN*.gro")))


@Project.label
def is_shear_25nN_grompp(job):
    return job.isfile("shear_25nN.tpr")


@Project.label
def shear_25nN_completed(job):
    return bool(glob.glob(job.fn("shear_25nN*.gro")))


@Project.label
def friction_calculated(job):
    is_calculated = list()
    for shear in [5, 15, 25]:
        is_calculated.append(job.isfile("friction_{}nN.txt".format(shear)))
    return all(is_calculated)


@Project.label
def cof_calculated(job):
    "if the job document has a listing for COF and intercept"
    is_calculated = (job.document.get("COF"), job.document.get("intercept"))
    return all(is_calculated)


@Project.operation
@Project.post.isfile("init.top")
@Project.post.isfile("init.gro")
@Project.post.isfile("init.lammps")
@Project.post.isfile("init.ndx")
@flow.with_job
def initialize_system(job):
    """ Generate the monolayer surfaces, parametrize, save LAMMPS, GRO, TOP.
    """
    """
    ---------------------------
    Read statepoint information
    ---------------------------
    """
    chainlength_a = job.statepoint()["chainlength_a"]
    chainlength_b = job.statepoint()["chainlength_b"]
    chainlength_c = job.statepoint()["chainlength_c"]
    chainlength_d = job.statepoint()["chainlength_d"]
    seed = job.statepoint()["seed"]
    pattern_type = job.statepoint()["pattern_type"]
    terminal_groups = job.statepoint()["terminal_groups"]
    a_fraction = job.statepoint()["fraction_a"]
    c_fraction = job.statepoint()["fraction_c"]
    num_chains = job.statepoint()["n"]
    backbone_1 = job.statepoint()["backbone_1"]
    backbone_2 = job.statepoint()["backbone_2"]
    locations = job.statepoint()["locations"]
    """
    -----------------------------------
    Generate amorphous silica interface
    -----------------------------------
    """
    surface_a = SilicaInterface(thickness=1.2, seed=seed)
    surface_b = SilicaInterface(thickness=1.2, seed=seed)

    """
    ------------------------------------------------------
    Generate prototype of functionalized alkylsilane chain
    ------------------------------------------------------
    """
    chain_prototype_a = AlkylsilaneInternal(
        chain_length=chainlength_a, internal_group=backbone_1, locations=locations, terminal_group=terminal_groups[0]
    )
    chain_prototype_b = AlkylsilaneInternal(
        chain_length=chainlength_b, internal_group=backbone_2, locations=locations, terminal_group=terminal_groups[1]
    )
    chain_prototype_c = AlkylsilaneInternal(
        chain_length=chainlength_c, internal_group=backbone_1, locations=locations, terminal_group=terminal_groups[2]
    )
    chain_prototype_d = AlkylsilaneInternal(
	chain_length=chainlength_d, internal_group=backbone_2, locations=locations, terminal_group=terminal_groups[3]
    )
    """
    ----------------------------------------------------------
    Create monolayer on surface, backfilled with hydrogen caps
    ----------------------------------------------------------
    """
    # bottom monolayer is backfilled with the other terminal group
    # num_chains = num_chains * a_fraction
    monolayer_a = SurfaceMonolayer(
        surface=surface_a,
        chains=[chain_prototype_a, chain_prototype_b],
        n_chains=num_chains,
        seed=seed,
        backfill=H(),
        rotate=False,
        fractions=[a_fraction, 1.0 - a_fraction],
    )
    monolayer_a.name = "Bottom"
    monolayer_b = SurfaceMonolayer(
        surface=surface_b,
        chains=[chain_prototype_c, chain_prototype_d],
        n_chains=num_chains,
        seed=seed,
        backfill=H(),
        rotate=False,
	fractions=[c_fraction, 1.0 - c_fraction],
    )
    monolayer_b.name = "Top"

    """
    ------------------------------------------
    Duplicate to yield two opposing monolayers
    ------------------------------------------
    """
    dual_monolayer = DualSurface(
        bottom=monolayer_a, top=monolayer_b, separation=2.0
    )

    """
    --------------------------------------------------------
    Make sure box is elongated in z to be pseudo-2D periodic
    --------------------------------------------------------
    """
    box = dual_monolayer.boundingbox
    dual_monolayer.periodicity += np.array([0, 0, 5.0 * box.lengths[2]])

    """
    -------------------------------------------------------------------
    - Save to .GRO, .TOP, and .LAMMPS formats
    - Atom-type the system using Foyer, with parameters from the OPLS
    force field obtained from GROMACS. Parameters are located in a
    Foyer XML file in the `atools` git repo, with references provided
    as well as notes where parameters have been added or altered to
    reflect the literature.
    -------------------------------------------------------------------
    """
    # path for project root dir
    proj = signac.get_project()
    forcefield_filepath = pathlib.Path(
        proj.root_directory() + "/src/util/forcefield/oplsaa.xml"
    )
    # change into job directoryA
    _switch_dir(job)
    logging.info("at dir: {}".format(job.ws))
    for p in dual_monolayer.particles():
        if p.name == "OS":
            p.name = "O"
    dual_monolayer.save("init.gro", residues=["Top", "Bottom"], overwrite=True)

    
            
            
    if not (
        job.isfile("init.top")
        and job.isfile("init.lammps")
        and job.isfile("init.gro")
    ):

        structure = dual_monolayer.to_parmed(
            box=None, residues=["Top", "Bottom"]
        )
        ff = Forcefield(forcefield_files=forcefield_filepath.as_posix())
        structure = ff.apply(structure)
        structure.combining_rule = "geometric"

        structure.save("init.top", overwrite=True)
        write_lammpsdata(filename="init.lammps", structure=structure, detect_forcefield_style=False)

        """
        --------------------------------------
        Specify index groups and write to file
        --------------------------------------
        """
    index_groups = generate_index_groups(
        system=dual_monolayer,
        terminal_groups=terminal_groups,
        freeze_thickness=0.5,
    )
    write_monolayer_ndx(rigid_groups=index_groups, filename="init.ndx")


def _switch_dir(job):
    p = pathlib.Path(job.workspace())
    os.chdir(str(p.absolute()))


@Project.operation
@Project.pre.after(initialize_system)
@Project.post.isfile("minimize.xtc")
@flow.cmd
def fix_overlaps(job):
    cmds = pathlib.Path(
        signac.get_project().root_directory() + "/src/util/mdp_files"
    )
    return "cd {}; srun -n 1 lmp_cori -in {}/in.minimize -log {}/minimize.log".format(
        job.workspace(), str(cmds.absolute()), job.workspace()
    )


@Project.operation
@Project.pre.after(fix_overlaps)
@Project.post.isfile("minimize.gro")
@flow.cmd
def lmp_to_gmx(job):
    msg = "cd {}; echo 0 | mpirun -np 1 gmx_sp trjconv -s init.gro -f minimize.xtc -o minimize.gro  -b 1.0 -e 1.0".format(job.workspace())
    return msg


@Project.operation
@Project.pre.after(lmp_to_gmx)
@Project.post.isfile("init2.ndx")
@flow.cmd
def make_ndx(job):
    '''Potentially integrate into another step'''
    msg = "cd {}; gmx_sp make_ndx -f minimize.gro -n init.ndx -o init2.ndx < ../../src/second-ndx.txt".format(job.workspace())
    return msg


@Project.operation
@Project.pre.after(make_ndx)
@Project.post.isfile("em.tpr")
@flow.cmd
def em_grompp(job):
    em_mdp_path = pathlib.Path(
        signac.get_project().root_directory() + "/src/util/mdp_files/em.mdp"
    )
    msg = "cd {}; gmx_sp grompp -f {} -c minimize.gro -p init.top -n init2.ndx -o em.tpr -maxwarn 1".format(job.workspace(), em_mdp_path)
    return msg


@Project.operation
@Project.pre.after(em_grompp)
@Project.post.isfile("em.gro")
@flow.cmd
def mdrun_em(job):
    msg = _mdrun_str("em")
    return "cd {}; srun -n 1 {} -ntomp 64".format(job.workspace(), msg)


@Project.operation
@Project.pre.after(mdrun_em)
@Project.post.isfile("nvt.tpr")
@flow.cmd
def nvt_equil_grompp(job):
    nvt_mdp_path = pathlib.Path(
        signac.get_project().root_directory() + "/src/util/mdp_files/nvt.mdp"
    )
    msg = "cd {}; gmx_sp grompp -f {} -c {} -p {} -n {} -o {} -maxwarn 1".format(
        job.workspace(),
        nvt_mdp_path,
        "em.gro",
        "init.top",
        "init2.ndx",
        "nvt.tpr",
    )
    return msg


@Project.operation
@Project.pre.after(nvt_equil_grompp)
@Project.post.isfile("nvt.gro")
@flow.cmd
def mdrun_nvt(job):
    msg = _mdrun_str("nvt")
    return "cd {}; srun -n 1 {} -ntomp 64".format(job.workspace(), msg)


@Project.operation
@Project.pre.after(mdrun_nvt)
@Project.post.isfile("compress.tpr")
@flow.cmd
def compress_grompp(job):
    compress_mdp_path = pathlib.Path(
        signac.get_project().root_directory()
        + "/src/util/mdp_files/test-compress.mdp"
    )
    msg = "cd {}; gmx_sp grompp -f {} -c {} -p {} -n {} -o {} -maxwarn 2".format(
        job.workspace(),
        compress_mdp_path,
        "nvt.gro",
        "init.top",
        "init2.ndx",
        "compress.tpr",
    )
    return msg


@Project.operation
@Project.pre.after(compress_grompp)
@Project.post.isfile("compress.gro")
@flow.cmd
def mdrun_compress(job):
    msg = _mdrun_str("compress")
    return "cd {}; srun -n 1 -v {} -ntomp 64".format(job.workspace(), msg)


@Project.operation
@Project.pre.after(mdrun_compress)
@Project.post.isfile("shear_5nN.tpr")
@flow.cmd
def shear_5nN_grompp(job):
    shear_5nN_mdp_path = pathlib.Path(
        signac.get_project().root_directory()
        + "/src/util/mdp_files/shear_5nN.mdp"
    )
    msg = "cd {}; mpirun -np 1 gmx_sp grompp -f {} -c {} -p {} -n {} -o {} -maxwarn 1".format(
    	job.workspace(),
        shear_5nN_mdp_path,
        "compress.gro",
        "init.top",
        "init2.ndx",
        "shear_5nN.tpr",
    )
    return msg


@Project.operation
@Project.pre.after(shear_5nN_grompp)
@Project.post(shear_5nN_completed)
@flow.cmd
def mdrun_shear_5nN(job):
    return "cd {}; srun -n 1 -v mdrun_mpi_sp -v -s shear_5nN.tpr -deffnm shear_5nN -cpi shear_5nN.cpt -cpo shear_5nN.cpt -noappend -ntomp 64".format(
        job.workspace()
    )


@Project.operation
@Project.pre.after(mdrun_compress)
@Project.post.isfile("shear_15nN.tpr")
@flow.cmd
def shear_15nN_grompp(job):
    shear_15nN_mdp_path = pathlib.Path(
        signac.get_project().root_directory()
        + "/src/util/mdp_files/shear_15nN.mdp"
    )
    msg = "cd {}; mpirun -np 1 gmx_sp grompp -f {} -c {} -p {} -n {} -o {} -maxwarn 1".format(
        job.workspace(),
        shear_15nN_mdp_path,
        "compress.gro",
        "init.top",
        "init2.ndx",
        "shear_15nN.tpr",
    )
    return msg


@Project.operation
@Project.pre.after(shear_15nN_grompp)
@Project.post(shear_15nN_completed)
@flow.cmd
def mdrun_shear_15nN(job):
    # msg = _mdrun_str("shear_15nN")
    return "cd {}; srun -n 1 mdrun_mpi_sp -s shear_15nN.tpr -deffnm shear_15nN -cpi shear_15nN.cpt -cpo shear_15nN.cpt -noappend -ntomp 64".format(
        job.workspace()
    )


@Project.operation
@Project.pre.after(mdrun_compress)
@Project.post.isfile("shear_25nN.tpr")
@flow.cmd
def shear_25nN_grompp(job):
    shear_25nN_mdp_path = pathlib.Path(
        signac.get_project().root_directory()
        + "/src/util/mdp_files/shear_25nN.mdp"
    )
    msg = "cd {}; mpirun -np 1 gmx_sp grompp -f {} -c {} -p {} -n {} -o {} -maxwarn 1".format(
        job.workspace(),
        shear_25nN_mdp_path,
        "compress.gro",
        "init.top",
        "init2.ndx",
        "shear_25nN.tpr",
    )
    return msg


@Project.operation
@Project.pre.after(shear_25nN_grompp)
@Project.post(shear_25nN_completed)
@flow.cmd
def mdrun_shear_25nN(job):
    # msg = _mdrun_str('shear_25nN')
    return "cd {}; srun -n 1 mdrun_mpi_sp -s shear_25nN.tpr -deffnm shear_25nN -cpi shear_25nN.cpt -cpo shear_25nN.cpt -noappend -ntomp 64".format(
        job.workspace()
    )


@Project.operation
@Project.pre.after(mdrun_shear_5nN)
@Project.post.isfile("shear_5nN_combined.xtc")
@flow.cmd
def shear_5nN_xtc_concat(job):
    file_list = glob.glob(job.workspace() + "/shear_5nN*.xtc")
    executable = __trjcat_naive_shear(5, file_list, "xtc")
    return _wraprun_add_tasks_helper(1, job.workspace(), executable)


@Project.operation
@Project.pre.after(mdrun_shear_5nN)
@Project.post.isfile("shear_5nN_combined.trr")
@flow.cmd
def shear_5nN_trr_concat(job):
    file_list = glob.glob(job.workspace() + "/shear_5nN*.trr")
    executable = __trjcat_naive_shear(5, file_list, "trr")
    return _wraprun_add_tasks_helper(1, job.workspace(), executable)


@Project.operation
@Project.pre.after(mdrun_shear_15nN)
@Project.post.isfile("shear_15nN_combined.xtc")
@flow.cmd
def shear_15nN_xtc_concat(job):
    file_list = glob.glob(job.workspace() + "/shear_15nN*.xtc")
    executable = __trjcat_naive_shear(15, file_list, "xtc")
    return _wraprun_add_tasks_helper(1, job.workspace(), executable)


@Project.operation
@Project.pre.after(mdrun_shear_15nN)
@Project.post.isfile("shear_15nN_combined.trr")
@flow.cmd
def shear_15nN_trr_concat(job):
    file_list = glob.glob(job.workspace() + "/shear_15nN*.trr")
    executable = __trjcat_naive_shear(15, file_list, "trr")
    return _wraprun_add_tasks_helper(1, job.workspace(), executable)


@Project.operation
@Project.pre.after(mdrun_shear_25nN)
@Project.post.isfile("shear_25nN_combined.xtc")
@flow.cmd
def shear_25nN_xtc_concat(job):
    file_list = glob.glob(job.workspace() + "/shear_25nN*.xtc")
    executable = __trjcat_naive_shear(25, file_list, "xtc")
    return _wraprun_add_tasks_helper(1, job.workspace(), executable)


@Project.operation
@Project.pre.after(mdrun_shear_25nN)
@Project.post.isfile("shear_25nN_combined.trr")
@flow.cmd
def shear_25nN_trr_concat(job):
    file_list = glob.glob(job.workspace() + "/shear_25nN*.trr")
    executable = __trjcat_naive_shear(25, file_list, "trr")
    return _wraprun_add_tasks_helper(1, job.workspace(), executable)


@Project.operation
@Project.pre.after(shear_5nN_trr_concat)
@Project.pre.after(shear_5nN_xtc_concat)
@Project.pre.after(shear_15nN_trr_concat)
@Project.pre.after(shear_15nN_xtc_concat)
@Project.pre.after(shear_25nN_trr_concat)
@Project.pre.after(shear_25nN_xtc_concat)
@Project.post(friction_calculated)
def calc_friction_system(job):
    from atools.structure_mixed import calc_friction

    for load in [5, 15, 25]:
        trr_file = "{}/shear_{}nN_combined.trr".format(job.workspace(), load)
        out_file = "{}/friction_{}nN.txt".format(job.workspace(), load)
        ndx_file = "{}/init2.ndx".format(job.workspace())
        calc_friction(
            trr_filename=trr_file,
            output_filename=out_file,
            ndx_filename=ndx_file,
        )


@Project.operation
@Project.pre.after(calc_friction_system)
@Project.post(cof_calculated)
def calc_cof(job):
    from scipy import stats

    loads = [5, 15, 25]
    friction_forces = []
    for load in loads:
        friction_data = np.loadtxt(
            "{}/friction_{}nN.txt".format(job.workspace(), load)
        )
        friction_data = friction_data[int(len(friction_data) * 0.6) :, 1]
        friction_forces.append(np.mean(friction_data))
    cof, intercept, r, p, stderr = stats.linregress(loads, friction_forces)
    job.document["COF"] = cof
    job.document["intercept"] = intercept


if __name__ == "__main__":
    Project().main()

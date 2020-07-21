{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash -l
#SBATCH --job-name={{ name }}
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time={{ walltime|format_timedelta }}
#SBATCH -C knl
#SBATCH --error=error.err
#SBATCH -q {{ queue }}
#SBATCH --mail-type={{ mail }}

module unload gromacs lammps
module load gromacs/2020.1.knl lammps/2018.12.12-knl openmpi
conda activate mosdef36
PYTHONIOENCODING=UTF-8


{% if operation.directives.omp_num_threads %}
export OMP_NUM_THREADS={{ operation.directives.omp_num_threads }}
{% endif %}
{{ mpi_prefix }}{{ cmd_prefix }}{{ operation.cmd }}{{ cmd_suffix }}
{% endfor %}
{% endblock %}
{% block footer %}
{% if parallel %}
wait
{% endif %}
{% endblock %}

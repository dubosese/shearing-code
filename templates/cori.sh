{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash -l
#SBATCH --job-name=screening-{{ id }}
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time={{ walltime|format_timedelta }}
#SBATCH -C knl
#SBATCH --error=error.err
#SBATCH -q regular
#SBATCH --mail-type=ALL
date +"%a-%R" > starttime.txt

module unload gromacs lammps
module load gromacs/2020.1.knl lammps/2018.12.12-knl openmpi

{% endblock %}
{% block body %}
{% set cmd_suffix = cmd_suffix|default('') ~ (' &' if parallel else '') %}
{% for operation in operations %}
{% if operation.directives.nranks and not mpi_prefix %}
{% set mpi_prefix = "" %}
{% endif %}

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

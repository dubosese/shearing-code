{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash -l
#PBS -j oe
#PBS -l nodes=1:ppn=16
#PBS -l walltime={{ walltime|format_timedelta }}
#PBS -q standard 

module load gromacs/2020 lammps
{% endblock %}
{% block body %}
{% set cmd_suffix = cmd_suffix|default('') ~ (' & \nwait' if parallel else '') %}
{% for operation in operations %}
{% if operation.directives.nranks and not mpi_prefix %}
{% set mpi_prefix = " " %}
{% endif %} 
{% if operation.directives.omp_num_threads %}
export OMP_NUM_THREADS={{ operation.directives.omp_num_threads }}
{% endif %}
{{ mpi_prefix }}{{ cmd_prefix }}{{ operation.cmd }}{{ cmd_suffix }}
{% endfor %}
{% endblock %}
{% block footer %}
{% endblock %}

{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash -l
#PBS -j oe
#PBS -l nodes={{ nn }}:ppn=16
#PBS -l walltime={{ walltime|format_timedelta }}
#PBS -q standard 
#PBS -m abe

module load gromacs/5.1.4
{% endblock %}
{% block body %}
{% set cmd_suffix = cmd_suffix|default('') ~ (' &' if parallel else '') %}
{% for operation in operations %}
{% if operation.directives.nranks and not mpi_prefix %}
{% set mpi_prefix = "" %}
{% endif %}

# {{ "%s"|format(operation) }}
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

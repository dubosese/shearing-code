{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash
#SBATCH --job-name=screening
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time={{ walltime|format_timedelta }}
#SBATCH -C haswell
#SBATCH --error=error.err
#SBATCH -q regular

echo working > test.txt

source activate mosdef36
module load gromacs
module load lammps

{% endblock %}

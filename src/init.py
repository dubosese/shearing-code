#!/usr/bin/env python
"""Initialize the project's data space.

Iterates over all defined state points and initializes
the associated job workspace directories."""
import argparse
import logging
import datetime
from copy import deepcopy
from hashlib import sha1

import signac

logging.basicConfig(filename='init.log', filemode='w', level=logging.INFO)


terminal_groups_A = ['acetyl','amino','benzoicacid','biphenyl','carboxyl',
                     'cyano','cyclopropyl','difluoromethyl','dihydroxyphenyl', 
                     'ethylene','fluorophenyl','formyl','hydroxyl', 'isopropyl',
                     'isopropylbenzene','methoxy', 'methyl', 'nitro', 'nitrophenyl',
                     'pentafluorophenyl','perfluoromethyl','phenol','phenyl',
                     'pyrrole', 'toluene']

terminal_groups_B = ['acetyl','amino','benzoicacid','biphenyl','carboxyl',
                     'cyano','cyclopropyl','difluoromethyl','dihydroxyphenyl', 
                     'ethylene','fluorophenyl','formyl','hydroxyl', 'isopropyl',
                     'isopropylbenzene','methoxy', 'methyl', 'nitro', 'nitrophenyl',
                     'pentafluorophenyl','perfluoromethyl','phenol','phenyl',
                     'pyrrole', 'toluene']

terminal_groups_C = ['acetyl','amino','benzoicacid','biphenyl','carboxyl',
                     'cyano','cyclopropyl','difluoromethyl','dihydroxyphenyl', 
                     'ethylene','fluorophenyl','formyl','hydroxyl', 'isopropyl',
                     'isopropylbenzene','methoxy', 'methyl', 'nitro', 'nitrophenyl',
                     'pentafluorophenyl','perfluoromethyl','phenol','phenyl',
                     'pyrrole', 'toluene']

terminal_groups_D = ['acetyl','amino','benzoicacid','biphenyl','carboxyl',
                     'cyano','cyclopropyl','difluoromethyl','dihydroxyphenyl', 
                     'ethylene','fluorophenyl','formyl','hydroxyl', 'isopropyl',
                     'isopropylbenzene','methoxy', 'methyl', 'nitro', 'nitrophenyl',
                     'pentafluorophenyl','perfluoromethyl','phenol','phenyl',
                     'pyrrole', 'toluene']


# Initialize the project
def main(args, random_seed):
    project = signac.init_project(args.project_name)
    logging.info("Init begin: {}".format(datetime.datetime.today()))
    logging.info("Initialized project name")
    statepoints = list()
    # generate the new top and bottom monolayers
    for replication_index in range(args.num_replicas):
        for terminal_group_a in terminal_groups_A:
            for terminal_group_b in terminal_groups_B:
                for terminal_group_c in terminal_groups_C:
                    for terminal_group_d in terminal_groups_D:
                    	if (terminal_group_a == terminal_group_b) and (terminal_group_a == terminal_group_c) and (terminal_group_a == terminal_group_d):
                        	the_statepoint = dict(
                                	# Carbon backbone A length
                                	chainlength_a = args.chain_length_a,
                                	# Carbon backbone B length
                                	chainlength_b = args.chain_length_b,
                                	# Carbon backbone C length
                                	chainlength_c = args.chain_length_c,
					# Carbon backbone D length
					chainlength_d = args.chain_length_d,
                        	        # backbone type 1 
                      	        	backbone_1 = args.backbone_1,
					# backbone type 2
					backbone_2 = args.backbone_2,
                        	        # Number of monolayer chains
       		                        n = args.num_chains,
                          	        # fraction of a vs b on bottom monolayer
                                        fraction_a = args.fraction_a,
					# fraction of c vs d on top monolayer
					fraction_c = args.fraction_c,
                        	        # monolayer pattern type
               		        	pattern_type = args.pattern_type,
                                	# Random seed
                                	seed = random_seed*(replication_index + 1),
                                	# Terminal group chemistries
                                	terminal_groups = tuple((terminal_group_a,
                                        terminal_group_b, terminal_group_c, terminal_group_d)),
                                	# location of atoms in backbone
                                	locations = args.locations)
                        	project.open_job(statepoint=the_statepoint).init()
                        	statepoints.append(the_statepoint)
                        	logging.info(msg="At the statepoint: {}".format(the_statepoint))
    # write statepoints to signac statepoint file
    project.write_statepoints(statepoints=statepoints)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Initialize the data space.")
    parser.add_argument(
	'-project_name',
	type=str,
	help="Signac project name")
    parser.add_argument(
        'random',
        type=str,
        help="A string to generate a random seed.")
    parser.add_argument(
        '-n', '--num-replicas',
        type=int,
        default=1,
        help="Initialize multiple replications.")
    parser.add_argument(
        '-c_a', '--chain-length-a',
        type=int,
        default=17,
        help="Backbone length of chain A (bottom).")
    parser.add_argument(
        '-c_b', '--chain-length-b',
        type=int,
        default=17,
        help="Backbone length of chain B (bottom).")
    parser.add_argument(
        '-c_c', '--chain-length-c',
        type=int,
        default=17,
        help="Backbone length of chain C (top).")
    parser.add_argument(
	'-c_d', '--chain-length-d',
	type=int,
	default=17,
	help="Backbone length of chain D (top).")
    parser.add_argument(
        '-t', '--pattern_type',
        type=str,
        default='random',
        help="Name of pattern, default \'random\'.")
    parser.add_argument(
	'-b_1', '--backbone_1',
	type=str,
	default='methylene',
	help="Type of backbone for backbone 1.")
    parser.add_argument(
	'-b_2', '--backbone_2',
	type=str,
	default='methylene',
	help="Type of backbone for backbone 2.")
    parser.add_argument(
        '-f_a', '--fraction_a',
        type=float,
        default=0.5,
        help="Fraction of chain A on bottom surface vs chain B.")
    parser.add_argument(
	'-f_c', '--fraction_c',
	type=float,
	default=0.5,
	help="Fraction of chain C on top surface vs chain D.")
    parser.add_argument(
	'--num_chains',
	type=int,
	default=100,
	help="Number of chains on each surface.")
    parser.add_argument(
	'-l', '--locations',
	type=list,
	default=[2,5,8,11,14],
	help="Location of backbone heteroatoms.")
    args = parser.parse_args()

    # Generate an integer from the random str.
    try:
        random_seed = int(args.random)
    except ValueError:
        random_seed = int(sha1(args.random.encode()).hexdigest(), 16) % (10 ** 8)

    logging.info("Params:\
                random: {}\
                num-replicas: {}\
                chain_length A: {}\
                chain_length B: {}\
                chain_length C: {}\
		chain_length D: {}\
                pattern_type: {}\
		backbone_1: {}\
		backbone_2: {}".format(
                    random_seed,
                    args.num_replicas,
                    args.chain_length_a,
                    args.chain_length_b,
                    args.chain_length_c,
		    args.chain_length_d,
                    args.pattern_type,
		    args.backbone_1,
		    args.backbone_2
                    ))

    main(args, random_seed)

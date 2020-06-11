#!/bin/bash

$(echo 0 |  gmx_mpi trjconv -s init.gro -f minimize.xtc -o minimize.gro  -b 1.0 -e 1.0  )

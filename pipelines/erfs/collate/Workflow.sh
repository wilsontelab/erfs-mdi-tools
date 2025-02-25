#!/bin/bash

# initialize shared variables
source ${MODULES_DIR}/genome/set_genome_vars.sh

# initialize temporary directory
source ${MODULES_DIR}/utilities/shell/create_shm_dir.sh

# establish binned read counts and insert size distributions per sample
runWorkflowStep 1 collate collate.sh

# assemble external call data
runWorkflowStep 2 calls calls.sh

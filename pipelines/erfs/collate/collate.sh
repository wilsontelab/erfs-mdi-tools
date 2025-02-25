#!/bin/bash

# check for file inputs
echo "checking file inputs"
if [ ! -f ${BRCA2_CTL} ]; then
    echo "ERROR: ${BRCA2_CTL} not found"
    exit 1
fi
if [ ! -f ${BRCA2_SHRNA} ]; then
    echo "ERROR: ${BRCA2_SHRNA} not found"
    exit 1
fi
if [ ! -f ${RAD51_CTL} ]; then
    echo "ERROR: ${RAD51_CTL} not found"
    exit 1
fi
if [ ! -f ${RAD51_SIRNA} ]; then
    echo "ERROR: ${RAD51_SIRNA} not found"
    exit 1
fi

# ensure bam files are indexed
echo "checking bam indices"
if [ ! -f ${BRCA2_CTL}.bai ]; then
    echo "indexing ${BRCA2_CTL}"
    samtools index ${BRCA2_CTL}
    checkPipe
fi
if [ ! -f ${BRCA2_SHRNA}.bai ]; then
    echo "indexing ${BRCA2_SHRNA}"
    samtools index ${BRCA2_SHRNA}
    checkPipe
fi
if [ ! -f ${RAD51_CTL}.bai ]; then
    echo "indexing ${RAD51_CTL}"
    samtools index ${RAD51_CTL}
    checkPipe
fi
if [ ! -f ${RAD51_SIRNA}.bai ]; then
    echo "indexing ${RAD51_SIRNA}"
    samtools index ${RAD51_SIRNA}
    checkPipe
fi

Rscript ${ACTION_DIR}/collate.R
checkPipe

#!/bin/bash

projectDir=/home/boss_lab/Projects/Boss_Projects/Bagde/BagdeHiC_2023-01-11/
dataDir=${projectDir}/data/

for s in 532 536; do
  for c in CD4 CD8; do
    for st in Stim control; do
      sampleDir=${dataDir}/${c}_${s}_${st}/
      mkdir -p ${sampleDir}
      for r in 1 2; do
        gunzip -c $(ls *_${s}-${c}-${st}_S*_R${r}_001.fastq.gz) > ${sampleDir}/${c}_${s}_${st}_R${r}.fastq &
      done
    done
  done
done

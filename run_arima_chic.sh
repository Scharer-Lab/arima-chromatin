#!/bin/bash

projectDir= # Replace, should end in a '/'
genome=hg38 # Replace, options: just hg38 until someone downloads mm10
res=1kb # Replace, options: 1kb, 3kb, 5kb
map=1 # Replace, should we map or just call loops on pre-mapped data? 1 for yes, 0 for no
email="" # Replace
threads=32 # Replace

# Reference paths
pchic=/home/boss_lab/Apps/genomes/captureHiC/
hicref=${pchic}/${genome}/reference/${genome}
digest=${pchic}/${genome}/reference/Digest_${genome}_Arima.txt
chicago_capture_probes=${pchic}/${genome}/${genome}.chicago_probes.bed
chicago_design=${pchic}/${genome}/${res}/
chicago_rmap=${chicago_design}/${genome}_chicago_input_${res}.rmap
chicago_baitmap=${chicago_design}/${genome}_chicago_input_${res}.baitmap

# Tools paths
hicup_path=/home/boss_lab/Apps/HiCUP-0.8.0/
chicago_path=/home/boss_lab/Apps/CHiC/chicagoTools/
bowtie=/home/boss_lab/Apps/anaconda3/bin/bowtie2
PATH=/home/boss_lab/Apps/bedtools-2.25.0/bin/:$PATH # IMPORTANT NOTE: bedtools v2.26 or later isn't compatible with CHiCAGO! Please use v2.25 instead!!!

# Replace with sample names
for sample in SAMPLE1 SAMPLE2 ...
do

    # Make the output directory
    covDir=${projectDir}/coverage/${sample}/
    mkdir -p ${covDir}

    start_time=`date +%s`

    sampleLog=${covDir}/${sample}_sample.${start_time}.log

    fastqDir=${projectDir}/data/${sample}

    # Run the software
    arima-chic -W ${map} -Y 1 -Z 1 -A ${bowtie} -X ${hicref} -d ${digest} \
               -H ${hicup_path} -C ${chicago_path} -b ${chicago_capture_probes} \
               -I ${fastqDir}/${sample}_R1.fastq,${fastqDir}/${sample}_R2.fastq \
               -o ${covDir} -p ${sample} -R ${chicago_rmap} \
               -B ${chicago_baitmap} -D ${chicago_design} \
               -O ${genome} -r ${res} -P 1 -t ${threads} &> ${sampleLog}


    end_time=`date +%s`
    total_time=$(( ${end_time} - ${start_time} ))

    # TODO: Add qc to this email?
    outmessage="Sample ${sample} finished in ${total_time} seconds. Find the results at: ${covDir}"
    if [ ! ${email} -eq ""] && [ ${total_time} -gt 30 ]; then
        mail -s "Arima Capture HiC sample ${sample} has Finished" ${email} <<< ${outmessage}
    fi
done # Sample loop


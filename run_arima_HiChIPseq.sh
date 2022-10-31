#!/bin/bash

projectDir=/PROJECT/DIR/ # Replace, should end in a '/'
downloadDir=/home/boss_lab/Projects/data_downloads/ # Replace, should end in a '/'
oseqDir=/ORIGINAL/PROJECT/DIR/ # Replace
organism="hg38" # Replace
threads=1 # An odd python message occurs when threads > 1. It is likely easily solved...
peaks=/home/boss_lab/Apps/genomes/encode_peaks/${organism}.H3K27ac.encodePeaks.bed # Replace, if necessary
genome=/home/boss_lab/Apps/genomes/bwa/HG38.nohap/HG38.nohap.fa # Replace, if necessary
email=""

for sample in SAMPLE1 SAMPLE2 ... # Replace
do
    start_time=`date +%s`
    logname=pipeline.${sample}.${start_time}.log

    # If this is a resequencing project, we need to combine the fastq files together
    if [ -d ${oseqDir} ]; then
        reseq="_reseq"
        if [ ${sample} == "" ]; then # Replace, if necessary
            match="" # This is what will match the 'X' column of a manifest file
        else
            echo "ERROR: No match found for resequenced sample ${sample}. Exiting"
            exit
        fi
        echo "Sample ${sample} matches with id ${match}." >> ${logname}
    fi

    # Make a sample data directory
    dataDir="${projectDir}/data/${sample}/"
    mkdir -p ${dataDir}  
    
    # Move sample fastq files to the sample's data directory, concatenate all lanes of the forward and backward reads together and then decompress the resulting file
    # TODO: Make sure that the lanes will always be in the same order -- add a debug printout?
    echo "Ensuring fastq files are in their proper spot." >> ${logname}
    if [ ! -f ${dataDir}/${sample}_R1.fastq ] || [ ! -f ${dataDir}/${sample}_R2.fastq ]; then
        echo "Concatenating the forward read files `ls ${downloadDir}/*${match}*R1*.fastq.gz` to ${dataDir}/${sample}_R1.fastq.gz" >> ${logname}
        cat ${downloadDir}/*${match}*R1*.fastq.gz > ${dataDir}/${sample}_R1.fastq.gz # Replace, if necessary because of a difference of naming schema
        echo "Concatenating the backward read files `ls ${downloadDir}/*${match}*R2*.fastq.gz` to ${dataDir}/${sample}_R2.fastq.gz" >> ${logname}
        cat ${downloadDir}/*${match}*R2*.fastq.gz > ${dataDir}/${sample}_R2.fastq.gz # Replace, if necessary because of a difference of naming schema

        echo "Unzipping `ls ${dataDir}/${sample}_R?.fastq.gz`"
        gunzip ${dataDir}/${sample}_R?.fastq.gz # Remember that '?' is the single character wildcard in filename substitution
    fi

    if [ -z ${reseq} ]; then
        echo "Merging old and new libraries" >> ${logname}

        fq1="${dataDir}/${sample}${reseq}_R1.fastq"
        if [ ! -f ${fq1} ]; then
            cat ${oseqDir}/data/${sample}/${sample}_R1.fastq ${dataDir}/${sample}_R1.fastq > ${fq1}
        fi

        fq2="${dataDir}/${sample}${reseq}_R2.fastq"
        if [ ! -f ${fq2} ]; then
            cat ${oseqDir}/data/${sample}/${sample}_R2.fastq ${dataDir}/${sample}_R2.fastq > ${fq2}
        fi
    fi

    # Make the coverage directory
    coverageDir="${projectDir}/coverage/${sample}/"
    mkdir -p ${coverageDir}

    arima_start_time=`date +%s`
    echo "Arima start time: ${arima_start_time}" >> ${logname}

    # Run arima pipeline
    # TODO: Must we use f=0?
    arima -I ${dataDir}/${sample}${reseq} -O ${coverageDir} -m ${peaks} -o ${organism} -b ${genome} -t ${threads} -f 0 &>> ${logname} # Magic -f 0

    # Finish up with logging
    end_time=`date +%s`
    echo "Arima end time: ${end_time}" >> ${logname}
    total_pipeline_time=$(( ${end_time} - ${pipeline_start_time} ))
    echo "Arima time total: ${total_pipeline_time}" >> ${logname}
    total_time=$(( ${end_time} - ${start_time} ))
    echo "Sample time total: ${total_time}" >> ${logname}

    # TODO: More descriptive email? Add qc?
    outmessage="Sample ${sample} has finished in ${total_time} seconds, ${total_pipeline_time} of it running arima itself. Find the results at: ${projectDir}/coverage/"
    echo ${outmessage} >> ${logname}

    if [ ${email} != "" ]; then
        mail -s "Arima sample ${sample} has finished" ${email} <<< ${outmessage}
    else
        write ${USER} <<< ${outmessage}
    fi
done # Sample loop


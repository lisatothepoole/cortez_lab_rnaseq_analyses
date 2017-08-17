#!/usr/bin/python
"""
RNAseq analysis for Cortez Lab
Authors: Lisa Poole and James Pino

This pipeline is optimized for the Lab Mac with IP address 10.105.17.158 on the desk by the glass window.
(Can be optimized for other computers by changing general information and setup under pc.)

Requirements for this experiment - sequencing data files must be in a folder entitled
"fasta" within a folder corresponding to the experiment name; modify general information
in the first part of the script (starting with number of samples).
This protocol is also optimized to be used after gene_expression_pipeline (such that the fasta files have already been
trimmed).

Important: The file names need to be a common sample base followed by a hyphen then sequential numbers (such as LP-1,
LP-2)


"""

import os
import subprocess

# General Information - modify appropriately for each experiment
number_of_samples = '4'  # Valid options are integer values indicating the number of sequenced samples
number_of_samples_add_1 = '5'  # Add 1 to the number_of_samples
species = 'human'  # Valid options for this pipeline are 'mouse' or 'human'
read_type = 'PE'  # Valid options are 'PE' or paired-end sequencing or 'SE' for single-end sequencing
read_length = 150  # Valid options are integer values for the length of reads ordered (50bp, 75bp, 150bp, etc)
sample_suffix = 'fastq'  # suffix second to last in sequencing file, usually fastq, fasta, fa
n_cpus = 8  # Number of threads used to run analysis, more = faster
pc = 'cortez_mac'
experiment_name = '20160816_rnaseq'  # Insert name of experiment that is also the name of folder
desired_search = 'indel'  # Valid options include 'indel' or 'snp' or 'both'

# Program setup and files needed for analysis - should not need to change unless you add an additional organism
if pc == 'cortez_mac':
    # Setting up programs
    samtools = '/usr/local/bin/samtools'
    samstat = '/usr/local/bin/samstat'
    gatk = '/Users/temporary/Sources/GenomeAnalysisTK.jar'
    bwa = '/usr/local/bin/bwa-0.7.15/bwa'
    picard = '/Users/temporary/Sources/picard.jar'

    # experiment specific information
    output_directory = "/Users/temporary/projects/{}".format(experiment_name)
    fasta_directory = "/Users/temporary/projects/{}/fastq".format(experiment_name)

    # Reference files
    if species == 'mouse':
        bwa_index_location = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/version0.7.15/genome_indel'
        reference_genome = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome_indel.fa'
        indel_vcf_file = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/version0.7.15/Mills_and_1000G_gold_standard.indels.hg38.vcf'

    elif species == 'human':
        bwa_index_location = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome'
        reference_genome = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa'
        indel_vcf_file = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/Mills_and_1000G_gold_standard.indels.hg38.vcf'
        snp_vcf_file = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/1000G_phase1.snps.high_confidence.hg38.vcf'

    else:
        print("Error - invalid species. Valid arguments are 'human' or 'mouse'")
        quit()


# Alignment to the genome using BWA aligner
def bwa_alignment(sample_base):
    # -T excludes alignments with score lower than given interval in output
    # -t number of threads
    # -M marks shorter split hits as secondary
    print("Starting BWA alignment for {}".format(sample_base))

    # Creation of folder for aligned files
    if not os.path.exists('{}/BWA_BAM_files'.format(output_directory)):
        os.mkdir('{}/BWA_BAM_files'.format(output_directory))

    path_to_executable = '{} mem'.format(bwa)
    if read_type == 'SE':
        reads = "{}/{}-trimmed.{}".format(fasta_directory, sample_base, sample_suffix)
    if read_type == 'PE':
        reads = '{0}/{1}-trimmed_1.{2} {0}/{1}-trimmed_2.{2}'.format(fasta_directory, sample_base, sample_suffix)
    important_options = '-T 15 -M -t {}'.format(n_cpus)
    path_to_reference = reference_genome
    export_to_file = '> {}/BWA_BAM_files/{}.sam'.format(output_directory, sample_base)
    command = [path_to_executable, important_options, path_to_reference, reads, export_to_file]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
            # if output:
            #         # global_output += output
            # print output.strip()
        rc = process.poll()

    print("Finished aligning {}".format(sample_base))


# Addition of reads groups to SAM file
def bwa_read_group(sample_base):
    print("Starting read group addition for {}".format(sample_base))

    path_to_executable = 'java -jar {}'.format(picard)
    picard_program = "AddOrReplaceReadGroups"
    input_files = 'I={}/BWA_BAM_files/{}.sam'.format(output_directory, sample_base)
    output_files = 'O={}/BWA_BAM_files/{}.rg.sam'.format(output_directory, sample_base)
    necessary_parameters = "RGID={0} RGLB={0} RGPL=ILLUMINA RGPU=ILLUMINA RGSM={0}".format(sample_base)
    command = [path_to_executable, picard_program, input_files, output_files, necessary_parameters]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()

    print("Finished adding read group {}".format(sample_base))
    # Removes SAM file (conserves space)
    os.remove('{}/BWA_BAM_files/{}.sam'.format(output_directory, sample_base))


# Convert to BAM file
def bwa_sam_to_bam(sample_base):
    # -S input is SAM file
    # -b output to BAM file
    print("Starting sam to bam conversion for {}".format(sample_base))

    path_to_executable = '{} view'.format(samtools)
    path_to_samples = '-S -b {}/BWA_BAM_files/{}.rg.sam'.format(output_directory, sample_base)
    output_filename = '-o {}/BWA_BAM_files/{}.bam'.format(output_directory, sample_base)
    threads = '--threads {}'.format(n_cpus)
    command = [path_to_executable, path_to_samples, threads, output_filename]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()

    print("Finished converting sam to bam conversion for {}".format(sample_base))
    # Removes SAM file (conserves space)
    os.remove('{}/BWA_BAM_files/{}.rg.sam'.format(output_directory, sample_base))


# Sort BAM file
def bwa_bam_sort(sample_base):
    print("Start sorting {}".format(sample_base))

    path_to_executable = '{} sort'.format(samtools)
    path_to_samples = '{}/BWA_BAM_files/{}.bam'.format(output_directory, sample_base)
    output_filename = '-o {}/BWA_BAM_files/{}.sorted.bam'.format(output_directory, sample_base)
    threads = '--threads {}'.format(n_cpus)
    command = [path_to_executable, threads, output_filename, path_to_samples]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()

    print("Finished sorting {}".format(sample_base))
    # Remove BAM file
    os.remove('{}/BWA_BAM_files/{}.bam'.format(output_directory, sample_base))


# Index BWA BAM file for viewing in IGV
def bwa_index(sample_base):
    print("Starting index for {}".format(sample_base))

    path_to_executable = '{} index'.format(samtools)
    path_to_samples = '{}/BWA_BAM_files/{}.sorted.bam'.format(output_directory, sample_base)
    command = [path_to_executable, path_to_samples]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()

    print("Finished with index for {}".format(sample_base))


# Quality Control with SAMSTAT
def bwa_samstat_analysis(sample_base):
    print("Starting samstat analysis for {}".format(sample_base))

    path_to_executable = samstat
    path_to_samples = '{}/BWA_BAM_files/{}.sorted.bam'.format(output_directory, sample_base)
    command = [path_to_executable, path_to_samples]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()

    print("Done with samstat analysis for {}".format(sample_base))
    # Move SAMSTAT analysis to quality control folder
    os.rename('{}/BWA_BAM_files/{}.sorted.bam.samstat.html'.format(output_directory, sample_base),
              '{}/quality_control/{}.bwa.sorted.bam.samstat.html'.format(output_directory, sample_base))


# INDEL/SNP specific alignments
def gatk_intervals(sample_base, entity_searched):
    print("Starting gatk intervals for {}".format(sample_base))

    path_to_executable = "java -jar {}".format(gatk)
    gatk_program = '-T RealignerTargetCreator'
    path_to_reference = "-R {}".format(reference_genome)
    input_files = '-I {}/BWA_BAM_files/{}.sorted.bam'.format(output_directory, sample_base)
    output_file = '-o {}/BWA_BAM_files/{}.{}.intervals'.format(output_directory, entity_searched, sample_base)

    if species == 'human':
        if entity_searched == 'indel':
            path_to_vcf = "--known {}".format(indel_vcf_file)
        elif entity_searched == 'snp':
            path_to_vcf = "--known {}".format(snp_vcf_file)

    elif species == 'mouse':
        if entity_searched == 'indel':
            path_to_vcf = "--known {}".format(indel_vcf_file)
        else:
            path_to_vcf = ''

    command = [path_to_executable, gatk_program, path_to_reference, input_files, path_to_vcf, output_file]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Finished with gatk intervals for {}".format(sample_base))


def gatk_realignment(sample_base, entity_searched):
    print("Starting gatk realignment for {}".format(sample_base))

    path_to_executable = "java -jar {}".format(gatk)
    gatk_program = '-T IndelRealigner'
    path_to_reference = "-R {}".format(reference_genome)
    options = '--maxReadsForRealignment 999999 --maxReadsInMemory 999999'
    input_files = '-I {}/BWA_BAM_files/{}.sorted.bam'.format(output_directory, sample_base)
    intervals = '-targetIntervals {}/BWA_BAM_files/{}.{}.intervals'.format(output_directory, sample_base,
                                                                           entity_searched)
    output_file = '-o {}/BWA_BAM_files/{}.{}.realigned.bam'.format(output_directory, sample_base, entity_searched)

    command = [path_to_executable, gatk_program, path_to_reference, input_files, intervals, options, output_file]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Finished with gatk realignment for {}".format(sample_base))


def gatk_recalibration(sample_base, entity_searched):
    print("Starting gatk recalibration for {}".format(sample_base))

    path_to_executable = "java -jar {}".format(gatk)
    gatk_program = '-T BaseRecalibrator'
    path_to_reference = "-R {}".format(reference_genome)
    options = '-l INFO'
    input_files = '-I {}/BWA_BAM_files/{}.{}.realigned.bam'.format(output_directory, sample_base, entity_searched)
    output_file = '-o {}/BWA_BAM_files/{}.{}.recal.table'.format(output_directory, sample_base, entity_searched)

    if species == 'human':
        if entity_searched == 'indel':
            path_to_vcf = "-knownSites {}".format(indel_vcf_file)
        elif entity_searched == 'snp':
            path_to_vcf = "--knownSites {}".format(snp_vcf_file)

    elif species == 'mouse':
        if entity_searched == 'indel':
            path_to_vcf = "-knownSites {}".format(indel_vcf_file)
        else:
            path_to_vcf = ''

    command = [path_to_executable, gatk_program, path_to_reference, input_files, options, path_to_vcf, output_file]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()

    print("Finished with gatk indel recalibration for {}".format(sample_base))


def gatk_realign_recal(sample_base, entity_searched):
    print("Starting gatk realign recal for {}".format(sample_base))

    path_to_executable = "java -jar {}".format(gatk)
    gatk_program = '-T PrintReads'
    path_to_reference = "-R {}".format(reference_genome)
    input_files = '-I {}/BWA_BAM_files/{}.{}.realigned.bam'.format(output_directory, sample_base, entity_searched)
    options = '-BQSR {}/BWA_BAM_files/{}.{}.recal.table'.format(output_directory, sample_base, entity_searched)
    output_file = '-o {}/BWA_BAM_files/{}.{}.realigned.recal.bam'.format(output_directory, sample_base, entity_searched)

    command = [path_to_executable, gatk_program, path_to_reference, input_files, options, output_file]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()

    print("Finished with gatk realign recal for {}".format(sample_base))


def mark_dup(sample_base, entity_searched):
    print("Starting mark dup for {}".format(sample_base))

    path_to_executable = "java -jar {}".format(picard)
    picard_program = 'MarkDuplicates'
    input_files = 'I={}/BWA_BAM_files/{}.{}.realigned.recal.bam'.format(output_directory, sample_base, entity_searched)
    output_file = 'O={}/BWA_BAM_files/{}.{}.realigned.recal.dupmarked.bam'.format(output_directory, sample_base,
                                                                                  entity_searched)
    options = 'M={}/BWA_BAM_files/{}-{}-marked_dup_metrics.txt'.format(output_directory, sample_base, entity_searched)

    command = [path_to_executable, picard_program, input_files, output_file, options]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()

    print("Finished with mark dup for {}".format(sample_base))


def dup_index(sample_base, entity_searched):
    print("Starting indel dup index for {}".format(sample_base))
    path_to_executable = '{} index'.format(samtools)
    path_to_samples = '{}/BWA_BAM_files/{}.{}.realigned.recal.dupmarked.bam'.format(output_directory, sample_base,
                                                                                    entity_searched)

    command = [path_to_executable, path_to_samples]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Done with dup index for {}".format(sample_base))


def final_entity_search(sample_base, entity_searched):
    print("Starting final search for {}".format(sample_base))
    path_to_executable = "java -jar {}".format(gatk)
    gatk_program = '-T UnifiedGenotyper -l INFO'
    path_to_reference = "-R {}".format(reference_genome)
    options = '-A Coverage -A AlleleBalance -G Standard -stand_call_conf 50.0 -stand_emit_conf 10.0 -mbq 20 -deletions 0.05 -dcov 1000'
    input_files = '-I {}/BWA_BAM_files/{}.{}.realigned.recal.dupmarked.bam'.format(output_directory, sample_base,
                                                                                   entity_searched)
    output_file = '--out {0}/{1}.{2}.vcf -metrics {0}/{1}.{2}.outmetrics.txt'.format(output_directory, sample_base,
                                                                                     entity_searched)

    if entity_searched == 'indel':
        search_item = '-glm INDEL'
    elif entity_searched == 'snp':
        search_item = '-glm SNP'

    if species == 'human':
        if entity_searched == 'indel':
            path_to_vcf = "-D {}".format(indel_vcf_file)
        elif entity_searched == 'snp':
            path_to_vcf = "-D {}".format(snp_vcf_file)

    elif species == 'mouse':
        if entity_searched == 'indel':
            path_to_vcf = "-D {}".format(indel_vcf_file)
        else:
            path_to_vcf = ''

    command = [path_to_executable, gatk_program, path_to_reference, input_files, path_to_vcf, output_file, options,
               search_item]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Finished with final search for {}".format(sample_base))


def bwa_setup(sample_base):
    bwa_alignment(sample_base)
    bwa_read_group(sample_base)
    bwa_sam_to_bam(sample_base)
    bwa_bam_sort(sample_base)
    bwa_index(sample_base)
    bwa_samstat_analysis(sample_base)


def entity_search(sample_base, entity_searched):
    gatk_intervals(sample_base, entity_searched)
    gatk_realignment(sample_base, entity_searched)
    gatk_recalibration(sample_base, entity_searched)
    gatk_realign_recal(sample_base, entity_searched)
    mark_dup(sample_base, entity_searched)
    dup_index(sample_base, entity_searched)
    final_entity_search(sample_base, entity_searched)


def bwa_entity_search(sample_base):
    if desired_search == 'indel':
        entity_search(sample_base, 'indel')

    elif bwa_entity_search == 'snp':
        entity_search(sample_base, 'snp')

    elif bwa_entity_search == 'both':
        entity_search(sample_base, 'indel')
        entity_search(sample_base, 'snp')


def final_bwa_search(sample_base):
    bwa_setup(sample_base)
    bwa_entity_search(sample_base)

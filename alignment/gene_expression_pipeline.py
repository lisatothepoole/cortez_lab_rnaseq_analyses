#!/usr/bin/python

import os
import subprocess

# General Information - modify appropriately for each experiment
number_of_samples = 10  # Valid options are integer values indicating the number of sequenced samples
number_of_samples_add_1 = 11  # Add 1 to the number_of_samples
species = 'mouse'  # Valid options for this pipeline are 'mouse' or 'human'
# read_type = 'PE'  # Valid options are 'PE' or paired-end sequencing or 'SE' for single-end sequencing
# read_length = 75  # Valid options are integer values for the length of reads ordered (50bp, 75bp, 150bp, etc)
# sample_suffix = 'fastq'  # suffix second to last in sequencing file, usually fastq, fasta, fa
# compression_suffix = 'gz'  # suffix at the end of sequencing file
# n_cpus = 8  # Number of threads used to run analysis, more = faster
pc = 'cortez_mac' # should be cortez_mac for the computer by the window, modify if a new computer is added
# experiment_name = 'alex_rnaseq'  # Insert name of experiment that is also the name of folder

# Program setup and files needed for analysis - should not need to change unless you add an additional organism
if pc == 'cortez_mac':
    # Setting up programs
    flexbar = '/usr/local/bin/flexbar_v2.5_macosx/flexbar'
    hisat2 = '/usr/local/bin/hisat2-2.0.4/hisat2'
    samtools = '/usr/local/bin/samtools'
    featurecounts = '/usr/local/bin/featureCounts'
    samstat = '/usr/local/bin/samstat'
    cufflinks = '/Users/temporary/Sources/cufflinks-2.1.1.OSX_x86_64/cufflinks'
    adaptors = '/Users/temporary/genomes/adapters/illumina_truseq.fasta'

    # # experiment specific information
    # output_directory = "/Users/temporary/projects/{}".format(experiment_name)
    # fasta_directory = "/Users/temporary/projects/{}/fastq".format(experiment_name)

    # Reference files
    if species == 'mouse':
        transcripts = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
        hisat2_index = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/HISAT2_index/mm10.genome'
        reference_genome = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome_indel.fa'

    elif species == 'human':
        transcripts = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'
        hisat2_index = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/HISAT2Index/hisat_genome'
        reference_genome = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa'

    else:
        print("Error - invalid species. Valid arguments are 'human' or 'mouse'")
        quit()


# Sequencing Adaptor Trimming via Flexbar
def flexbar_trim(sample_base, read_type, read_length, adaptors, sample_suffix, compression_suffix, n_cpus, fasta_directory):
    # Flexbar options
    # '-a' = adaptor sequences (fasta format)
    # '-n' = number of threads
    # '-u' = max uncalled bases for each read to pass filtering
    # '-m' = min read length to remain after filtering/trimming
    # '-t' = prefix for output file names
    # '-ao' = adapter min overlap
    # '-ae' = adapter trim end
    # Options chosen below are based on suggestions from the Pietenpol Lab, but feel free to change
    print("Start trimming {}".format(sample_base))

    path_to_executable = flexbar
    suffix_for_output = '-t {}/{}-trimmed'.format(fasta_directory, sample_base)
    adaptor_trim_end = '-ae ANY'
    adaptor_overlap = '-ao 5'
    path_to_adaptors = "-a {}".format(adaptors)
    number_max_uncalled_bases_pass = '-u {}'.format(read_length)
    main_read_length_to_remain = '-m 18'
    threads = '-n {}'.format(n_cpus)

    if read_type not in ('SE', 'PE'):
        print("Error. Invalid read type. Valid options are SE or PE Aborting.")
        quit()

    elif read_type == 'SE':
        reads = '-r {}/{}.{}.{}'.format(fasta_directory, sample_base, sample_suffix, compression_suffix)

    elif read_type == 'PE':
        reads = ' -r {0}/{1}_1.{2}.{3} -p {0}/{1}_2.{2}.{3}'.format(fasta_directory, sample_base, sample_suffix,
                                                                    compression_suffix)

    command = [path_to_executable, reads, path_to_adaptors, threads, suffix_for_output, adaptor_overlap,
               adaptor_trim_end, number_max_uncalled_bases_pass, main_read_length_to_remain]
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
    print("Done trimming {}".format(sample_base))


# Alignment to the genome using HISAT2
def hisat2_alignment(sample_base, n_cpus, read_type, sample_suffix, output_directory, fasta_directory):
    print("Starting alignment of {}".format(sample_base))

    # Creation of folder to hold the BAM files
    if not os.path.exists('{}/BAM_files'.format(output_directory)):
        os.mkdir('{}/BAM_files'.format(output_directory))

    # Options for HISAT2 alignment
    path_to_executable = hisat2
    output_name = '-S {}/BAM_files/{}.sam'.format(output_directory, sample_base)
    threads = '-p {}'.format(n_cpus)
    indices = '-x {}'.format(hisat2_index)

    if read_type == 'SE':
        reads = '-U {}/{}-trimmed.{}'.format(fasta_directory,
                                             sample_base, sample_suffix)
    elif read_type == 'PE':
        reads = '-1 {0}/{1}-trimmed_1.{2} -2 {0}/{1}-trimmed_2.{2}'.format(fasta_directory, sample_base, sample_suffix)

    command = [path_to_executable, threads, indices, reads, output_name]
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
    print("Finished aligning {}".format(sample_base))


# Conversion from SAM file to BAM file
def sam_to_bam(sample_base, n_cpus, output_directory):
    # -S input is SAM file
    # -b output to BAM file
    print("Start SAM to BAM conversion for {}".format(sample_base))

    path_to_executable = '{} view'.format(samtools)
    path_to_samples = '-S -b {}/BAM_files/{}.sam'.format(output_directory, sample_base)
    output_filename = '-o {}/BAM_files/{}.bam'.format(output_directory, sample_base)
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
    print("Finished SAM to BAM conversion for {}".format(sample_base))


# Sorting BAM file
def bam_sort(sample_base, n_cpus, output_directory):
    print("Start sorting {}".format(sample_base))

    path_to_executable = '{} sort'.format(samtools)
    path_to_samples = '{}/BAM_files/{}.bam'.format(output_directory, sample_base)
    output_filename = '-o {}/BAM_files/{}.sorted.bam'.format(output_directory, sample_base)
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


# Indexing BAM file for viewing with IGV
def bam_index(sample_base, n_cpus, output_directory):
    print("Start indexing {}".format(sample_base))

    path_to_executable = '{} index'.format(samtools)
    path_to_samples = '{}/BAM_files/{}.sorted.bam'.format(output_directory, sample_base)
    threads = '--threads {}'.format(n_cpus)
    command = [path_to_executable, threads, path_to_samples]
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
    print("Finished indexing {}".format(sample_base))


# SAMSTAT quality check of alignment
def samstat_analysis(sample_base, output_directory):
    print("Start SAMSTAT check for {}".format(sample_base))

    # Creation of folder to contain quality control analyses
    if not os.path.exists('{}/quality_control'.format(output_directory)):
        os.mkdir('{}/quality_control'.format(output_directory))

    path_to_executable = samstat
    path_to_samples = '{}/BAM_files/{}.sorted.bam'.format(output_directory, sample_base)
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
    print("Done with SAMSTAT check for {}".format(sample_base))

    # Move SAMSTAT files to quality control folder
    os.rename('{}/BAM_files/{}.sorted.bam.samstat.html'.format(output_directory, sample_base),
              '{}/SAMSTAT_analysis/{}.hisat2.sorted.bam.samstat.html'.format(output_directory, sample_base))


# Delete excessive files produced during pipeline
def excess_file_cleanup(sample_base, read_type, sample_suffix, output_directory, fasta_directory):
    os.remove('{}/BAM_files/{}.sam'.format(output_directory, sample_base))
    os.remove('{}/BAM_files/{}.bam'.format(output_directory, sample_base))

    if read_type == 'SE':
        os.remove('{}/{}-trimmed.{}'.format(fasta_directory,
                                             sample_base, sample_suffix))
    elif read_type == 'PE':
        os.remove('{0}/{1}-trimmed_1.{2}'.format(fasta_directory, sample_base, sample_suffix))
        os.remove('{0}/{1}-trimmed_2.{2}'.format(fasta_directory, sample_base, sample_suffix))


# Completion of all steps of Aligning to RNAseq data
def rnaseq_expression_alignment():
    flexbar_trim()
    hisat2_alignment()
    sam_to_bam()
    bam_sort()
    bam_index()
    samstat_analysis()


# FeatureCounts - align reads to genes to give expression level
def featurecounts_analysis(sample, experiment_name, output_directory, n_cpus, read_type):
    # FeatureCounts options
    # -a annotation file
    # -o name of output file with read counts
    # -p fragments counted instead of reads (paired end specific)
    # -B count only read pairs that have both ends successfully aligned only
    # -C don't count reads that have two ends mapping to different chromosomes or mapping to same chromosome but on different strands
    # -Q minimum mapping quality score a read must satisfy in order to be counted
    # -t feature type in GTF file, default is "exon"
    # -g specify attribute in GTF file, default is gene_id
    print('Begin featureCounts analysis of all samples')

    path_to_executable = featurecounts
    annotation_file = "-a {}".format(transcripts)
    output_name = "-o {}_gene_counts.txt".format(experiment_name)
    gtf_feature = '-t exon'
    gtf_attibute = '-g gene_id'
    quality_score = '-Q 30'
    out_string = ''
    for i in range(1, number_of_samples_add_1):
        input_files = '{}/BAM_files/{}-{}.sorted.bam'.format(output_directory, sample, i)
        out_string += input_files
    if read_type == 'SE':
        important_options = "-T {}".format(n_cpus)
    elif read_type == 'PE':
        important_options = "-p -B -C -T {}".format(n_cpus)


    command = [path_to_executable, important_options, gtf_feature, gtf_attibute, quality_score, annotation_file,
               output_name, out_string]
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
    print("Done running {}".format(sample))


def fpkm_tpm_calculation(sample_base, output_directory, n_cpus, read_type):
    print('Begin cufflinks FPKM and TPM calculations')
    # Creation of folder to contain fpkm_analysis
    if not os.path.exists('{}/quantification'.format(output_directory)):
        os.mkdir('{}/quantification'.format(output_directory))

    os.mkdir('{}/quantification/{}'.format(output_directory, sample_base))

    path_to_executable = cufflinks
    output = "-o {}/quantification/{}".format(output_directory, sample_base)
    threads = '-p {}'.format(n_cpus)
    annotation_file = "-g {}".format(transcripts)
    if read_type == 'SE':
        library_type = "--library-type fr-unstranded"
    elif read_type == 'PE':
        library_type = "--library-type fr-firststrand"

    input_file = '{}/BAM_files/{}.sorted.bam'.format(output_directory, sample_base)

    command = [path_to_executable, output, threads, annotation_file, library_type, input_file]
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
    print("Done running {}".format(sample_base))


def run_all(sample):
    for i in range(1, number_of_samples_add_1):
        rnaseq_expression_alignment('{}-{}'.format(sample, i))
        fpkm_tpm_calculation('{}-{}'.format(sample, i))
    featurecounts_analysis(sample)
    for i in range(1, number_of_samples_add_1):
        excess_file_cleanup('{}-{}'.format(sample, i))
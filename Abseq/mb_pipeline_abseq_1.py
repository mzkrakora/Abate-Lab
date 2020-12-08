'''

master script for running pipeline to export mission bio data to sqlite3 db
mission bio single-cell abseq experiment #1
written by ben 9.17.2018

requirements:
-paired-end sequencing for both Ab tags and panel
-all Ab tag data has associated panel data
-concatenated Ab/panel files

'''

import os
import subprocess
import sys
import sqlite3

class AbseqSample(object):
    # class for storing metadata for each mission bio abseq sample (sample = ab tags + panel)

    def __init__(self,
                 sample_name,
                 db,
                 sample_num,
                 panel_r1,
                 panel_r2,
                 panel_r1_filtered,
                 panel_r2_filtered,
                 ab_r1,
                 ab_r2,
                 ab_r1_filtered,
                 ab_r2_filtered):

        self.sample_name = sample_name          # sample base name
        self.sample_num = sample_num            # number identifying sample (or tube)

        self.db = db                            # sample database path

        self.panel_r1 = panel_r1                # panel R1 fastq path
        self.panel_r2 = panel_r2                # panel R2 fastq path
        self.panel_r1_filtered = panel_r1_filtered      # filtered panel R1 fastq path
        self.panel_r2_filtered = panel_r2_filtered      # filtered panel R2 fastq path

        self.ab_r1 = ab_r1                      # ab R1 fastq path
        self.ab_r2 = ab_r2                      # ab R2 fastq path
        self.ab_r1_filtered = ab_r1_filtered    # filtered ab R1 fastq path
        self.ab_r2_filtered = ab_r2_filtered    # filtered ab R2 fastq path


    def filter_valid_reads(self, r1, r2, r1_filtered, r2_filtered):
        # from R1 files, extract barcode sequences and export to db

        # filter R1 and R2 files to only keep reads with correct barcode structure in R1
        filter_cmd = 'cutadapt -g %s %s %s -o %s -p %s --wildcard-file=%s --discard-untrimmed --quiet' \
                     % (cell_barcode_structure,
                        r1,
                        r2,
                        r1_filtered,
                        r2_filtered,
                        os.path.expanduser(wildcard_file + str(self.sample_num)))

        process = subprocess.Popen(filter_cmd, shell=True)

        return process

    def extract_cell_barcodes(self, table_name):
        # from cutadapt wildcard file, extract valid barcodes and export to db

        barcodes_cmd = 'python %simport_cell_barcodes_db.py %s %s %s -table %s -bar_suffix %s' % \
                       (code_dir,
                        self.db,
                        cell_barcode_csv,
                        os.path.expanduser(wildcard_file + str(self.sample_num)),
                        table_name,
                        self.sample_num)

        process = subprocess.Popen(barcodes_cmd, shell=True)

        return process

    def align_panel(self):
        # align the panel to the bowtie2 index and export to db
        # only output alignments with MAPQ >= 2

        align_cmd = '(bowtie2 -x %s -1 %s -2 %s -p 32 --no-mixed --no-discordant) 2>%s%s_bowtiestats.txt' \
                    ' | samtools view -bq2 -@ 24' \
                    ' | samtools sort -n -@ 24' \
                    ' | python ~/code/missionbio/import_aln_db.py --stdin -db %s -table panel_alignments' \
                     % (os.path.expanduser(bt2_ref),
                        self.panel_r1_filtered,
                        self.panel_r2_filtered,
                        os.path.expanduser(bt2_out_dir),
                        self.sample_name,
                        self.db)

        process = subprocess.Popen(align_cmd, shell=True)

        return process


    def extract_abs(self):
        # extract ab tag barcodes and umis and save to db

        extract_cmd = 'python %simport_ab_tags_db.py %s %s %s %s %s %d %d %d %d %d %d -table ab_tags' \
                      % (code_dir,
                         self.ab_r2_filtered,
                         self.db,
                         ab_5_handle,
                         ab_3_handle,
                         ab_barcode_csv,
                         ab_bar_coord[0],
                         ab_bar_coord[1],
                         ab_umi1_coord[0],
                         ab_umi1_coord[1],
                         ab_umi2_coord[0],
                         ab_umi2_coord[1])

        process = subprocess.Popen(extract_cmd, shell=True)

        return process


    def count_umis(self, n_children, umi_grouping_methods):
        # count umis using umi-tools and save to db

        count_cmd = 'python %sgroup_umis_by_cell_db.py %s -methods %s ' \
                    '-t %d -in-table ab_tags -out-table umi_counts' \
                    % (code_dir,
                       self.db,
                       ' '.join(umi_grouping_methods),
                       n_children)

        process = subprocess.Popen(count_cmd, shell=True)

        return process

    def clean_db(self):
        # performs a series of joins on the resulting db to create summary tables

        sql = sqlite3.connect(self.db, isolation_level='Exclusive')
        s = sql.cursor()  # create cursor

        # join tables to create panel table
        s.execute(
            '''create table panel
                as select
                panel_cell_barcodes.read_id, panel_cell_barcodes.corr_cell_barcode as cell_barcode, ref_1, ref_2, is_diff_ref
                from panel_cell_barcodes
                left join
                panel_alignments
                on
                panel_cell_barcodes.read_id = panel_alignments.read_id''')

        # close the db connection
        sql.commit()  # commit db changes to disk one last time
        sql.close()  # close the db

def file_summary(samples):
    # displays experiment files and checks for existing db files

    # print summary of all samples identified, find existing db files
    existing_dbs = []
    for sample in samples:

        if os.path.isfile(sample.db):
            existing_dbs += [sample.db]
        else:
            subprocess.Popen('touch %s' % sample.db, shell=True)

        s = vars(sample)

        for item in s:
            print item,': ' ,s[item]
        print '\n'

    print '%d samples identified.\n' % len(samples)

    # print existing db files (if any exist)
    if len(existing_dbs) > 0:

        print 'Found existing database files:'

        for d in existing_dbs:
            print d
        print '\n'

        # ask user to delete existing files
        yes = {'yes', 'y'}
        no = {'no', 'n'}

        while True:

            choice = raw_input('Do you want to delete existing db files? ("yes" to delete, "no" to exit) ').lower()

            if choice in yes:
                print 'Deleting existing files...'
                for d in existing_dbs:
                    os.remove(d)
                    subprocess.Popen('touch %s' % d, shell=True)
                break

            elif choice in no:
                print 'Exiting...'
                raise SystemExit

            else:
                sys.stdout.write('''Please respond with 'yes' or 'no'.\n''')

def get_fastq_filenames(path_to_fastq, paired=True):
    # gets fastq filenames in a given directory and runs some simple checks
    # assumes files are compressed with .gz extensions

    R1_files = []
    R2_files = []

    for file in os.listdir(path_to_fastq):

        if file.endswith('.fastq.gz'):
            R = file.split('_')[-2]  # R1 or R2

            if R == 'R1':
                R1_files += [path_to_fastq + file]

            elif R == 'R2' and paired:
                R2_files += [path_to_fastq + file]

            else:
                print 'Unexpected filename structure. Exiting...'
                raise SystemExit

    if len(R1_files) != len(R2_files) and paired:
        print 'Unpaired FASTQ files exist! Check input files.'
        raise SystemExit

    R1_files.sort()
    R2_files.sort()

    return R1_files, R2_files

def merge_dbs(samples):
    # merges all sample dbs into a single db with concatenated records

    merged_db_name = os.path.expanduser(db_out_dir) + sample_basename + '.db'   # new db name

    # merge dbs using external script
    to_merge = [sample.db for sample in samples]
    merge_cmd = 'python %sdb_merge_script.py %s -other_dbs %s' % (code_dir, to_merge[0], ' '.join(to_merge[1:]))
    process = subprocess.Popen(merge_cmd, shell=True)
    process.communicate()

    # rename new db and remove old ones
    os.rename(to_merge[0], merged_db_name)
    for db in to_merge[1:]:
        os.remove(db)

def wait(processes):
    # waits for processes to finish
    return [process.communicate() for process in processes]


if __name__ == "__main__":

    # global experiment variables - modify for each experiment

    sample_basename = 'abseq1'
    base = '~/missionbio/' + sample_basename + '/'

    fastq_dir = base + 'fastq/'  # raw fastq file directory (files should be .gz zipped)
    db_out_dir = base + 'dbs/'  # directory to output sqlite db files
    bt2_out_dir = base + 'bowtiestats/'  # directory to output bowtiestats files for panel alignment

    code_dir = '~/code/missionbio/'     # code directory

    bt2_ref = '~/missionbio/bt2/AML_amplicons'  # bowtie2 index location

    ab_barcode_csv = base + 'barcodes/ab_barcodes_abseq1.csv'  # ab barcodes csv path
    cell_barcode_csv = base + 'barcodes/mb_cell_barcodes.csv'  # cell barcodes csv path

    # structure of cell barcode to pass to cutadapt and location of barcode
    cell_barcode_structure = '^NNNNNNNNGAGTGATTGCTTGTGACGCCTTNNNNNNNNCGATGACG'
    wildcard_file = base + 'fastq/temp/wildcard.temp'     # temporary wildcard file

    # locations of umi and ab barcode in trimmed ab read
    ab_5_handle, ab_3_handle = 'GTAAGTGCTGATCTTGG', 'AAGCTTGTTTCTGTGCACTGAG'  # ab tag handle sequences
    ab_bar_coord = [6, 14]
    ab_umi1_coord = [0, 6]
    ab_umi2_coord = [14, 20]

    umi_grouping_methods = ['all']      # methods used to count umis

    print '''

                                Welcome to AbSeq!

          ___                         ___           ___                   
         /  /\         _____         /  /\         /  /\          ___     
        /  /::\       /  /::\       /  /:/_       /  /:/_        /  /\    
       /  /:/\:\     /  /:/\:\     /  /:/ /\     /  /:/ /\      /  /::\   
      /  /:/~/::\   /  /:/~/::\   /  /:/ /::\   /  /:/ /:/_    /  /:/\:\  
     /__/:/ /:/\:\ /__/:/ /:/\:| /__/:/ /:/\:\ /__/:/ /:/ /\  /  /:/~/::\ 
     \  \:\/:/__\/ \  \:\/:/~/:/ \  \:\/:/~/:/ \  \:\/:/ /:/ /__/:/ /:/\:\ 
      \  \::/       \  \::/ /:/   \  \::/ /:/   \  \::/ /:/  \  \:\/:/__\/
       \  \:\        \  \:\/:/     \__\/ /:/     \  \:\/:/    \  \::/     
        \  \:\        \  \::/        /__/:/       \  \::/      \__\/      
         \__\/         \__\/         \__\/         \__\/             


Beginning proteogenomic pipeline...

####################################################################################
# Step 1: identify samples
####################################################################################
    '''

    # get all fastq filenames
    R1_files, R2_files = get_fastq_filenames(os.path.expanduser(fastq_dir))

    # store sample info in AbseqSample objects
    samples = []

    for i in range(0, len(R1_files)/2):

        # assign filenames to sample types
        # note: using alphabetization pattern which may not exist in future
        
        ab_r1 = R1_files[i]
        ab_r2 = R2_files[i]
        panel_r1 = R1_files[i + len(R1_files)/2]
        panel_r2 = R2_files[i + len(R1_files)/2]

        sample_num = i + 1

        # set file locations and append to sample object
        db = os.path.expanduser(db_out_dir) + sample_basename + '-' + str(sample_num) + '.db'

        panel_r1_filtered = panel_r1.split('.fastq.gz')[0] + '_filtered.fastq.gz'
        panel_r1_filtered = '/'.join(panel_r1_filtered.split('/')[:-1] + ['filtered'] + [panel_r1_filtered.split('/')[-1]])

        panel_r2_filtered = panel_r2.split('.fastq.gz')[0] + '_filtered.fastq.gz'
        panel_r2_filtered = '/'.join(panel_r2_filtered.split('/')[:-1] + ['filtered'] + [panel_r2_filtered.split('/')[-1]])

        ab_r1_filtered = ab_r1.split('.fastq.gz')[0] + '_filtered.fastq.gz'
        ab_r1_filtered = '/'.join(ab_r1_filtered.split('/')[:-1] + ['filtered'] + [ab_r1_filtered.split('/')[-1]])

        ab_r2_filtered = ab_r2.split('.fastq.gz')[0] + '_filtered.fastq.gz'
        ab_r2_filtered = '/'.join(ab_r2_filtered.split('/')[:-1] + ['filtered'] + [ab_r2_filtered.split('/')[-1]])

        samples.append(AbseqSample(sample_basename,
                                 db,
                                 sample_num,
                                 panel_r1,
                                 panel_r2,
                                 panel_r1_filtered,
                                 panel_r2_filtered,
                                 ab_r1,
                                 ab_r2,
                                 ab_r1_filtered,
                                 ab_r2_filtered))

    # display sample summary and check files for existing databases
    file_summary(samples)

    print '''
####################################################################################
# Step 2: extract cell barcodes from panel reads
####################################################################################
'''

    # for panel reads, filter reads with valid read structure
    filter_reads = [
        sample.filter_valid_reads(sample.panel_r1, sample.panel_r2, sample.panel_r1_filtered, sample.panel_r2_filtered)
        for sample in samples]
    # wait for all processes to finish before continuing
    wait(filter_reads)

    # extract valid panel cell barcodes and export to db
    cell_barcodes = [sample.extract_cell_barcodes('panel_cell_barcodes') for sample in samples]
    # wait for all processes to finish before continuing
    wait(cell_barcodes)

    # delete temporary files
    for sample in samples:
        os.remove(os.path.expanduser(wildcard_file + str(sample.sample_num)))

    print '''
####################################################################################
# Step 3: extract cell barcodes from antibody reads
####################################################################################
'''

    # for ab reads, filter reads with valid read structure
    filter_reads = [
        sample.filter_valid_reads(sample.ab_r1, sample.ab_r2, sample.ab_r1_filtered, sample.ab_r2_filtered)
        for sample in samples]
    # wait for all processes to finish before continuing
    wait(filter_reads)

    # extract valid ab cell barcodes and export to db
    cell_barcodes = [sample.extract_cell_barcodes('ab_cell_barcodes') for sample in samples]
    # wait for all processes to finish before continuing
    wait(cell_barcodes)

    # delete temporary files
    for sample in samples:
        os.remove(os.path.expanduser(wildcard_file + str(sample.sample_num)))

    print '''
####################################################################################
# Step 4: align amplicon panel to reference genome
####################################################################################
'''

    # align panel and export to db
    align_panel_samples = [sample.align_panel() for sample in samples]
    # wait for all processes to finish before continuing
    wait(align_panel_samples)

    print 'All panel samples aligned to reference and saved to database.'

    print '''
####################################################################################
# Step 5: extract antibody barcodes and UMIs
####################################################################################
'''

    # next, extract all ab barcodes and umis and save to db
    # includes error correction for ab barcodes
    extract_ab_tags = [sample.extract_abs() for sample in samples]
    # wait for all processes to finish before continuing
    wait(extract_ab_tags)

    print 'Ab barcodes and umis extracted from all samples and saved to database.'

    print '''
####################################################################################
# Step 6: count UMIs using selected clustering methods
####################################################################################
'''

    # lastly, group umis using umi-tools and import counts into db
    # each process spawns n_children processes, so ensure (n_samples * n_children) < available hardware threads
    n_children = 24
    count_umis_mp = [sample.count_umis(n_children, umi_grouping_methods) for sample in samples]
    # wait for all processes to finish before continuing
    wait(count_umis_mp)

    print '''
####################################################################################
# Step 7: cleanup and merge databases
####################################################################################
'''

    print 'Cleaning up databases...\n'

    # join tables to create summary tables
    clean_dbs = [sample.clean_db() for sample in samples]

    print 'Merging databases...\n'

    # merge databases containing different samples
    merge_dbs(samples)

    print '\nPipeline complete! All databases saved in %s.\n' % db_out_dir



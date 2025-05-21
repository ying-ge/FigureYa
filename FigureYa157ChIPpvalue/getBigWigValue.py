#!/usr/bin/env python3

# ------------------------------
# Import Modules
# ------------------------------
import os, sys
from optparse import OptionParser
import subprocess
from multiprocessing import cpu_count

# ------------------------------
# Sub process
# ------------------------------
signal_capture = '''#! /usr/bin/env python3
import os, sys
import subprocess
from distutils.spawn import find_executable

scripy_name = os.path.basename(sys.argv[0])
USAGE = f'{scripy_name} <name> <bedfile> <datapoints> <strand> <bigwigfiles+>'

datapoints = int(sys.argv[3])
strand_specific = int(sys.argv[4])
bwfhds = [x for x in sys.argv[5:]]

rfhds = [
    open(sys.argv[1] + '_siteprof' + str(i + 1), 'w')
    for i in range(len(bwfhds))
]
bed = {}
with open(sys.argv[2]) as bedFhd:
    for line in bedFhd:
        line = line.strip().split()
        chrom, start, end = line[0], int(line[1]), int(line[2])
        if strand_specific:
            strand = line[5]
        for i in range(len(bwfhds)):
            sp = subprocess.Popen([
                f'{find_executable("bigWigSummary")}'
                if find_executable("bigWigSummary") else '../bigWigSummary',
                f'{bwfhds[i]}', f'{chrom}', f'{max(0, start)}', f'{end}',
                f'{datapoints}'
            ],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  encoding='utf-8')
            out, err = sp.communicate()
            if out:
                output_values = [
                    '0' if value == 'n/a' else value
                    for value in out.strip().split()
                ]
            else:
                output_values = ['0' for datapoint in range(datapoints)]
            if strand_specific and strand == '-':
                output_values.reverse()
            rfhds[i].write('\\t'.join(output_values) + '\\n')
'''

coverageCapture = '''#! /usr/bin/env python3
import os, sys
import subprocess
from distutils.spawn import find_executable

scripy_name = os.path.basename(sys.argv[0])
USAGE = f'{scripy_name} <name> <bedfile> <datapoints> <strand> <bigwigfiles+>'

datapoints = int(sys.argv[3])
strand_specific = int(sys.argv[4])
bwfhds = [x for x in sys.argv[5:]]

rfhds = [
    open(sys.argv[1] + '_siteprof' + str(i + 1), 'w')
    for i in range(len(bwfhds))
]
bed = {}
with open(sys.argv[2]) as bedFhd:
    for line in bedFhd:
        line = line.strip().split()
        chrom, start, end = line[0], int(line[1]), int(line[2])
        if strand_specific:
            strand = line[5]
        for i in range(len(bwfhds)):
            sp = subprocess.Popen([
                f'{find_executable("bigWigSummary")}'
                if find_executable("bigWigSummary") else '../bigWigSummary',
                '-type=coverage', f'{bwfhds[i]}', f'{chrom}',
                f'{max(0, start)}', f'{end}', f'{datapoints}'
            ],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  encoding='utf-8')
            out, err = sp.communicate()
            if out:
                output_values = [
                    '0' if value == 'n/a' else value
                    for value in out.strip().split()
                ]
            else:
                output_values = ['0' for datapoint in range(datapoints)]
            if strand_specific and strand == '-':
                output_values.reverse()
            rfhds[i].write('\\t'.join(output_values) + '\\n')
'''


# ------------------------------
# Sub Funcitons
# ------------------------------
def prepare_optparser():
    '''
    Prepare optparser object. New options will be added in thisfunction first.
    '''
    scripy_name = os.path.basename(sys.argv[0])
    usage = f'USAGE: {scripy_name} <-b bedFile> <-w bigWigFiles+> [-n name] [-p process] [-s datapoints] [-m model]'
    description = 'bigWigSignalCapture -- Capture signal form bigWigFiles.'
    
    # option processor
    optparser = OptionParser(version=f'{scripy_name} 0.1',
                             description=description,
                             usage=usage,
                             add_help_option=False)
    
    # basic setting
    optparser.add_option('-h',
                         '--help',
                         action='help',
                         help='Show this help message and exit.')
    optparser.add_option('-n', '--name',dest='name', type='string',\
                         help='Name of this run. If not given, the time will be used.')
    optparser.add_option('-p', '--process',dest='process', type='int',\
                         help='Number of subprocess. If not given, number of computer processor will be used.',default=cpu_count())
    optparser.add_option('-s', '--datapoints',dest='datapoints', type='int',\
                         help='The capture region will broken into dataPoints equal parts. Default is 1.', default=1)
    optparser.add_option('-b', '--bedFile',dest='bedFile', type='string',\
                         help='BedFile for capture.')
    optparser.add_option('-w', '--bigWig',dest='bigWigFiles', type='string', action='append',\
                         help='bigWigFiles for capture.')
    optparser.add_option('-m', '--model',dest='model', type='string',\
                         help='Capture model: security or speed. Default is Security',default='security')
    optparser.add_option('--strand',dest='strand', action='store_true',\
                         help='If Strand specific')
    
    return optparser


def opt_validate(optparser):
    '''Validate options from a OptParser object.
    Ret: Validated options object.
    '''
    (options, args) = optparser.parse_args()
    
    # input bigwig and bed files must be given
    if not (options.bedFile and options.bigWigFiles):
        sys.stdout.write('Input bigwig and bed files must be given!\n')
        optparser.print_help()
        sys.exit(1)
    
    # check each input file
    for i in range(len(options.bigWigFiles)):
        if not os.path.isfile(options.bigWigFiles[i]):
            sys.stdout.write(
                f'Can not find the bigwig file: {options.bigWigFiles[i]}.\n')
            optparser.print_help()
            sys.exit(1)
    
    if not os.path.isfile(options.bedFile):
        sys.stdout.write(
            f'Can not find the bigwig file: {options.bedFile}.\n'.format)
        optparser.print_help()
        sys.exit(1)
    
    # check capture model
    if not (options.model == 'security' or options.model == 'speed'):
        sys.stdout.write('Invalid capture model! Must be security or speed.\n')
        optparser.print_help()
        sys.exit(1)
    
    # get name
    if not options.name:
        options.name = f'bigWigSignalCapture_{time.strftime("%Y.%b.%d.%H-%M-%S", time.localtime())}'
    
    return options


def load_bed(bedFile):
    with open(bedFile) as fhd:
        bed = [line for line in fhd]
    return bed


def prepare_capture(options):
    # check subprocess file
    with open('signalCapture.py', 'w') as fhd:
        fhd.write(signal_capture)
    with open('coverageCapture.py', 'w') as fhd:
        fhd.write(coverageCapture)
    os.system('chmod 755 signalCapture.py coverageCapture.py')
    
    bed = load_bed(options.bedFile)
    # creat temp dirctory
    CMD = f'mkdir {options.name} && cd {options.name}\n'
    for i in range(len(options.bigWigFiles)):
        CMD += f'ln -s ../{options.bigWigFiles[i]} .\n'
    os.system(CMD)
    # write temp bed files
    sub_sizes = int(len(bed) / options.process)
    sub_index = []
    for i in range(options.process - 1):
        sub_index.append((i * sub_sizes, i * sub_sizes + sub_sizes))
    sub_index.append((options.process * sub_sizes - sub_sizes, len(bed)))
    for i in range(options.process):
        with open(f'{options.name}/temp_{options.name}_{i+1}.bed', 'w') as fhd:
            fhd.write(''.join(bed[sub_index[i][0]:sub_index[i][1]]))


def capture_signal(options):
    strand_parameter = 1 if options.strand else 0
    bw_files_parameter = ' '.join(
        os.path.basename(filename) for filename in options.bigWigFiles)
    # Caputre signal
    processes = []
    for i in range(options.process):
        CMD = f'cd {options.name} && python3 ../signalCapture.py temp_signal_{options.name}_{i+1} temp_{options.name}_{i+1}.bed {options.datapoints} {strand_parameter} {bw_files_parameter}'
        # print('Sub Command: ', CMD)
        processes.append(subprocess.Popen(CMD, shell=True))
    stats = []
    for i in range(options.process):
        stats.append(os.waitpid(processes[i].pid, 0))
    CMD = ''
    for i in range(len(options.bigWigFiles)):
        CMD += 'cat '
        for j in range(options.process):
            CMD += f'{options.name}/temp_signal_{options.name}_{j+1}_siteprof{i+1} '
        CMD += f'| gzip - > signal_{options.name}_siteprof{i+1}.gz && \\\n'
    CMD = CMD[:-5]
    with open('tem.sh', 'w') as fhd:
        fhd.write(CMD)
    # print('Cat Command: ', CMD)
    p = subprocess.Popen('bash tem.sh', shell=True)
    try:
        os.waitpid(p.pid, 0)
    except OSError:
        pass
    # Coverage capture
    if options.model == 'security':
        processes = []
        for i in range(options.process):
            CMD = f'cd {options.name} && python3 ../coverageCapture.py temp_coverage_{options.name}_{i+1} temp_{options.name}_{i+1}.bed {options.datapoints} {strand_parameter} {bw_files_parameter}'
            # print('Sub Command: ', CMD)
            processes.append(subprocess.Popen(CMD, shell=True))
        stats = []
        for i in range(options.process):
            stats.append(os.waitpid(processes[i].pid, 0))
        CMD = ''
        for i in range(len(options.bigWigFiles)):
            CMD += 'cat '
            for j in range(options.process):
                CMD += f'{options.name}/temp_coverage_{options.name}_{j+1}_siteprof{i+1} '
            CMD += f'| gzip - > coverage_{options.name}_siteprof{i+1}.gz && \\\n'
        CMD = CMD[:-5]
        with open('tem.sh', 'w') as fhd:
            fhd.write(CMD)
        # print('Cat Command: ', CMD)
        p = subprocess.Popen('bash tem.sh', shell=True)
        try:
            os.waitpid(p.pid, 0)
        except OSError:
            pass
    p = subprocess.Popen(
        f'rm -r {options.name} tem.sh signalCapture.py coverageCapture.py',
        shell=True)
    try:
        os.waitpid(p.pid, 0)
    except OSError:
        pass


# ------------------------------
# Main Funcitons
# ------------------------------
def main():
    
    # read the options and validate them
    options = opt_validate(prepare_optparser())
    
    # prepare capture
    prepare_capture(options)
    
    # capture signal
    capture_signal(options)


# ------------------------------
# Program Running
# ------------------------------

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        Info('User interrupts me! ;-) See you!')
        sys.exit(0)
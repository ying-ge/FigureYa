#!/usr/bin/env python2

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
signal_capture = """#! /usr/bin/env python
import os, sys
import subprocess

USAGE = '{} <name> <bedfile> <datapoints> <bigwigfiles>'.format(os.path.basename(sys.argv[0]))

datapoints = int(sys.argv[3])
bwfhds = [x for x in sys.argv[4:]]

rfhds = [open(sys.argv[1]+'_siteprof'+str(i+1), 'w') for i in range(len(bwfhds))]
bed = {}
with open(sys.argv[2]) as bedFhd:
    for line in bedFhd:
        line = line.strip().split()
        chrom, start, end = line[0], int(line[1]), int(line[2])
        for i in range(len(bwfhds)):
            sp = subprocess.Popen(['../bigWigSummary', '{}'.format(bwfhds[i]), '{}'.format(chrom), '{}'.format(max(0, start)), '{}'.format(end), '{}'.format(datapoints)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = sp.communicate()
            if out:
                rfhds[i].write('\\t'.join(['0' if value == 'n/a' else value for value in out.strip().split()]) + '\\n')
            else:
                rfhds[i].write('\\t'.join(['0' for i in range(datapoints)]) + '\\n')
"""

coverageCapture = """#! /usr/bin/env python
import os, sys
import subprocess

USAGE = '{} <name> <bedfile> <datapoints> <bigwigfiles>'.format(os.path.basename(sys.argv[0]))

datapoints = int(sys.argv[3])
bwfhds = [x for x in sys.argv[4:]]

rfhds = [open(sys.argv[1]+'_siteprof'+str(i+1), 'w') for i in range(len(bwfhds))]
bed = {}
with open(sys.argv[2]) as bedFhd:
    for line in bedFhd:
        line = line.strip().split()
        chrom, start, end = line[0], int(line[1]), int(line[2])
        for i in range(len(bwfhds)):
            sp = subprocess.Popen(['../bigWigSummary', "-type=coverage", '{}'.format(bwfhds[i]), '{}'.format(chrom), '{}'.format(max(0, start)), '{}'.format(end), '{}'.format(datapoints)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = sp.communicate()
            if out:
                rfhds[i].write('\\t'.join(['0' if value == 'n/a' else value for value in out.strip().split()]) + '\\n')
            else:
                rfhds[i].write('\\t'.join(['0' for i in range(datapoints)]) + '\\n')
"""

# ------------------------------
# Sub Funcitons
# ------------------------------
def prepare_optparser():
    
    """
    Prepare optparser object. New options will be added in thisfunction first.
    """
    usage = "USAGE: {} <-b bedFile> <-w bigWigFiles+> [-n name] [-p process] [-s datapoints] [-m model]".format(os.path.basename(sys.argv[0]))
    description = "bigWigSignalCapture -- Capture signal form bigWigFiles."
    
    # option processor
    optparser = OptionParser(version="{} 0.1".format(os.path.basename(sys.argv[0])), description=description, usage=usage, add_help_option=False)
    
    # basic setting
    optparser.add_option("-h","--help", action="help", help="Show this help message and exit.")
    optparser.add_option("-n", "--name",dest="name", type="string",\
                         help="Name of this run. If not given, the time will be used.")
    optparser.add_option("-p", "--process",dest="process", type="int",\
                         help="Number of subprocess. If not given, number of computer processor will be used.",default=cpu_count())
    optparser.add_option("-s", "--datapoints",dest="datapoints", type="int",\
                         help="The capture region will broken into dataPoints equal parts. Default is 1.", default=1)
    optparser.add_option("-b", "--bedFile",dest="bedFile", type="string",\
                         help="BedFile for capture.")
    optparser.add_option("-w", "--bigWig",dest="bigWigFiles", type="string", action="append",\
                         help="bigWigFiles for capture.")
    optparser.add_option("-m", "--model",dest="model", type="string",\
                         help="Capture model: security or speed. Default is Security",default="security")
    
    return optparser

def opt_validate(optparser):
    
    """Validate options from a OptParser object.
    Ret: Validated options object.
    """
    (options,args) = optparser.parse_args()
    
    # input bigwig and bed files must be given 
    if not (options.bedFile and options.bigWigFiles):
        sys.stdout.write("Input bigwig and bed files must be given!\n")
        optparser.print_help()
        sys.exit(1)
    
    # check each input file
    for i in range(len(options.bigWigFiles)):
        if not os.path.isfile(options.bigWigFiles[i]):
            sys.stdout.write("Can not find the bigwig file: {}.\n".format(options.bigWigFiles[i]))
            optparser.print_help()
            sys.exit(1)
    
    if not os.path.isfile(options.bedFile):
        sys.stdout.write("Can not find the bigwig file: {}.\n".format(options.bedFile))
        optparser.print_help()
        sys.exit(1)
    
    # check capture model
    if not (options.model == "security" or options.model == "speed"):
        sys.stdout.write("Invalid capture model! Must be security or speed.\n")
        optparser.print_help()
        sys.exit(1)
    
    # get name
    if not options.name:
        options.name="bigWigSignalCapture_%s" %(time.strftime("%Y.%b.%d.%H-%M-%S", time.localtime()))
    
    return options

def load_bed(bedFile):
    with open(bedFile) as fhd:
        bed = [line for line in fhd]
    return bed

def prepare_capture(options):
    # check subprocess file
    if not os.path.isfile("signalCapture.py"):
        with open("signalCapture.py", "w") as fhd:
            fhd.write(signal_capture)
    if not os.path.isfile("coverageCapture.py"):
        with open("coverageCapture.py", "w") as fhd:
            fhd.write(coverageCapture)
    os.system("chmod 755 signalCapture.py coverageCapture.py")
    
    bed = load_bed(options.bedFile)
    # creat temp dirctory
    CMD = "mkdir {} && cd {}\n".format(options.name, options.name)
    for i in range(len(options.bigWigFiles)):
        CMD += "ln -s ../{} .\n".format(options.bigWigFiles[i])
    os.system(CMD)
    # write temp bed files
    sub_sizes = len(bed) / options.process
    sub_index = []
    for i in range(options.process-1):
        sub_index.append((i*sub_sizes, i*sub_sizes+sub_sizes))
    sub_index.append((options.process*sub_sizes-sub_sizes,len(bed)))
    for i in range(options.process):
        with open("{}/temp_{}_{}.bed".format(options.name, options.name, i+1), "w") as fhd:
            fhd.write("".join(bed[sub_index[i][0]:sub_index[i][1]]))

def capture_signal(options):
    # Caputre signal
    processes = []
    for i in range(options.process):
        CMD = "cd {} && python ../signalCapture.py temp_signal_{}_{} temp_{}_{}.bed {} {}".format(options.name, options.name, i+1, options.name, i+1, options.datapoints, " ".join(options.bigWigFiles))
        # print("Sub Command: ", CMD)
        processes.append(subprocess.Popen(CMD, shell=True))
    stats = []
    for i in range(options.process):
        stats.append(os.waitpid(processes[i].pid,0))
    CMD = ''
    for i in range(len(options.bigWigFiles)):
        CMD += "cat "
        for j in range(options.process):
            CMD += "{}/temp_signal_{}_{}_siteprof{} ".format(options.name, options.name, j+1, i+1)
        CMD += "| gzip - > signal_{}_siteprof{}.gz && \\\n".format(options.name, i+1)
    CMD = CMD[:-5]
    with open('tem.sh', 'w') as fhd:
        fhd.write(CMD)
    # print("Cat Command: ", CMD)
    p = subprocess.Popen('bash tem.sh', shell = True)
    try:
        os.waitpid(p.pid, 0)
    except OSError:
        pass
    # Coverage capture
    if options.model == "security":
        processes = []
        for i in range(options.process):
            CMD = "cd {} && python ../coverageCapture.py temp_coverage_{}_{} temp_{}_{}.bed {} {}".format(options.name, options.name, i+1, options.name, i+1, options.datapoints, " ".join(options.bigWigFiles))
            # print("Sub Command: ", CMD)
            processes.append(subprocess.Popen(CMD, shell=True))
        stats = []
        for i in range(options.process):
            stats.append(os.waitpid(processes[i].pid,0))
        CMD = ''
        for i in range(len(options.bigWigFiles)):
            CMD += "cat "
            for j in range(options.process):
                CMD += "{}/temp_coverage_{}_{}_siteprof{} ".format(options.name, options.name, j+1, i+1)
            CMD += "| gzip - > coverage_{}_siteprof{}.gz && \\\n".format(options.name, i+1)
        CMD = CMD[:-5]
        with open('tem.sh', 'w') as fhd:
            fhd.write(CMD)
        # print("Cat Command: ", CMD)
        p = subprocess.Popen('bash tem.sh', shell = True)
        try:
            os.waitpid(p.pid, 0)
        except OSError:
            pass
    p = subprocess.Popen('rm -r {} tem.sh signalCapture.py coverageCapture.py'.format(options.name), shell = True)
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
        Info("User interrupts me! ;-) See you!")
        sys.exit(0)

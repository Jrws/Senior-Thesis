import os
import argparse
import fileinput
import numpy as np
from neo import io

# To run:
# > python asciiToNix.py "SSR files/Nov1.ASC"

# Parse file name argument by running script from command line
parser = argparse.ArgumentParser("Convert large ASCII SSR files to individual NIX files")
parser.add_argument('file', help="Specifies ASCII file to be processed")
args = parser.parse_args()      # only arg is path to file
if not os.path.isfile(args.file):
    raise SystemExit(f'File {args.file} not found!')

# create folders for ascii, binary, and nix subfiles
asc_dir = '.'.join(args.file.split('.')[:-1])   # name of .ASC file without file extension as dirname
for folder in ['','ascii','bin','nix','plt']:
    if not os.path.isdir(os.path.join(asc_dir,folder)):
        os.mkdir(os.path.join(asc_dir,folder))

# Process file line by line
is_SSR, f, signalNums, sampleRates = False, None, [], []
for line in fileinput.input(args.file):
    if line[0] == ';':                # if line is commented out; contains metadata
        if line.endswith('-3\n'):     # ... SigX-3 (denotes start of SSR daa)
            if f != None and not f.closed:
                f.close()
            sigNum = line.split(' ')[-1][3:-3]          # extract signal number from this metadata line
            signalNums.append(sigNum)
            currentASC = os.path.join(asc_dir, 'ascii', sigNum + '.asc')   # file to write ascii data is signal number + .ASC, placed in ascii folder
            f = open(currentASC, 'w')
            f.write(line)
            is_SSR = True
            continue
        elif is_SSR and line.endswith(('-4\n','-D\n')):   # ... SigX-4 (EAG), SigX-D (Digital): start of irrelevant data
            f.close()
            is_SSR = False
    if is_SSR:
        f.write(line)
if f != None and not f.closed:      # closing any open file at end
    f.close()

for sigNum in signalNums:       # file conversions
    trace_ascii = os.path.join(asc_dir, 'ascii', sigNum + '.asc')
    trace_binary = os.path.join(asc_dir, 'bin', sigNum + '.bin')
    trace_nix = os.path.join(asc_dir, 'nix', sigNum + '.nix')

    sampleRate = None
    with open(trace_ascii, 'r') as a, open(trace_binary, 'wb') as b:    # convery ASCII to binary file
        np.genfromtxt(a, skip_header=4)[:,1].tofile(b)
        a.seek(0)
        for line in a:
            if line.startswith('; Sample rate'):
                sampleRate = float(line.split(' ')[-1])
                break

    # convert binary to neo NIX file
    r = io.RawBinarySignalIO(filename=trace_binary, dtype=float, nb_channel=1, sampling_rate=sampleRate)
    writer = io.NixIO(trace_nix, mode="ow")
    writer.write(r.read_block())
    writer.close()

with open(os.path.join(asc_dir, 'signals.txt'),'w') as f:
    f.write('\n'.join(signalNums))

open(os.path.join(asc_dir, 'signalsToTest.txt'),'w').close()

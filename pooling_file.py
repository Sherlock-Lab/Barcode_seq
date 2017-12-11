"""create a single file of barcode counts per condition across all timepoints"""
# this script will create a single file of barcode counts per
# condition across all timepoints

# the output is a file named "Evol_generation_condition.txt" with the
# BC \t count form the early to the late generation comma separated.

# to make an all csv file only you'll have to apply the
# changing_file.py but because only a subset of all those files are
# intersting you'll need to put them in a new directory evolution_row

import sys
import os

import os.path

# a single command line argument is supplied, which is the name of the
# condition, which allows it to know in which directories it needs to
# look to find the data for each timepoint

CONDITION = sys.argv[1]

COUNTS = {}
GENERATIONS = []

# we're going to loop over all the files, in generation order, and
# store the data in a two dimensional dictionary, then output the data
# by iterating over each identifier

LAST_GENERATION = 248
TIME_INTERVAL = 8

for generation in range(0, LAST_GENERATION, TIME_INTERVAL):

    Timepoint_file = 'T%s_%s/cluster_final.txt' % (generation, CONDITION)

    if os.path.exists(Timepoint_file):

        print Timepoint_file
        print generation

        GENERATIONS.append(generation)

        Timepoint_file_handle = open(Timepoint_file, 'r')

        # discard the header line

        next(Timepoint_file_handle)

        for line in Timepoint_file_handle:

            (clust_ID, BC, freq) = line.rstrip('\n').split('\t')

            if BC not in COUNTS:
                COUNTS[BC] = {}

            COUNTS[BC][generation] = freq


# now generate the output file

EVOLUTION_FILE = open('pool_file/Evol_%s_%s.txt' % (GENERATIONS[-1], CONDITION), 'w')

# print out a header

EVOLUTION_FILE.write("BARCODE,")
EVOLUTION_FILE.write(",".join(map(str, GENERATIONS)) + "\n")

# now we simply go over our dictionary, and print out its contents

for BC in COUNTS:

    EVOLUTION_FILE.write(BC)

    for generation in GENERATIONS:

        EVOLUTION_FILE.write(",")

        if generation in COUNTS[BC]:

            EVOLUTION_FILE.write(str(COUNTS[BC][generation]))

        else:

            EVOLUTION_FILE.write("0")

    EVOLUTION_FILE.write("\n")

EVOLUTION_FILE.close()

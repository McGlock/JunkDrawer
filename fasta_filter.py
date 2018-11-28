#!/usr/bin/env python
import sys


# No bioinformatic software would be complete without a contribution from Heng Li.
# Adapted from his readfq generator
def generate_fasta(fasta_handler):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:
            for line in fasta_handler:  # search for the start of the next record
                if line[0] == '>':  # fasta header line
                    last = line[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:], [], None
        for line in fasta_handler:  # read the sequence
            if line[0] == '>':
                last = line[:-1]
                break
            seqs.append(line[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs)  # yield a fasta record
            if not last:
                break
        else:
            seq, seqs = ''.join(seqs), []
            for line in fasta_handler:  # read the quality
                seqs.append(line[:-1])
            if last:  # reach EOF before reading enough quality
                yield name, seq  # yield a fasta record instead
                break


def read_fasta_to_dict(fasta_file):
    fasta_dict = dict()
    try:
        fasta_handler = open(fasta_file, 'r')
    except IOError:
        logging.error("Unable to open " + fasta_file + " for reading!\n")
        sys.exit(5)
    for record in generate_fasta(fasta_handler):
        name, sequence = record
        fasta_dict[name] = sequence.upper()
    return fasta_dict


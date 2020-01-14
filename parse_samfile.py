mate_pairs = {}
QNAME_list = set(QNAME_list) # make the checking a tiny bit faster by making the list a set
with open(pathToSam, 'rb') as samFile:
    for line in samFile:
        QNAME = line.split('\t')[0]
        if QNAME in QNAME_list:
            try:
                mate_pairs[QNAME].append(line)
            except KeyError:
                mate_pairs[QNAME] = [line]
with open(parsed_SAM_file_path, "wb") as outfile:
    for read,mate in mate_pairs.values():
        outfile.write(read + '\n' + mate + '\n')

import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO


class SequenceAnalyzer:

    # Constructor
    def __init__(self, program, database, format_type, valid_query_coverage):
        self.program = program
        self.record = ""
        self.database = database
        self.format_type = format_type
        self.valid_query_coverage = valid_query_coverage

    # Method to call blast search, we save result to a file in case we can check what happened
    def quickBlastSearch(self, query, with_file = "qblast_search.xml"):
        HITLIST = 1000
        blast_output = NCBIWWW.qblast(self.program,
            self.database,
            self.record.seq,
            entrez_query= query,
            format_type=self.format_type,
            hitlist_size=HITLIST)
        save_file = open(with_file, "w")
        save_file.write(blast_output.read())
        return

    # Method to load main sequence
    def loadSequenceFromPath(self, path, format="fasta"):
        self.record = SeqIO.read(open(path), format=format)

    # Method to return result of qblast
    # we call quickBlastSearch and specify a file to save xml result to
    def qBlastResult(self, query, with_file = "qblast_search.xml"):
        self.quickBlastSearch(query, with_file)
        result_handle = open(with_file)
        result = NCBIXML.read(result_handle)
        return result;

    # Returns only those records with query coverage over valid given query coverage
    def qBlastResultWithFilter(self, query, with_file):
        brecord = self.qBlastResult(query, with_file)
        not_valid_aligments = []
        query_letters = brecord.query_letters

        # Calculate each alignment query coverage and add invalid records to not_valid_alignments
        for alignment in brecord.alignments:
            match_len_sum = 0
            for hsp in alignment.hsps:
                match_len_sum += hsp.align_length
            val = (match_len_sum*100)/query_letters
            if val < self.valid_query_coverage:
                not_valid_aligments.append(alignment)

        # Remove invalid alignments
        for alignment in not_valid_aligments:
            brecord.alignments.remove(alignment)

        # Return only valid records
        return brecord

    # Method to save blast record in FASTA format
    def saveBlastRecordAsFasta(self, brecord, save_to = "blast_record.fasta"):
        fasta_list = []
        for alignment in brecord.alignments:
            for hsp in alignment.hsps:
                fasta_list.append('>')
                fasta_list.append(alignment.title)
                fasta_list.append('\n')
                fasta_list.append(hsp.sbjct)
                fasta_list.append('\n')
                fasta_list.append('\n')
        with open(save_to,"w") as temp_file:
            temp_file.writelines(fasta_list)

    # Align sequences using mafft algorithm
    # MAFFT must be installed, and for os environmental path must be specified to folder containing mafft
    def align_mafft(self, with_file = "blast_record.fasta"):
        os.system('mafft --quiet ' + with_file + ' > mafft.fasta')
        return

    # Removes identical sequences from file
    def remove_identical(self, input_file, output_file):
        os.system('cd_hit -i ' + input_file + ' -o ' + output_file + ' -c 1')
        return

    # Method to calculate mistakes in a part starting from some index
    def mistakes_in_part(self, main_sequence, probe, start_index, part_length):
        mistakes_count = 0;
        i = 0
        # count mistakes in main_sequence compared to that sequence that we took in a range
        for j in range (start_index, start_index + part_length):
            if (probe[i] != main_sequence[j]):
                mistakes_count += 1
                i += 1
        return mistakes_count

    # Simply returns all sequences from the file in fasta format
    def return_sequences(self, with_file = "mafft.fasta"):
        #Reading source file to get sequences
        all_seqs = []
        for record in SeqIO.parse(open(with_file), "fasta"):
            all_seqs.append(record.seq)
        return all_seqs


# Function which appends one file to another
def append_file_to_file(f, file):
    r = open(file, "r")
    f.write(r.read())
    r.close()
    return

# counts and return all letters in a sequence
def countLetter(seq):
    counter = 0
    for letter in seq:
        if letter != '-':
            counter += 1
    return counter


# Main function
def main():
    program = "blastn"
    main_record = "main.fasta"
    database = "nr"
    format_type = "XML"
    valid_query_coverage = 60

    LIST_DANGEROUS = [16, 18, 31, 33, 35, 51, 52]
    LIST_HARMLESS = [6, 11, 40, 42, 43, 44, 57, 81]

    sequence_analyzer = SequenceAnalyzer(program, database, format_type, valid_query_coverage)

    sequence_analyzer.loadSequenceFromPath(main_record)

    filename_merged = "merged_sequences.fasta"
    
    f = open(filename_merged, "w")

    # Firstly working with dangerous types and only then with harmless
    # 1. Perform blastn select and filter it by valid query coverage
    # 2. Convert from XML to fasta format
    # 3. Remove identical lines using external program cd-hit
    # 4. Put everything into one file. REMEMBER: firstly dangerous and only then harmless types
    for number in LIST_DANGEROUS:
        entrez_query = '"papillomavirus"[Organism] AND ( *"type ' + str(number) + '"*[title] AND *human*[title])'
        filename = 'dangerous_type_' + str(number)
        blast_record = sequence_analyzer.qBlastResultWithFilter(entrez_query, filename + '.xml')
        sequence_analyzer.saveBlastRecordAsFasta(blast_record, filename + '.fasta')
        sequence_analyzer.remove_identical(filename + '.fasta', filename + '_after.fasta')
        append_file_to_file(f, filename + '_after.fasta')

    for number in LIST_HARMLESS:
        entrez_query = '"papillomavirus"[Organism] AND ( *"type ' + str(number) + '"*[title] AND *human*[title])'
        filename = 'harmless_type_' + str(number)
        blast_record = sequence_analyzer.qBlastResultWithFilter(entrez_query, filename + '.xml')
        sequence_analyzer.saveBlastRecordAsFasta(blast_record, filename + '.fasta')
        sequence_analyzer.remove_identical(filename + '.fasta', filename + '_after.fasta')
        append_file_to_file(f, filename + '_after.fasta')

    f.close()

    # align sequences with external program called - mafft
    sequence_analyzer.align_mafft(filename_merged)

    # calculate how many dangerous sequences we have
    counter = 0
    for number in LIST_DANGEROUS:
        filename = 'dangerous_type_' + str(number) + '_after.fasta'
        for record in SeqIO.parse(open(filename), "fasta"):
            counter += 1

    # calculate how many harmless sequences we have
    counter_all = 0
    filenameMafft = 'mafft.fasta'
    for record in SeqIO.parse(open(filenameMafft), "fasta"):
        counter_all += 1

    # define variable to store all sequences and initiate it
    all_sequences = []
    all_sequences = sequence_analyzer.return_sequences()

    probe_system = []
    length = all_sequences[1].__len__()

    # find all possible probes. Probes must be 25-35BP. And u must take from the region in range of 60BP
    # probes must be identified from dangerous types

    # for all dangerous types
    for i in range(0, counter):
        # define current sequences, and region approbation
        curr_sequence = all_sequences[i]
        startSearch = 0
        endSearch = length - 60
        # for all signs in a sequence
        for j in range(startSearch, endSearch):
            # search in a region
            for region in range(j, j + 60):
                # simple check not to go out of bounds
                if (j + 60) > length:
                    break
                # approbations for a probe
                minProbeL = 25
                maxProbeL = 35
                x = minProbeL
                y = 60
                # for 25-60. we know we can't make probe with lower that 25BP
                for z in range(x, y):
                    # get substring
                    substring = curr_sequence[j:j+z]
                    # compare with probe approbations and append it to whole probe system
                    if (countLetter(substring) >= minProbeL) and (countLetter(substring) <= maxProbeL):
                        probe_system.append(substring)

    # list of probes what are not supposed to be in probe system
    bad_probes = []
    # filter with harmless types.
    # Probe must have minimum 3 differences between harmless types

    # for all harmless sequences
    for i in range(counter, counter_all):
        # take that sequence
        curr_sequence = all_sequences[i]
        # for all probes in the system
        for probe in probe_system:
            # for number in sequence length(all sequences same length, because of mafft alignment
            for j in range (0, length):
                # simple check not to go out of bounds
                if (probe.__len__() + j) >= length:
                    break
                # calculate mistakes and identify bad probes
                if sequence_analyzer.mistakes_in_part(curr_sequence, probe, j, probe.__len__) < 3:
                    bad_probes.append(probe)
                    break

    # remove bad probes from original probe list
    for probe in bad_probes:
        probe_system.remove(probe)

    # Last filter. We must form minimum length probe system, so for each sequence in dangerous types
    # we find exactly one probe to match that sequence with no more than two mistakes

    final_probe_list = []

    # for all dangerous types
    for i in range(0, counter): # for all dangerous sequences
        # identify current sequence
        curr_sequence = all_sequences[i]
        # for all probes in a filtered system
        for probe in probe_system:
            # for all symbols in sequence
            for j in range (0, length):
                # check if not out of bounds
                if (probe.__len__() + j) >= length:
                    break
                # calculate mistakes and append correct probes
                if sequence_analyzer.mistakes_in_part(curr_sequence, probe, j, probe.__len__) <= 2:
                    final_probe_list.append(probe)
                    break

    # Prints only filtered minimum probe system
    for probe in final_probe_list:
        print probe
# Run main
main()
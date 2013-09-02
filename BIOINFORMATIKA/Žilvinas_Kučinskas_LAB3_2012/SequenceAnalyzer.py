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
        blast_output = NCBIWWW.qblast(self.program,
            self.database,
            self.record.seq,
            entrez_query= query,
            format_type=self.format_type)
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

# Function which appends one file to another
def append_file_to_file(f, file):
    r = open(file, "r")
    f.write(r.read())
    r.close()
    return

# Main function
def main():
    program = "blastn"
    main_record = "main.fasta"
    database = "nr"
    format_type = "XML"
    valid_query_coverage = 60

    # define dangerous type xml files
    filename_dangerous_type_16_xml = 'dangerous_type_16.xml'
    filename_dangerous_type_18_xml = 'dangerous_type_18.xml'
    filename_dangerous_type_31_xml = 'dangerous_type_31.xml'
    filename_dangerous_type_33_xml = 'dangerous_type_33.xml'
    filename_dangerous_type_35_xml = 'dangerous_type_35.xml'
    filename_dangerous_type_51_xml = 'dangerous_type_51.xml'
    filename_dangerous_type_52_xml = 'dangerous_type_52.xml'

    # define dangerous type fasta files
    filename_dangerous_type_16_fasta = 'dangerous_type_16.fasta'
    filename_dangerous_type_18_fasta = 'dangerous_type_18.fasta'
    filename_dangerous_type_31_fasta = 'dangerous_type_31.fasta'
    filename_dangerous_type_33_fasta = 'dangerous_type_33.fasta'
    filename_dangerous_type_35_fasta = 'dangerous_type_35.fasta'
    filename_dangerous_type_51_fasta = 'dangerous_type_51.fasta'
    filename_dangerous_type_52_fasta = 'dangerous_type_52.fasta'

    # define dangerous type entrez queries
    entrez_query_dangerous_type_16 = '"papillomavirus"[Organism] AND ( *"type 16"*[title] AND *human*[title])'
    entrez_query_dangerous_type_18 = '"papillomavirus"[Organism] AND ( *"type 18"*[title] AND *human*[title])'
    entrez_query_dangerous_type_31 = '"papillomavirus"[Organism] AND ( *"type 31"*[title] AND *human*[title])'
    entrez_query_dangerous_type_33 = '"papillomavirus"[Organism] AND ( *"type 33"*[title] AND *human*[title])'
    entrez_query_dangerous_type_35 = '"papillomavirus"[Organism] AND ( *"type 35"*[title] AND *human*[title])'
    entrez_query_dangerous_type_51 = '"papillomavirus"[Organism] AND ( *"type 51"*[title] AND *human*[title])'
    entrez_query_dangerous_type_52 = '"papillomavirus"[Organism] AND ( *"type 52"*[title] AND *human*[title])'

    # define harmless type xml files
    filename_harmless_type_6_xml = 'harmless_type_6.xml'
    filename_harmless_type_11_xml = 'harmless_type_11.xml'
    filename_harmless_type_40_xml = 'harmless_type_40.xml'
    filename_harmless_type_42_xml = 'harmless_type_42.xml'
    filename_harmless_type_43_xml = 'harmless_type_43.xml'
    filename_harmless_type_44_xml = 'harmless_type_44.xml'
    filename_harmless_type_57_xml = 'harmless_type_57.xml'
    filename_harmless_type_81_xml = 'harmless_type_81.xml'

    # define harmless type fasta files
    filename_harmless_type_6_fasta = 'harmless_type_6.fasta'
    filename_harmless_type_11_fasta = 'harmless_type_11.fasta'
    filename_harmless_type_40_fasta = 'harmless_type_40.fasta'
    filename_harmless_type_42_fasta = 'harmless_type_42.fasta'
    filename_harmless_type_43_fasta = 'harmless_type_43.fasta'
    filename_harmless_type_44_fasta = 'harmless_type_44.fasta'
    filename_harmless_type_57_fasta = 'harmless_type_57.fasta'
    filename_harmless_type_81_fasta = 'harmless_type_81.fasta'

    # define harmless type entrez queries
    entrez_query_harmless_type_6 = '"papillomavirus"[Organism] AND ( *"type 6"*[title] AND *human*[title])'
    entrez_query_harmless_type_11 = '"papillomavirus"[Organism] AND ( *"type 11"*[title] AND *human*[title])'
    entrez_query_harmless_type_40 = '"papillomavirus"[Organism] AND ( *"type 40"*[title] AND *human*[title])'
    entrez_query_harmless_type_42 = '"papillomavirus"[Organism] AND ( *"type 42"*[title] AND *human*[title])'
    entrez_query_harmless_type_43 = '"papillomavirus"[Organism] AND ( *"type 43"*[title] AND *human*[title])'
    entrez_query_harmless_type_44 = '"papillomavirus"[Organism] AND ( *"type 44"*[title] AND *human*[title])'
    entrez_query_harmless_type_57 = '"papillomavirus"[Organism] AND ( *"type 57"*[title] AND *human*[title])'
    entrez_query_harmless_type_81 = '"papillomavirus"[Organism] AND ( *"type 81"*[title] AND *human*[title])'


    sequence_analyzer = SequenceAnalyzer(program, database, format_type, valid_query_coverage)
    sequence_analyzer.loadSequenceFromPath(main_record)

    # form xml files with specifies query coverage
    blast_record_dangerous_type16 = sequence_analyzer.qBlastResultWithFilter(entrez_query_dangerous_type_16,
                                                                   filename_dangerous_type_16_xml)
    blast_record_dangerous_type18 = sequence_analyzer.qBlastResultWithFilter(entrez_query_dangerous_type_18,
                                                                             filename_dangerous_type_18_xml)
    blast_record_dangerous_type31 = sequence_analyzer.qBlastResultWithFilter(entrez_query_dangerous_type_31,
                                                                             filename_dangerous_type_31_xml)
    blast_record_dangerous_type33 = sequence_analyzer.qBlastResultWithFilter(entrez_query_dangerous_type_33,
                                                                             filename_dangerous_type_33_xml)
    blast_record_dangerous_type35 = sequence_analyzer.qBlastResultWithFilter(entrez_query_dangerous_type_35,
                                                                             filename_dangerous_type_35_xml)
    blast_record_dangerous_type51 = sequence_analyzer.qBlastResultWithFilter(entrez_query_dangerous_type_51,
                                                                             filename_dangerous_type_51_xml)
    blast_record_dangerous_type52 = sequence_analyzer.qBlastResultWithFilter(entrez_query_dangerous_type_52,
                                                                             filename_dangerous_type_52_xml)

    blast_record_harmless_type6  = sequence_analyzer.qBlastResultWithFilter(entrez_query_harmless_type_6,
                                                                            filename_harmless_type_6_xml)
    blast_record_harmless_type11 = sequence_analyzer.qBlastResultWithFilter(entrez_query_harmless_type_11,
                                                                           filename_harmless_type_11_xml)
    blast_record_harmless_type40 = sequence_analyzer.qBlastResultWithFilter(entrez_query_harmless_type_40,
                                                                            filename_harmless_type_40_xml)
    blast_record_harmless_type42 = sequence_analyzer.qBlastResultWithFilter(entrez_query_harmless_type_42,
                                                                            filename_harmless_type_42_xml)
    blast_record_harmless_type43 = sequence_analyzer.qBlastResultWithFilter(entrez_query_harmless_type_43,
                                                                            filename_harmless_type_43_xml)
    blast_record_harmless_type44 = sequence_analyzer.qBlastResultWithFilter(entrez_query_harmless_type_44,
                                                                            filename_harmless_type_44_xml)
    blast_record_harmless_type57 = sequence_analyzer.qBlastResultWithFilter(entrez_query_harmless_type_57,
                                                                            filename_harmless_type_57_xml)
    blast_record_harmless_type81 = sequence_analyzer.qBlastResultWithFilter(entrez_query_harmless_type_81,
                                                                            filename_harmless_type_81_xml)

    # Converting from xml to fasta
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_dangerous_type16, filename_dangerous_type_16_fasta)
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_dangerous_type18, filename_dangerous_type_18_fasta)
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_dangerous_type31, filename_dangerous_type_31_fasta)
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_dangerous_type33, filename_dangerous_type_33_fasta)
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_dangerous_type35, filename_dangerous_type_35_fasta)
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_dangerous_type51, filename_dangerous_type_51_fasta)
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_dangerous_type52, filename_dangerous_type_52_fasta)

    sequence_analyzer.saveBlastRecordAsFasta(blast_record_harmless_type6, filename_harmless_type_6_fasta)
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_harmless_type11, filename_harmless_type_11_fasta)
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_harmless_type40, filename_harmless_type_40_fasta)
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_harmless_type42, filename_harmless_type_42_fasta)
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_harmless_type43, filename_harmless_type_43_fasta)
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_harmless_type44, filename_harmless_type_44_fasta)
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_harmless_type57, filename_harmless_type_57_fasta)
    sequence_analyzer.saveBlastRecordAsFasta(blast_record_harmless_type81, filename_harmless_type_81_fasta)

    # define dangerous type files after cd-hit processing
    filename_dangerous_type_16_fasta_after = 'dangerous_type_16_after.fasta'
    filename_dangerous_type_18_fasta_after = 'dangerous_type_18_after.fasta'
    filename_dangerous_type_31_fasta_after = 'dangerous_type_31_after.fasta'
    filename_dangerous_type_33_fasta_after = 'dangerous_type_33_after.fasta'
    filename_dangerous_type_35_fasta_after = 'dangerous_type_35_after.fasta'
    filename_dangerous_type_51_fasta_after = 'dangerous_type_51_after.fasta'
    filename_dangerous_type_52_fasta_after = 'dangerous_type_52_after.fasta'

    # remove identical sequences with cd-hit
    sequence_analyzer.remove_identical(filename_dangerous_type_16_fasta, filename_dangerous_type_16_fasta_after)
    sequence_analyzer.remove_identical(filename_dangerous_type_18_fasta, filename_dangerous_type_18_fasta_after)
    sequence_analyzer.remove_identical(filename_dangerous_type_31_fasta, filename_dangerous_type_31_fasta_after)
    sequence_analyzer.remove_identical(filename_dangerous_type_33_fasta, filename_dangerous_type_33_fasta_after)
    sequence_analyzer.remove_identical(filename_dangerous_type_35_fasta, filename_dangerous_type_35_fasta_after)
    sequence_analyzer.remove_identical(filename_dangerous_type_51_fasta, filename_dangerous_type_51_fasta_after)
    sequence_analyzer.remove_identical(filename_dangerous_type_52_fasta, filename_dangerous_type_52_fasta_after)

    # define harmless type files after cd-hit processing
    filename_harmless_type_6_fasta_after = 'harmless_type_6_after.fasta'
    filename_harmless_type_11_fasta_after = 'harmless_type_11_after.fasta'
    filename_harmless_type_40_fasta_after = 'harmless_type_40_after.fasta'
    filename_harmless_type_42_fasta_after = 'harmless_type_42_after.fasta'
    filename_harmless_type_43_fasta_after = 'harmless_type_43_after.fasta'
    filename_harmless_type_44_fasta_after = 'harmless_type_44_after.fasta'
    filename_harmless_type_57_fasta_after = 'harmless_type_57_after.fasta'
    filename_harmless_type_81_fasta_after = 'harmless_type_81_after.fasta'

    # remove identical sequences with cd-hit
    sequence_analyzer.remove_identical(filename_harmless_type_6_fasta, filename_harmless_type_6_fasta_after)
    sequence_analyzer.remove_identical(filename_harmless_type_11_fasta, filename_harmless_type_11_fasta_after)
    sequence_analyzer.remove_identical(filename_harmless_type_40_fasta, filename_harmless_type_40_fasta_after)
    sequence_analyzer.remove_identical(filename_harmless_type_42_fasta, filename_harmless_type_42_fasta_after)
    sequence_analyzer.remove_identical(filename_harmless_type_43_fasta, filename_harmless_type_43_fasta_after)
    sequence_analyzer.remove_identical(filename_harmless_type_44_fasta, filename_harmless_type_44_fasta_after)
    sequence_analyzer.remove_identical(filename_harmless_type_57_fasta, filename_harmless_type_57_fasta_after)
    sequence_analyzer.remove_identical(filename_harmless_type_81_fasta, filename_harmless_type_81_fasta_after)

    # define merged dangerous and harmless types filename
    filename_merged = "merged_sequences.fasta"
    f = open(filename_merged, "w")

    # first append dangerous and only then harmless
    append_file_to_file(f, filename_dangerous_type_16_fasta_after)
    append_file_to_file(f, filename_dangerous_type_18_fasta_after)
    append_file_to_file(f, filename_dangerous_type_31_fasta_after)
    append_file_to_file(f, filename_dangerous_type_33_fasta_after)
    append_file_to_file(f, filename_dangerous_type_35_fasta_after)
    append_file_to_file(f, filename_dangerous_type_51_fasta_after)
    append_file_to_file(f, filename_dangerous_type_52_fasta_after)

    append_file_to_file(f, filename_harmless_type_6_fasta_after)
    append_file_to_file(f, filename_harmless_type_11_fasta_after)
    append_file_to_file(f, filename_harmless_type_40_fasta_after)
    append_file_to_file(f, filename_harmless_type_42_fasta_after)
    append_file_to_file(f, filename_harmless_type_43_fasta_after)
    append_file_to_file(f, filename_harmless_type_44_fasta_after)
    append_file_to_file(f, filename_harmless_type_57_fasta_after)
    append_file_to_file(f, filename_harmless_type_81_fasta_after)

    f.close()
    # align with mafft
    sequence_analyzer.align_mafft(filename_merged)

# Run main
main()
from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO
from Bio.SeqUtils import GC
import numpy as np
import pylab
import Levenshtein as levenshtein
import io


class BioSeqFinder:
    input_sequences = []
    save_file_name = ""
    start_codons = ["ATG"]
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []
    proteins = []
    input_file_format = "fasta"
    min_bp_length = 100
    table = 11
    other_orfs = []
    translated_proteins = []
    muscle_exe_location = "muscle.exe"
    muscle_aligned_output_file_name = "msa_MUSCLE_aligned_proteins.fasta"
    muscle_aligned_output_file_name_phy = "msa_MUSCLE_aligned_proteins.phy"
    align = ''
    gc = []

    # constructor that accepts input sequence files
    def __init__(self, input_data_files, protein_save_file):
        self.save_file_name = protein_save_file
        for file in input_data_files:
            self.get_sequences_from_file(file, self.input_file_format)

    def get_sequences_from_file(self, file_name, name):
        for seq_record in SeqIO.parse(file_name, name):
            self.input_sequences.append(seq_record)

    def write_proteins_to_file(self):
        seq_records = []
        i = 1
        for [j, protein] in self.proteins:
            description = "Baltymas " + str(i)
            seq_records.append(
                SeqRecord(protein, id=str(i), description=description))
            i += 1

        if len(self.proteins) > 0:
            SeqIO.write(seq_records, self.save_file_name, "fasta")

    def split_into_codones(self, seq):
        codones = [seq[i:i + 3] for i in range(0, len(seq), 3)]
        if len(codones[-1]) < 3:
            codones.pop()

        return codones

    def get_proteins(self):
        proteins = []
        i = 0
        for orf in self.orfs:
            codones = self.split_into_codones(orf)
            protein = ''
            for codone in codones:
                if codone in self.stop_codons and protein != '':
                    proteins.append([i, protein.translate(11)])
                    protein = ''
                elif protein != '' or codone in self.start_codons:
                    protein += codone
            i += 1

        self.proteins = proteins
        return proteins

    def filter_proteins_by_bp(self):
        if (len(self.proteins) > 0):
            self.proteins = [[i, protein] for [i, protein] in self.proteins if len(
                protein) >= self.min_bp_length]

    def get_all_orfs(self):
        for seq in self.input_sequences:
            self.orfs += self.__get_seq_ORFs(seq)

    def __get_seq_ORFs(self, sequence):
        orfs = []
        for strand, nuc in [(+1, sequence.seq), (-1, sequence.seq.reverse_complement())]:
            for frame in range(3):
                length = 3 * ((len(sequence)-frame) // 3)  # Multiple of three
                for orf in nuc[frame:frame+length].split("*"):
                    orfs.append(orf)
        return orfs

    def perform_MSA_on_proteins(self):
        cline = MuscleCommandline(
            self.muscle_exe_location, input=self.save_file_name, out=self.muscle_aligned_output_file_name)
        cline()
        count = AlignIO.convert(
            self.muscle_aligned_output_file_name, "fasta", self.muscle_aligned_output_file_name_phy, "phylip")
        align = AlignIO.read(
            self.muscle_aligned_output_file_name_phy, "phylip")
        self.align = align

    def calculate_distance_matrix(self):
        calculator = DistanceCalculator('blosum62')
        dm = calculator.get_distance(self.align)
        return dm

    def draw_phylo_tree_according_MSA_alignment(self):
        calculator = DistanceCalculator('blosum62')
        constructor = DistanceTreeConstructor(calculator)
        tree = constructor.build_tree(self.align)
        Phylo.draw(tree)

    def calculate_GC_percentages(self):
        gc = [GC(orf) for orf in self.orfs]
        self.gc = gc
        return gc

    def draw_histogram_for_GC(self):
        pylab.xlabel("ORF")
        pylab.xlim(0, 12)
        pylab.xticks(np.arange(0, len(self.orfs), 1.0))
        pylab.ylabel("Percentage of GC")
        pylab.grid()
        pylab.plot(self.gc)
        protein_orf_nr = [i for [i, protein] in self.proteins]
        pylab.plot(protein_orf_nr)
        pylab.savefig("gc.png")


input_data_files = ["plazmide.fasta", "test.fasta"]
protein_output_file = "proteins.fasta"

seq_analyze_instance = BioSeqFinder(input_data_files, protein_output_file)
seq_analyze_instance.get_all_orfs()
seq_analyze_instance.get_proteins()
seq_analyze_instance.filter_proteins_by_bp()
seq_analyze_instance.write_proteins_to_file()
seq_analyze_instance.perform_MSA_on_proteins()
seq_analyze_instance.calculate_distance_matrix()
seq_analyze_instance.draw_phylo_tree_according_MSA_alignment()
seq_analyze_instance.calculate_GC_percentages()
seq_analyze_instance.draw_histogram_for_GC()
print("%d", len(seq_analyze_instance.orfs))

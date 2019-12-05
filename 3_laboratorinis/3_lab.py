from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Alphabet import IUPAC
from Bio.Alphabet import ProteinAlphabet
from Bio.Align import AlignInfo
from Bio.SubsMat import FreqTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import os.path
from os import path


class SerumAlbuminAnalyzer:
    homo_sapiens_serum_albumin = {}
    homo_sapiens_serum_albumin_path = ''
    _serum_albumins = []

    def __init__(self, albumin_file):
        self._homo_sapiens_serum_albumin_path = albumin_file
        self.load_homo_sapiens_serum_albumin()
        self.blast_result = []
        self.serum_albumins = []

    def load_homo_sapiens_serum_albumin(self):
        for seq_record in SeqIO.parse(self._homo_sapiens_serum_albumin_path, "fasta"):
            self._homo_sapiens_serum_albumin = seq_record

    def blast_search_similar_sequences(self):
        if not path.exists("blast.xml"):
            result_handler = NCBIWWW.qblast(
                program="blastp",
                database="swissprot",
                sequence=self._homo_sapiens_serum_albumin.seq,
                perc_ident=80,
                entrez_query='mammals[Organism]')

            with open(".\\blast.xml", "w") as out_handle:
                out_handle.write(result_handler.read())
            result_handler.close()

        self._blast_result = SearchIO.read("blast.xml", "blast-xml")

    def filter_and_save_blast_results(self):
        with open("similar_albumins.fasta", "w") as output:
            for hit in self._blast_result.hits:
                if "Serum albumin" in hit.description:
                    self._serum_albumins.append(hit)
                    record = hit.hsps[0].hit
                    SeqIO.write(record, output, "fasta")
        

    def do_mafft_alignment(self):
        mafft_exe = "mafft-win\mafft.bat"
        in_file = "similar_albumins.fasta"
        mafft_cline = MafftCommandline(mafft_exe, input=in_file)
        stdout, stderr = mafft_cline()
        with open("mafft_msa.fasta", "w") as output:
            output.write(stdout)
    
    def draw_phylo_tree(self):
        alignment = AlignIO.read("mafft_msa.fasta", "fasta")
        alignment._alphabet = ProteinAlphabet()
        calculator = DistanceCalculator('blosum62')
        constructor = DistanceTreeConstructor(calculator)
        tree = constructor.build_tree(alignment)
        Phylo.draw(tree)

    def calculate_information_content(self):
        alignment = AlignIO.read("mafft_msa.fasta", "fasta")
        alignment._alphabet = ProteinAlphabet()
        summary_align = AlignInfo.SummaryInfo(alignment)
        msa_length = len(alignment[0])
        max_content_start_index = -1
        min_content_start_index = -1
        max_information_content = -1
        min_information_content = 999999
        i = 0
        information_contents = []
        indexes = []
        while i < msa_length - 19:
            info_content = summary_align.information_content(i, i+19)
            information_contents.append(info_content)
            indexes.append(i)
            if info_content > max_information_content:
                max_information_content = info_content
                max_content_start_index = i
            if info_content < min_information_content:
                min_information_content = info_content
                min_content_start_index = i
            i = i + 1
        plt.plot(indexes,information_contents)
        plt.ylabel('information content')
        plt.xlabel('position')
        plt.show()
        print("Similar sequence: %s ", self._homo_sapiens_serum_albumin.seq[max_content_start_index:max_content_start_index+20])
        print("Differing sequence: %s ", self._homo_sapiens_serum_albumin.seq[min_content_start_index:min_content_start_index+20])

#inicijuojam klase analizavimui
analyzer = SerumAlbuminAnalyzer('serum_albumin_homo_sapiens.fasta')
#blast panasiu zinduoliu ZSA seku paieska
analyzer.blast_search_similar_sequences()
#filtruojame ir issaugome i faila
analyzer.filter_and_save_blast_results()
#MAFFT command line MSA palyginio generavimas
analyzer.do_mafft_alignment()
#Filogenetinio medzio brezimas
analyzer.draw_phylo_tree()
#Gauname ZSA seku fragmentus, kurie panasiausi ir labiausiai skiriasi
analyzer.calculate_information_content()

#!python
import re
import numpy as np

# song@mpipz.mpg.de

class Transcript:
    chromosome_name = ""
    if_orf_conserved = True
    meta_informaiton = ""
    name = ""
    start = 1000000000
    end = 0
    strand = ""
    genome_sequence = ""
    cds_sequence = ""
    exons = np.empty([0,2], int) #start and end position
    Cds = np.empty([0,2], int)
    def __init__(self, name, strand, chromosome_name):
        self.name = name
        self.strand = strand
        self.chromosome_name = chromosome_name

    def add_cds(self, start, end):
        self.Cds = np.append(self.Cds, np.array([[start, end]]), axis=0)
        # self.Cds.append([start, end])

    def add_exon(self, start, end):
        self.exons = np.append(self.exons, np.array([[start, end]]), axis=0)

    def updateCordinate(self):
        if ("+" == self.strand):
            self.Cds=np.sort(self.Cds, axis=0 )
            self.exons = np.sort(self.exons, axis=0)
        else:
            self.Cds=-np.sort(-self.Cds, axis=0)
            self.exons = -np.sort(-self.exons, axis=0)
        if( len(self.Cds) >0 ):
            if("+" == self.strand):
                self.start = self.Cds[0][0]
                self.end = self.Cds[len(self.Cds)-1][1]
            else:
                self.start = self.Cds[len(self.Cds) - 1][0]
                self.end = self.Cds[0][1]

    def __lt__(self, other):
        if (self.start < other.start):
            return True
        if (self.start == other.start) & (self.end > other.end):
            return True
        if (self.start == other.start) & (self.end == other.end):
            return False
        return False

    def __gt__(self, other):
        if (self.start < other.start):
            return False
        if (self.start == other.start) & (self.end > other.end):
            return False
        if (self.start == other.start) & (self.end == other.end):
            return False
        return True

    def __eq__(self, other):
        if (self.start == other.start) & (self.end == other.end):
            return True
        return False


class Gene:
    name = ""
    start = 100000000
    end = 0
    transcripts = np.empty([0,1], Transcript)  # start and end position
    length = 0
    strand = ""
    def __init__(self, name, strand):
        self.name = name
        self.strand = strand
    def add_transcript(self, transcript):
        self.transcripts=np.append(self.transcripts,[transcript])

    def updateCordinate(self):
        self.transcripts = np.sort(self.transcripts)
        for transcript in self.transcripts:
            self.start = self.transcripts[0].start
            if(transcript.end > self.end):
                self.end = transcript.end
        self.length = self.end - self.start

    def __lt__(self, other):
        if (self.start < other.start):
            return True
        if (self.start == other.start) & (self.end > other.end):
            return True
        if (self.start == other.start) & (self.end == other.end):
            return False
        return False

    def __gt__(self, other):
        if (self.start < other.start):
            return False
        if (self.start == other.start) & (self.end > other.end):
            return False
        if (self.start == other.start) & (self.end == other.end):
            return False
        return True

    def __eq__(self, other):
        if (self.start == other.start) & (self.end == other.end):
            return True
        return False
#end of gene class


def readGff(gffFilePath):
    # the input if the gff file
    # return a dictionary [key(chromosome_name)]=genes  genes is a np.empty list
    chromosome_gene_dict=dict()
    chromosome_transcript_dict = dict()
    chromosome_gene_list = dict()  # start and end position
    transcript_gene_map=dict()
    gene_chr_map = dict()
    with open(gffFilePath) as f:
        for line in f:
            m = re.search('^#', line)
            if( m == None ):
                m = re.search(r'^(\S+)\t(\S+)\tCDS\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*Parent=([a-zA-Z0-9._^+-]+)', line)
                if( m != None ):
                    chromosome_name = m.group(1)
                    if (not chromosome_name in chromosome_transcript_dict):
                        chromosome_transcript_dict[chromosome_name]=dict()

                    start = int(m.group(3))
                    end = int(m.group(4))
                    if( start > end ):
                        temp = start
                        start = end
                        end = temp

                    transcript_name = m.group(8)
                    if (not transcript_name in chromosome_transcript_dict[chromosome_name]):
                        strand = m.group(6)
                        transcript = Transcript(transcript_name, strand, chromosome_name)
                        (chromosome_transcript_dict[chromosome_name])[transcript_name] = transcript

                    chromosome_transcript_dict[chromosome_name][transcript_name].add_cds(start, end)

                m = re.search(r'^(\S+)\t(\S+)\tmRNA\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*ID=([a-zA-Z0-9._^+-]+);Parent=([a-zA-Z0-9._^+-]+);', line)
                if( m != None ):
                    transcript_name = m.group(8)
                    gene_name = m.group(9)
                    transcript_gene_map[transcript_name] = gene_name
#                    print(transcript_name + "\t" + gene_name)
                m = re.search(r'^(\S+)\t(\S+)\tmRNA\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*ID=([a-zA-Z0-9._^+-]+);Parent=([a-zA-Z0-9._^+-]+)$', line)
                if (m != None):
                    transcript_name = m.group(8)
                    gene_name = m.group(9)
                    transcript_gene_map[transcript_name] = gene_name
            #                    print(transcript_name + "\t" + gene_name)
    f.close()

    with open(gffFilePath) as f:
        for line in f:
            m = re.search('^#', line)
            if( m == None ):
                m = re.search(r'^(\S+)\t(\S+)\texon\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*Parent=([a-zA-Z0-9._^+-]+)', line)
                if( m != None ):
                    chromosome_name = m.group(1)
                    if (not chromosome_name in chromosome_transcript_dict):
                        chromosome_transcript_dict[chromosome_name]=dict()

                    start = int(m.group(3))
                    end = int(m.group(4))
                    if( start > end ):
                        temp = start
                        start = end
                        end = temp

                    transcript_name = m.group(8)

                    chromosome_transcript_dict[chromosome_name][transcript_name].add_exon(start, end)
    f.close()

    #update transcript, get the start end information
    for chromosome_name in chromosome_transcript_dict:
        if (not chromosome_name in chromosome_gene_dict):
            chromosome_gene_dict[chromosome_name] = dict()
        for transcript_name in chromosome_transcript_dict[chromosome_name]:
            transcript = chromosome_transcript_dict[chromosome_name][transcript_name]
            transcript.updateCordinate()
            gene_name = transcript_gene_map[transcript_name]
            gene_chr_map[gene_name] = chromosome_name
            if not gene_name in chromosome_gene_dict[chromosome_name]:
                chromosome_gene_dict[chromosome_name][gene_name] = Gene(gene_name, transcript.strand)

            chromosome_gene_dict[chromosome_name][gene_name].add_transcript(transcript)
    #print("transcript update done")
    for chromosome_name in chromosome_gene_dict:
        gene_list = np.empty([0, 1], Gene)
        for gene_name in chromosome_gene_dict[chromosome_name]:
            chromosome_gene_dict[chromosome_name][gene_name].updateCordinate()
            gene_list = np.append(gene_list, [chromosome_gene_dict[chromosome_name][gene_name]])

        gene_list = np.sort(gene_list)
        chromosome_gene_list[chromosome_name] = gene_list

    return chromosome_gene_dict, chromosome_gene_list, transcript_gene_map, gene_chr_map

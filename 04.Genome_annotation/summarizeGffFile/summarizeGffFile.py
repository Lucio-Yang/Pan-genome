#!python

import GffFile
import MyUtil
import argparse

# song@mpipz.mpg.de


def read_data(paspGffFile):
    chromosome_gene_dict, chromosome_gene_list, transcript_gene_map, gene_chr_map = GffFile.readGff(paspGffFile)
    numberOfGene = 0
    numberOfTranscript = 0
    totalTranscriptLength = 0  # CDA+intron
    totalCdsLength = 0
    totalCdsNumber = 0
    totalExonLength = 0
    totalExonNumber = 0
    totalIntronLength = 0
    totalIntronNumber = 0
    for chr in chromosome_gene_dict:
        for geneName in chromosome_gene_dict[chr]:
            gene = chromosome_gene_dict[chr][geneName]
            numberOfGene = numberOfGene + 1
            for transcript in gene.transcripts:
                if transcript.name.endswith(".mRNA1"):
                    numberOfTranscript = numberOfTranscript + 1
                    totalTranscriptLength = totalTranscriptLength + transcript.end - transcript.start + 1
                    thisCdsNumber = 0
                    thisCdsTotalLength = 0
                    for cds in transcript.Cds:
                        totalCdsNumber = totalCdsNumber + 1
                        totalCdsLength = totalCdsLength + cds[1] - cds[0] + 1
                        thisCdsNumber = thisCdsNumber + 1
                        thisCdsTotalLength = thisCdsTotalLength + cds[1] - cds[0] + 1
                    if ("+" == transcript.strand):
                        totalIntronLength = totalIntronLength + transcript.exons[len(transcript.exons) - 1][1] - transcript.exons[0][0] + 1 - thisCdsTotalLength
                    else:
                        totalIntronLength = totalIntronLength +  transcript.exons[0][1] - transcript.exons[len(transcript.exons) - 1][0] + 1 - thisCdsTotalLength

                    totalIntronNumber = totalIntronNumber + thisCdsNumber - 1

    print(str(numberOfGene) + "\t" + str(
        round(float(totalTranscriptLength) / float(numberOfTranscript), 2)) + "\t" + str(
        round(float(totalCdsLength) / float(numberOfTranscript), 2)) + "\t" + str(
        round(float(totalCdsLength) / float(totalCdsNumber), 2)) +
          "\t" + str(round(float(totalIntronLength) / float(numberOfTranscript), 2)) + "\t" + str(
        round(float(totalCdsNumber) / float(numberOfGene), 2)))

                    # thisExonNumber = 0
                    # thisExonTotalLength = 0
                    # for exon in transcript.exons:
                    #     thisExonNumber = thisExonNumber + 1
                    #     totalExonNumber = totalExonNumber + 1
                    #     totalExonLength = totalExonLength + exon[1] - exon[0] + 1
                    #     thisExonTotalLength = thisExonTotalLength + exon[1] - exon[0] + 1
                    # if ("+" == transcript.strand):
                    #     totalIntronLength = totalIntronLength + transcript.exons[len(transcript.exons) - 1][1] - transcript.exons[0][0] + 1 - thisExonTotalLength
                    # else:
                    #     totalIntronLength = totalIntronLength +  transcript.exons[0][1] - transcript.exons[len(transcript.exons) - 1][0] + 1 - thisExonTotalLength
                    # totalIntronNumber = totalIntronNumber + thisExonNumber - 1

#    print("Number of transcripts:" + str(numberOfTranscript))
#     print("Number of genes:" + str (numberOfGene))
#     print("CDS+intron len(avg, bp):" + str(round(float(totalTranscriptLength)/float(numberOfTranscript),2)))
#     print("total CDS length per gene(avg, bp):" + str(round(float(totalCdsLength) / float(numberOfTranscript),2)))
#     print("CDS length(avg, bp):" + str(round(float(totalCdsLength) / float(totalCdsNumber),2)))
#     print("intron length in CDS region per gene(avg):" + str(round(float(totalIntronLength) / float(numberOfTranscript),2)))
#     print("CDSs per gene(avg):" + str(round(float(totalCdsNumber) / float(numberOfGene),2)))


#
# Kronos	TW_t2
# IG99236	TW_t9
# PI294478	TW_t15
# XM000834	TW_834
# NU01905	TW_1905
# NU00021	TW_t21
# PI330553	TW_t7
# IG77365	TW_t5
# XM001097	TW_1097
# XM000719	TW_719
# NU01954	TW_1954
# Zavitan	TW_zavitan



print("Number of genes\t" + "CDS+intron len(avg, bp)\t" + "total CDS length per gene(avg, bp)" + "\tCDS length(avg, bp)\t" + "intron length in CDS region per gene(avg)\t" + "CDSs per gene(avg)")
#Kronos
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t2/TW_t2.HC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t2/TW_t2.LC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t2/TW_t2.TE.gff3")

#XM001097
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/1097/TW_1097.HC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/1097/TW_1097.LC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/1097/TW_1097.TE.gff3")

#XM000834
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/834/TW_834.HC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/834/TW_834.LC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/834/TW_834.TE.gff3")

#XM000719
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/719/TW_719.HC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/719/TW_719.LC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/719/TW_719.TE.gff3")


#NU00021
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t21/TW_t21.HC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t21/TW_t21.LC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t21/TW_t21.TE.gff3")

#IG77365
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t5/TW_t5.HC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t5/TW_t5.LC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t5/TW_t5.TE.gff3")


#PI294478
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t15/TW_t15.HC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t15/TW_t15.LC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t15/TW_t15.TE.gff3")


#PI330553
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t7/TW_t7.HC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t7/TW_t7.LC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t7/TW_t7.TE.gff3")


#IG99236
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t9/TW_t9.HC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t9/TW_t9.LC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/t9/TW_t9.TE.gff3")

#NU01954
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/1954/TW_1954.HC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/1954/TW_1954.LC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/1954/TW_1954.TE.gff3")

#NU01905
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/1905/TW_1905.HC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/1905/TW_1905.LC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/1905/TW_1905.TE.gff3")

#Zavitan
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/zavitan/TW_zavitan.HC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/zavitan/TW_zavitan.LC.gff3")
read_data("/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/07_diamondAgainstHomo/zavitan/TW_zavitan.TE.gff3")


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='Summarize GffFile')
#     parser.add_argument('-i', '--input', dest='input', type=str, help="The input gff file.", required=True)
#     parser.add_argument('-o', '--output', dest='output', type=str, help="The output summarize file.", required=True)
#     args = parser.parse_args()
#     read_data(args.input, args.input)

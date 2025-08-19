#!python

import GffFile
import MyUtil
import argparse

# song@mpipz.mpg.de

_buckets = []

def read_data(paspGffFile, transdecoderGff, deleteList, putbackList, pasa_busco, transdecoder_busco):
    pasa_chromosome_gene_dict, pasa_chromosome_gene_list, pasa_transcript_gene_map, pasa_gene_chr_map = GffFile.readGff( paspGffFile )
    transdecoder_chromosome_gene_dict, transdecoder_chromosome_gene_list, transdecoder_transcript_gene_map, transdecoder_gene_chr_map = GffFile.readGff(transdecoderGff)
    pasa_missing_busco_ids = set()

    with open(pasa_busco) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                pasa_missing_busco_ids.add(line)
            #print(line)
    f.close()

    put_back_pasa_missing_gene_ids = set()
    with open(transdecoder_busco) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                elements = line.split()
                if elements[1] != "Missing":
                    if elements[0] in pasa_missing_busco_ids:
                        if elements[2] in transdecoder_transcript_gene_map:
                            put_back_pasa_missing_gene_ids.add(transdecoder_transcript_gene_map[elements[2]])

    f.close()

    all_the_overlapped_genes = set()
    for gene in put_back_pasa_missing_gene_ids:
        overlapped_genes = MyUtil.overlap_certain_gene_with_genes(transdecoder_chromosome_gene_dict[transdecoder_gene_chr_map[gene]][gene], transdecoder_gene_chr_map[gene], pasa_chromosome_gene_dict)
        all_the_overlapped_genes.update(overlapped_genes)

    f = open(putbackList, "w")
    for gene in put_back_pasa_missing_gene_ids:
        f.write(gene + "\n")
    f.close()

    f = open(deleteList, "w")
    for gene in all_the_overlapped_genes:
        f.write(gene + "\n")
    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Put missing BUsco gene back')
    parser.add_argument('-p', '--pasa', dest='pasa_Gff_file', type=str, help="The gff file from two rounds of PASA", required=True)
    parser.add_argument('-t', '--transdecoder', dest='transdecoder_Gff_file', type=str, help="The gff file from transdecoder", required=True)
    parser.add_argument('-d', '--deleted', dest='delete_list_from_pasa', type=str, help="The gene list to be deleted from pasa", required=True)
    parser.add_argument('-b', '--putback', dest='putback_list_from_transdecoder', type=str, help="putback_list_from_transdecoder", required=True)
    parser.add_argument('-summarizeGffFile', '--pasa_busco', dest='pasa_busco', type=str, help="pasa busco results", required=True)
    #/home/bs674/HPC2/data/yangg/project/GeneAnnot/09.gene_classifying/01_predictResults/BUSCO_t2/run_poales_odb10/missing_busco_list.tsv
    parser.add_argument('-r', '--transdecoder_busco', dest='transdecoder_busco', type=str, help="transdecoder busco results", required=True)
    #/home/bs674/HPC2/data/yangg/project/GeneAnnot/04.transcript_evidence/01_merge_and_transdecoder.eval/TW_t2/BUSCO/run_poales_odb10/full_table.tsv

#    parser.add_argument('-summarizeGffFile', '--minSimilarity', dest='minSimilarity', type=float, default=0.5, help="The minimum similarity for an alignment to be used for upstream and downstream anchors.")
    args = parser.parse_args()
    read_data(args.pasa_Gff_file, args.transdecoder_Gff_file, args.delete_list_from_pasa, args.putback_list_from_transdecoder, args.pasa_busco, args.transdecoder_busco)



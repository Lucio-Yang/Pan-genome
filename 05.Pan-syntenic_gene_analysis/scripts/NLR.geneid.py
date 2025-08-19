import os
import re
import shutil
from collections import defaultdict

def setup_directories(base_dir):
    """创建所有需要的目录结构"""
    dirs = {
        'nlr_ann': os.path.join(base_dir, "0.NLR_Ann"),
        'gff3': os.path.join(base_dir, "0.gff3"),
        'nlr_id': os.path.join(base_dir, "1.NLR_id"),
        'nlr_id_sort': os.path.join(base_dir, "1.NLR_id_sort"),
    }
    
    for dir_path in dirs.values():
        os.makedirs(dir_path, exist_ok=True)
        print(f"确保目录存在: {dir_path}")
    
    return dirs

def log_step(step, message):
    """记录处理步骤"""
    print(f"\n{'='*50}")
    print(f"步骤 {step}: {message}")
    print(f"{'='*50}")

def parse_nlr_gff(nlr_gff_file):
    """解析NLR GFF文件"""
    nlr_data = []
    try:
        with open(nlr_gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                attributes = fields[8]
                match = re.search(r'nlrClass=([^;]+)', attributes)
                nlr_class = match.group(1) if match else "Unknown"
                nlr_data.append((chrom, start, end, nlr_class))
        print(f"成功解析 {len(nlr_data)} 条记录来自 {nlr_gff_file}")
    except Exception as e:
        print(f"解析 {nlr_gff_file} 时出错: {str(e)}")
    return nlr_data

def parse_species_gff(species_gff_file):
    """解析物种GFF3文件"""
    species_data = {}
    try:
        with open(species_gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'gene':
                    continue
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                attributes = fields[8]
                match = re.search(r'ID=(?:gene:)?([^;]+)', attributes)
                gene_id = match.group(1) if match else "Unknown"
                if chrom not in species_data:
                    species_data[chrom] = []
                species_data[chrom].append((start, end, gene_id))
        print(f"成功解析 {sum(len(v) for v in species_data.values())} 条基因记录来自 {species_gff_file}")
    except Exception as e:
        print(f"解析 {species_gff_file} 时出错: {str(e)}")
    return species_data

def find_matches(nlr_data, species_data):
    """查找匹配的基因"""
    results = []
    for chrom, start, end, nlr_class in nlr_data:
        if chrom in species_data:
            for mrna_start, mrna_end, gene_id in species_data[chrom]:
                if not (end < mrna_start or start > mrna_end):
                    results.append((gene_id, start, end, nlr_class))
    print(f"找到 {len(results)} 个匹配基因")
    return results

def process_gene_matching(dirs, species_list):
    """处理基因匹配步骤"""
    log_step(1, "处理基因匹配")
    
    for species in species_list:
        print(f"\n处理物种: {species}")
        nlr_gff = os.path.join(dirs['nlr_ann'], f"{species}.NLR.gff")
        species_gff = os.path.join(dirs['gff3'], f"{species}.HC.gff3")
        output_file = os.path.join(dirs['nlr_id'], f"{species}.idgene.txt")
        
        if not os.path.exists(nlr_gff):
            print(f"警告: 未找到NLR GFF文件: {nlr_gff}")
            continue
        if not os.path.exists(species_gff):
            print(f"警告: 未找到物种GFF文件: {species_gff}")
            continue
        
        nlr_data = parse_nlr_gff(nlr_gff)
        species_data = parse_species_gff(species_gff)
        matches = find_matches(nlr_data, species_data)
        
        with open(output_file, 'w') as f:
            f.write("Gene_ID\tStart\tEnd\tNLR_Class\n")
            for gene_id, start, end, nlr_class in matches:
                f.write(f"{gene_id}\t{start}\t{end}\t{nlr_class}\n")
        print(f"已写入匹配结果到: {output_file}")

def main():
    # 基础目录
    base_dir = "./"
    
    # 设置目录结构
    dirs = setup_directories(base_dir)
    
    # 物种列表
    species_list = ["sample1","sample2"]
    
    # 执行处理流程
    process_gene_matching(dirs, species_list)
    
    print("\n所有处理步骤已完成!")

if __name__ == "__main__":
    main()


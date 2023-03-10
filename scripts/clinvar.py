import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", '--num_variants', type=int, default=500)
args = parser.parse_args()

def get_mut_info(change, type):
    if type == "nt":
        match = re.match("c\.(.+)(\w)>(\w)", change)
        if not match:
            return None
        pos = match.group(1)
        ref = match.group(2)
        alt = match.group(3)
    elif type == "protein":
        match = re.match("p\.(\D+)(\d+)(\D+)", change)
        if not match:
            return None
        pos = match.group(2)
        ref = match.group(1)
        alt = match.group(3)
    return pos, ref, alt

def add_snv(variant_id, name, gene, build, chr, genomic_pos):
    try:
        match = re.match("(.+)\((.+)\)\:(.+).\((.+)\)", name)
        protein_change = match.group(4)
        protein_pos, ref_aa, alt_aa = get_mut_info(protein_change, "protein")
    except:
        match = re.match("(.+)\((.+)\)\:(.+)", name)
        protein_pos, ref_aa, alt_aa = "NA", "NA", "NA"
        if not match:
            return
    transcript = match.group(1)
    gene = match.group(2) 
    coding_change = match.group(3)
    try:
        coding_pos, ref_nt, alt_nt = get_mut_info(coding_change, "nt")
    except:
        return
    result = [variant_id, "clinvar", build, gene, transcript, chr, genomic_pos, coding_pos, ref_nt, alt_nt, protein_pos, ref_aa, alt_aa]
    return ",".join(result) + "\n"

def add_cna(variant_id, name, gene, build, chr, gain_or_loss):
    try:
        match = re.match("(.+)\s(.+)\((.+)\)x(\d+)", name)
        copies = match.group(4)
    except:
        match = re.match("(.+)\s(.+)\((.+)\)", name)
        copies = "NA"
        if not match:
            return
    cytogenetic_loc = match.group(2)
    genomic_loc = match.group(3)
    try:
        match = re.match("(.+):(.+)-(.+)", genomic_loc)
        genomic_start = match.group(2)
        genomic_end = match.group(3)
    except:
        return
    result = [variant_id, "clinvar", build, gene, copies, gain_or_loss, chr, genomic_start, genomic_end, cytogenetic_loc]
    return ",".join(result) + "\n"

def add_statement(phenotypes, clinsig, disease_outfile, statement_outfile, involves_outfile):
    global variant_id, disease_id, diseases, statement_id
    for pheno in phenotypes:
        if pheno == "not provided" or pheno == "not specified":
            continue
        if pheno in diseases:
            id = diseases.index(pheno)
        else:
            disease_outfile.write(",".join([str(disease_id), pheno]) + " \n")
            id = disease_id
            disease_id += 1
            diseases.append(pheno)
        statement_outfile.write(",".join([str(statement_id), clinsig]) + "\n")
        involves_outfile.write(",".join([str(statement_id), str(variant_id), str(id), ""]) + "\n")
        statement_id += 1

def main():
    global variant_id, disease_id, diseases, statement_id

    infile = open("/Users/ninaxiong/Desktop/clinvar_variant_summary.txt") # change this
    snv_outfile = open("../csv/snv.csv", "a")
    cna_outfile = open("../csv/cna.csv", "a")
    disease_outfile = open("../csv/disease.csv", "a")
    statement_outfile = open("../csv/statement.csv", "a")
    involves_outfile = open("../csv/involves.csv", "a")

    while variant_id < args.num_variants:
        line = infile.readline()
        line = line.strip().split("\t")
        # extract data
        type, name, geneid, gene = line[1], line[2], line[3], line[4]
        chr, start, stop = line[18], line[19], line[20]
        clinsig = line[6]
        phenotypes = line[13].strip().split("|")
        build = line[16]
        # add variants
        if build != "GRCh38":
            continue
        if type == "single nucleotide variant":
            snv = add_snv(str(variant_id), name, gene, build, chr, start)
            if snv:
                print(name)
                snv_outfile.write(snv)
                add_statement(phenotypes, clinsig, disease_outfile, statement_outfile, involves_outfile)
                variant_id += 1
        if type.startswith("copy number"):
            gain_or_loss = re.match("copy number (.+)", type).group(1)
            if geneid == "-1":
                gene = "-"
            cna = add_cna(str(variant_id), name, gene, build, chr, gain_or_loss)
            if cna:
                cna_outfile.write(cna)
                add_statement(phenotypes, clinsig, disease_outfile, statement_outfile, involves_outfile)
                variant_id += 1

    infile.close()
    snv_outfile.close()
    cna_outfile.close()
    disease_outfile.close()
    statement_outfile.close()
    involves_outfile.close()

variant_id = 0
disease_id = 0
diseases = []
statement_id = 0

if __name__ == "__main__":
    main()
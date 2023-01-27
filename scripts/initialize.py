def create_file(filename, fields):
    outfile = open(filename, "w")
    header = ",".join(fields) + "\n"
    outfile.write(header)
    return outfile

def main():
    files = {"../csv/snv.csv": ["variant_id", "source", "build", "gene", "transcript", "chromosome", "genomic_pos", "coding_pos", "ref_nt", 
    "alt_nt", "protein_pos", "ref_aa", "alt_aa"],
        "../csv/cna.csv": ["variant_id", "source", "build", "gene", "copies", "gain_or_loss", "chromosome", "genomic_start", "genomic_end", 
    "cytogenetic_loc"],
        "../csv/disease.csv": ["disease_id", "name"],
        "../csv/statement.csv": ["statement_id", "significance"],
        "../csv/involves.csv": ["statement_id", "variant_id", "disease_id", "drug_id"],
        "../csv/drug.csv": ["drug_id, name"]
    }
    for filename, fields in files.items():
        create_file(filename, fields)

if __name__ == "__main__":
    main()
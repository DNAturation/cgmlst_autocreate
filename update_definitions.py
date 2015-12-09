import json
import os
import argparse
from Bio import SeqIO

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--alleles', required = True, help = "Allele directory")
    parser.add_argument('-j', '--jsons', required = True, help = "JSON directory")
    parser.add_argument('-t', '--test', required = True, help = "Test name")

    return parser.parse_args()

def load_data(jsonpath):
    
    with open(jsonpath, 'r') as f:
        data = json.load(f)
    return data

def get_known_alleles(allele_dir):
    
    # strips path and extension
    get_name = lambda x: os.path.basename(os.path.splitext(x)[0])

    known_alleles = {}

    for (root, directory, fname) in os.walk(allele_dir):

        with open(fname, 'r') as f:
            alleles = [str(x) for x in SeqIO.parse(f, 'fasta')]

        known_alleles[get_name(fname)] = alleles

    return known_alleles

def update_alleles(known_alleles, allele_dir):
    """Overwrites allele files with updated definitions."""

    for gene in known_alleles:

        fname = os.path.join(allele_dir, gene + ".fasta")

        with open(fname, 'w') as f:
            out = ""
            counter = 0
            for allele in known_alleles[gene]:
                out += ">{}\n".format(counter + 1)
                out += known_alleles[gene][counter]

            f.write(out)
                
def update(json_dir, known_alleles, test):
    """Reads JSON output from MIST.
    
    If a new allele is discovered, it is appended to the list of known alleles
    and the JSON is updated accordingly.
    
    Returns dict of known alleles and JSON objects for later writing.
    """

    json_out = {}

    for (root, director, fname) in os.walk(json_dir):
        
        json_name = os.path.join(root, fname)
        data = load_data(json_name)

        genes = data["Results"][0]["TestResults"][test]

        for g in genes:
            
            gene = genes[g]

            if gene["BlastResults"] is not None and not gene["CorrectMarkerMatch"]:
                br = gene["BlastResults"]

                sj = br["SubjAln"].replace('-', '')
                
                br["QueryAln"] = sj 
                br["SubjAln"] = sj 
                
                known_alleles[g].append(sj)

                br["Mismatches"] = 0
                br["Gaps"] = 0
                br["PercentIdentity"] = 100.0
                gene["Mismatches"] = 0
                gene["BlastPercentIdentity"] = 100.0
                gene["CorrectMarkerMatch"] = True

                gene["MarkerCall"] = str(known_alleles[g].index(br["SubjAln"]) + 1)
                gene["AlleleMatch"] = gene["MarkerCall"]

            json_out[json_name] = data
   
   return known_alleles, json_out

def main():

    args = arguments()

    known_alleles = get_known_alleles(args.alleles)
    
    known_alleles, json_out = update(args.jsons, known_alleles, args.test)

    update_alleles(known_alleles, args.allele_dir)
    
    for j in json_out:
        with open(j, 'w') as f:
            json.dump(json_out[j], f)

if __name__ == '__main__':
    main()

from sys import argv


def main():
    records = {}
    input_file = argv[1]
    output_fname = argv[2]
    with open(input_file) as input_fhand:
        for line in input_fhand:
            if line.startswith("!"):
                continue
            else:
                line = line.rstrip().split()
                gene_name = line[1]
                go_term = line[4]
                if gene_name not in records:
                    records[gene_name] = [go_term]
                else:
                    records[gene_name].append(go_term)
    with open(output_fname, "w") as output_fhand:
        for gene, goterms in records.items():
            line = gene+"\t"+";".join(goterms)+"\n"
            output_fhand.write(line)
            output_fhand.flush()


if __name__ == "__main__":
    main()
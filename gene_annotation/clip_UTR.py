import argparse
import pandas as pd

def filter_gff3_remove_utrs(gff3):
    """
    Keep only gene, mRNA, and CDS features.
    Remove UTRs and other feature types.
    """
    rows = []
    with open(gff3, "r") as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            if feature_type not in {"gene", "mRNA", "CDS"}:
                continue

            attributes = parts[8]
            attr_dict = dict(
                item.split("=", 1) for item in attributes.split(";") if "=" in item
            )

            gene_id = None
            if feature_type == "gene":
                gene_id = attr_dict.get("ID")
            else:
                parent = attr_dict.get("Parent")
                gene_id = parent.split("-")[0] if parent else None

            rows.append(parts + [gene_id])

    return pd.DataFrame(
        rows,
        columns=[
            "seqid", "source", "type", "start", "end",
            "score", "strand", "phase", "attr", "gene"
        ],
    )


def correct_gene_and_mRNA_boundaries(gff3_df):
    df = gff3_df.copy()
    result_rows = []

    # clean attributes
    gene_mask = df["type"] == "gene"
    df.loc[gene_mask, "attr"] = df.loc[gene_mask, "attr"].apply(
        lambda x: ";".join(p for p in x.split(";") if not p.startswith("Name="))
    )

    cds_mask = df["type"] == "CDS"
    df.loc[cds_mask, "attr"] = df.loc[cds_mask, "attr"].str.replace(',""', '', regex=False)

    for gene_id, group in df.groupby("gene", sort=False):
        group = group.reset_index(drop=True)

        cds_rows = group[group["type"] == "CDS"]

        # duplicate CDS as exon
        exon_rows = cds_rows.copy()
        exon_rows["type"] = "exon"
        exon_rows["attr"] = exon_rows["attr"].str.replace(".cds", ".exon", regex=False)

        if not exon_rows.empty:
            strand = exon_rows["strand"].iloc[0]
            exon_rows = exon_rows.sort_values(
                by=["start", "end"], ascending=(strand == "-")
            ).reset_index(drop=True)

            for i in range(len(exon_rows)):
                exon_rows.loc[i, "attr"] = ";".join(
                    p.replace(".exon", f".exon{i+1}") if p.startswith("ID=") else p
                    for p in exon_rows.loc[i, "attr"].split(";")
                )

        group = pd.concat([group, exon_rows], ignore_index=True)

        if not cds_rows.empty:
            min_start = cds_rows["start"].astype(int).min()
            max_end = cds_rows["end"].astype(int).max()

            group.loc[group["type"] == "mRNA", ["start", "end"]] = [min_start, max_end]
            group.loc[group["type"] == "gene", ["start", "end"]] = [min_start, max_end]

        result_rows.append(group)

    return pd.concat(result_rows, ignore_index=True)


def main(gff3_file, output_file):
    filtered_df = filter_gff3_remove_utrs(gff3_file)
    corrected_df = correct_gene_and_mRNA_boundaries(filtered_df)
    corrected_df.iloc[:, :-1].to_csv(
        output_file,
        sep="\t",
        header=False,
        index=False,
        quoting=3,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove UTRs and update gene/mRNA coordinates using CDS boundaries."
    )
    parser.add_argument("-i", "--input", required=True, help="Input GFF3 file")
    parser.add_argument("-o", "--output", required=True, help="Output GFF3 file")

    args = parser.parse_args()
    main(args.input, args.output)

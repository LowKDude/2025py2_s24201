# File: s24201_2025-2.py
# Folder: 2025py2_s24201
# Author: s24201
# Purpose: Retrieve GenBank data from NCBI using a taxid, filter by sequence length,
# generate a CSV summary and a PNG chart of record lengths.

from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt

# === CONFIG ===
Entrez.tool = 'BioScriptEx10'

# === CLASS ===
class NCBIRetriever:
    def __init__(self, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key
        self.webenv = None
        self.query_key = None
        self.count = 0

    def search(self, taxid):
        try:
            search_term = f"txid{taxid}[Organism]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            result = Entrez.read(handle)
            self.webenv = result["WebEnv"]
            self.query_key = result["QueryKey"]
            self.count = int(result["Count"])
            return self.count
        except Exception as e:
            print("Search error:", e)
            return 0

    def fetch_batch(self, retstart, retmax):
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=retstart,
                retmax=retmax,
                webenv=self.webenv,
                query_key=self.query_key
            )
            return list(SeqIO.parse(handle, "genbank"))
        except Exception as e:
            print("Fetch error:", e)
            return []

# === MAIN FUNCTION ===
def main():
    email = input("Enter your NCBI email: ")
    api_key = input("Enter your NCBI API key: ")
    taxid = input("Enter the taxonomic ID (taxid): ")
    min_len = int(input("Enter minimum sequence length: "))
    max_len = int(input("Enter maximum sequence length: "))

    retriever = NCBIRetriever(email, api_key)
    count = retriever.search(taxid)
    if count == 0:
        print("No records found.")
        return

    print(f"Found {count} records. Downloading and filtering...")

    all_records = []
    for start in range(0, count, 500):
        batch = retriever.fetch_batch(start, 500)
        for record in batch:
            length = len(record.seq)
            if min_len <= length <= max_len:
                acc = record.id
                desc = record.description
                all_records.append((acc, length, desc))
        if len(all_records) >= 1000:  # limit for safety
            break

    if not all_records:
        print("No records in specified length range.")
        return

    df = pd.DataFrame(all_records, columns=["Accession", "Length", "Description"])
    csv_file = f"taxid_{taxid}_filtered.csv"
    df.to_csv(csv_file, index=False)
    print(f"CSV report saved to {csv_file}")

    # Sort and plot
    df_sorted = df.sort_values(by="Length", ascending=False)
    plt.figure(figsize=(12, 6))
    plt.plot(df_sorted["Accession"], df_sorted["Length"], marker='o')
    plt.xticks(rotation=90, fontsize=6)
    plt.xlabel("Accession")
    plt.ylabel("Sequence Length")
    plt.title(f"Sequence Lengths for TaxID {taxid}")
    plt.tight_layout()
    img_file = f"taxid_{taxid}_chart.png"
    plt.savefig(img_file)
    print(f"Chart saved to {img_file}")

if __name__ == "__main__":
    main()

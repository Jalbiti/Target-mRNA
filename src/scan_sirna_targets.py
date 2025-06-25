import os
import re
import subprocess
from Bio import SeqIO
import pandas as pd

# Parameters
WINDOW_SIZE = 21
RNAplfold_WINDOW = 80
RNAplfold_MAXBP = 40
RNAplfold_UP = 21

def run_rnaplfold(seq, transcript_id):
    with open("temp_rnaplfold_input.fa", "w") as f:
        f.write(f">{transcript_id}\n{seq}\n")

    subprocess.run(
        ["RNAplfold", "-W", str(RNAplfold_WINDOW), "-L", str(RNAplfold_MAXBP), "-u", str(RNAplfold_UP)],
        input=open("temp_rnaplfold_input.fa").read(),
        text=True,
        check=True,
    )

    lunp_file = f"{transcript_id}_lunp"
    if not os.path.exists(lunp_file):
        raise FileNotFoundError("RNAplfold output not found")

    accessibility_scores = {}
    with open(lunp_file) as fh:
        header = next(fh)                         # Skip the header line
        for line in fh:
            start, *probs = line.rstrip().split("\t")
            accessibility_scores[int(start)] = [float(p) for p in probs]

    return accessibility_scores

def extract_windows(seq, accessibility_scores):
    candidates = []
    for i in range(len(seq) - WINDOW_SIZE + 1):
        window_seq = seq[i:i+WINDOW_SIZE]
        gc_content = (window_seq.count("G") + window_seq.count("C")) / WINDOW_SIZE
        prob_unpaired = accessibility_scores.get(i + 1, [0.0])[WINDOW_SIZE - 1]
        candidates.append({
            "start": i,
            "sequence": window_seq,
            "gc_content": round(gc_content, 3),
            "unpaired_prob": round(prob_unpaired, 4)
        })
    return pd.DataFrame(candidates).sort_values(by="unpaired_prob", ascending=False)

def cleanup_temp():
    for f in os.listdir():
        if re.search(r'_lunp$|_dp\.ps|temp_rnaplfold_input\.fa', f):
            os.remove(f)

def main(fasta_path):
    record = list(SeqIO.parse(fasta_path, "fasta"))[0]
    seq = str(record.seq).upper()
    transcript_id = record.id

    scores = run_rnaplfold(seq, transcript_id)
    df = extract_windows(seq, scores)
    df.to_csv("sirna_candidates.csv", index=False)
    print(df.head(10))
    cleanup_temp()

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python scan_sirna_targets.py path/to/input.fasta")
    else:
        main(sys.argv[1])

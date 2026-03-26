# ============================================
# Sequence Alignment Tool
# Author: Padma Shree
# Phase 2 - Project 1
# Algorithms: Needleman-Wunsch (Global)
#             Smith-Waterman (Local)
# ============================================

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from datetime import datetime

# ─────────────────────────────────────────────
# SCORING SCHEME
# ─────────────────────────────────────────────

MATCH    =  2
MISMATCH = -1
GAP      = -2

# ─────────────────────────────────────────────
# NEEDLEMAN-WUNSCH — GLOBAL ALIGNMENT
# ─────────────────────────────────────────────

def needleman_wunsch(seq1, seq2):
    m, n   = len(seq1), len(seq2)
    matrix = np.zeros((m+1, n+1), dtype=int)

    # Initialize
    for i in range(m+1):
        matrix[i][0] = i * GAP
    for j in range(n+1):
        matrix[0][j] = j * GAP

    # Fill
    for i in range(1, m+1):
        for j in range(1, n+1):
            match  = matrix[i-1][j-1] + (MATCH if seq1[i-1] == seq2[j-1] else MISMATCH)
            delete = matrix[i-1][j] + GAP
            insert = matrix[i][j-1] + GAP
            matrix[i][j] = max(match, delete, insert)

    # Traceback
    aligned1, aligned2, matches = [], [], []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + (MATCH if seq1[i-1] == seq2[j-1] else MISMATCH):
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            matches.append("|" if seq1[i-1] == seq2[j-1] else ".")
            i -= 1; j -= 1
        elif i > 0 and matrix[i][j] == matrix[i-1][j] + GAP:
            aligned1.append(seq1[i-1])
            aligned2.append("-")
            matches.append(" ")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j-1])
            matches.append(" ")
            j -= 1

    a1 = "".join(reversed(aligned1))
    a2 = "".join(reversed(aligned2))
    ms = "".join(reversed(matches))

    score      = matrix[m][n]
    match_count   = ms.count("|")
    mismatch_count = ms.count(".")
    gap_count     = ms.count(" ")
    identity      = round((match_count / len(a1)) * 100, 2)
    similarity    = round(((match_count + mismatch_count) / len(a1)) * 100, 2)

    return {
        "type"       : "Global (Needleman-Wunsch)",
        "aligned1"   : a1,
        "aligned2"   : a2,
        "matches"    : ms,
        "score"      : score,
        "length"     : len(a1),
        "matches_count"  : match_count,
        "mismatches"     : mismatch_count,
        "gaps"           : gap_count,
        "identity"       : identity,
        "similarity"     : similarity,
        "matrix"         : matrix
    }

# ─────────────────────────────────────────────
# SMITH-WATERMAN — LOCAL ALIGNMENT
# ─────────────────────────────────────────────

def smith_waterman(seq1, seq2):
    m, n   = len(seq1), len(seq2)
    matrix = np.zeros((m+1, n+1), dtype=int)

    max_score = 0
    max_pos   = (0, 0)

    # Fill
    for i in range(1, m+1):
        for j in range(1, n+1):
            match  = matrix[i-1][j-1] + (MATCH if seq1[i-1] == seq2[j-1] else MISMATCH)
            delete = matrix[i-1][j] + GAP
            insert = matrix[i][j-1] + GAP
            matrix[i][j] = max(0, match, delete, insert)
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos   = (i, j)

    # Traceback from max score
    aligned1, aligned2, matches = [], [], []
    i, j = max_pos
    while i > 0 and j > 0 and matrix[i][j] > 0:
        if matrix[i][j] == matrix[i-1][j-1] + (MATCH if seq1[i-1] == seq2[j-1] else MISMATCH):
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            matches.append("|" if seq1[i-1] == seq2[j-1] else ".")
            i -= 1; j -= 1
        elif matrix[i][j] == matrix[i-1][j] + GAP:
            aligned1.append(seq1[i-1])
            aligned2.append("-")
            matches.append(" ")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j-1])
            matches.append(" ")
            j -= 1

    a1 = "".join(reversed(aligned1))
    a2 = "".join(reversed(aligned2))
    ms = "".join(reversed(matches))

    match_count    = ms.count("|")
    mismatch_count = ms.count(".")
    gap_count      = ms.count(" ")
    identity       = round((match_count / max(len(a1), 1)) * 100, 2)

    return {
        "type"          : "Local (Smith-Waterman)",
        "aligned1"      : a1,
        "aligned2"      : a2,
        "matches"       : ms,
        "score"         : max_score,
        "length"        : len(a1),
        "matches_count" : match_count,
        "mismatches"    : mismatch_count,
        "gaps"          : gap_count,
        "identity"      : identity,
        "similarity"    : identity,
        "matrix"        : matrix
    }

# ─────────────────────────────────────────────
# BLAST-LIKE OUTPUT
# ─────────────────────────────────────────────

def blast_output(result, seq1, seq2):
    print("\n" + "="*65)
    print(f"  Alignment Type : {result['type']}")
    print("="*65)
    print(f"  Score          : {result['score']}")
    print(f"  Length         : {result['length']}")
    print(f"  Identities     : {result['matches_count']}/{result['length']} ({result['identity']}%)")
    print(f"  Mismatches     : {result['mismatches']}/{result['length']}")
    print(f"  Gaps           : {result['gaps']}/{result['length']}")
    print(f"  Similarity     : {result['similarity']}%")
    print()

    # Print alignment in blocks of 60
    block = 60
    for start in range(0, result['length'], block):
        end   = min(start + block, result['length'])
        q_seq = result['aligned1'][start:end]
        m_seq = result['matches'][start:end]
        s_seq = result['aligned2'][start:end]
        print(f"  Query  {start+1:>4}  {q_seq}  {end}")
        print(f"              {m_seq}")
        print(f"  Sbjct  {start+1:>4}  {s_seq}  {end}")
        print()
    print("="*65)

# ─────────────────────────────────────────────
# VISUALIZATION
# ─────────────────────────────────────────────

def visualize_alignment(result, seq1, seq2, label):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle(f"Sequence Alignment — {result['type']}",
                fontsize=14, fontweight="bold")

    # Chart 1 — Scoring matrix heatmap
    ax1  = axes[0]
    matrix = result["matrix"]
    im   = ax1.imshow(matrix, cmap="YlOrRd", aspect="auto")
    ax1.set_xticks(range(len(seq2)+1))
    ax1.set_xticklabels(["-"] + list(seq2), fontsize=8)
    ax1.set_yticks(range(len(seq1)+1))
    ax1.set_yticklabels(["-"] + list(seq1), fontsize=8)
    ax1.set_title("Scoring Matrix", fontweight="bold")
    plt.colorbar(im, ax=ax1)

    # Annotate small matrices
    if len(seq1) <= 15 and len(seq2) <= 15:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                ax1.text(j, i, str(matrix[i][j]),
                        ha="center", va="center", fontsize=7, color="black")

    # Chart 2 — Alignment statistics
    ax2      = axes[1]
    labels   = ["Matches", "Mismatches", "Gaps"]
    values   = [result["matches_count"], result["mismatches"], result["gaps"]]
    colors   = ["#4CAF50", "#F44336", "#FF9800"]
    bars     = ax2.bar(labels, values, color=colors, edgecolor="black", width=0.5)
    for bar, val in zip(bars, values):
        ax2.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() + 0.1, str(val),
                ha="center", fontweight="bold")
    ax2.set_title("Alignment Statistics", fontweight="bold")
    ax2.set_ylabel("Count")

    # Add identity text
    ax2.text(0.5, 0.95, f"Identity: {result['identity']}%",
            transform=ax2.transAxes, ha="center",
            fontsize=12, fontweight="bold", color="#2196F3")

    plt.tight_layout()
    filename = f"alignment_{label.replace(' ', '_')}.png"
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n   Chart saved as: {filename}")

# ─────────────────────────────────────────────
# SAVE RESULTS
# ─────────────────────────────────────────────

def save_results(all_results, seq_pairs):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename  = f"alignment_results_{timestamp}.txt"
    with open(filename, "w") as f:
        f.write("SEQUENCE ALIGNMENT TOOL — RESULTS\n")
        f.write(f"Generated: {datetime.now().strftime('%d-%m-%Y %H:%M:%S')}\n")
        f.write("="*65 + "\n\n")
        for i, (result, (seq1, seq2)) in enumerate(zip(all_results, seq_pairs)):
            f.write(f"Pair {i+1}\n")
            f.write(f"Sequence 1     : {seq1}\n")
            f.write(f"Sequence 2     : {seq2}\n")
            f.write(f"Type           : {result['type']}\n")
            f.write(f"Score          : {result['score']}\n")
            f.write(f"Identity       : {result['identity']}%\n")
            f.write(f"Similarity     : {result['similarity']}%\n")
            f.write(f"Gaps           : {result['gaps']}\n")
            f.write(f"Aligned seq 1  : {result['aligned1']}\n")
            f.write(f"Match line     : {result['matches']}\n")
            f.write(f"Aligned seq 2  : {result['aligned2']}\n")
            f.write("-"*65 + "\n\n")
    print(f"\n   Results saved to: {filename}")

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

if __name__ == "__main__":
    # Install numpy if needed
    try:
        import numpy as np
    except ImportError:
        print("Installing numpy...")
        os.system("pip install numpy")
        import numpy as np

    print("\n" + "="*65)
    print("    Sequence Alignment Tool")
    print("    Author: Padma Shree | Phase 2 - Project 1")
    print("    Algorithms: Needleman-Wunsch | Smith-Waterman")
    print("="*65)

    all_results = []
    seq_pairs   = []
    pair_count  = 1

    while True:
        print(f"\n--- Sequence Pair #{pair_count} ---")
        print("Enter Sequence 1 (or 'done' to finish):")
        seq1 = input("Seq1 >>> ").strip().upper()
        if seq1 == "DONE":
            break

        print("Enter Sequence 2:")
        seq2 = input("Seq2 >>> ").strip().upper()

        print("\nAlignment type:")
        print("  1. Global (Needleman-Wunsch)")
        print("  2. Local  (Smith-Waterman)")
        print("  3. Both")
        choice = input("Choice (1/2/3): ").strip()

        label = f"Pair_{pair_count}"

        if choice in ("1", "3"):
            print("\n  Running Global Alignment...")
            result = needleman_wunsch(seq1, seq2)
            blast_output(result, seq1, seq2)
            visualize_alignment(result, seq1, seq2, f"Global_{label}")
            all_results.append(result)
            seq_pairs.append((seq1, seq2))

        if choice in ("2", "3"):
            print("\n  Running Local Alignment...")
            result = smith_waterman(seq1, seq2)
            blast_output(result, seq1, seq2)
            visualize_alignment(result, seq1, seq2, f"Local_{label}")
            all_results.append(result)
            seq_pairs.append((seq1, seq2))

        pair_count += 1

    if all_results:
        save_results(all_results, seq_pairs)
        print(f"\n  Total alignments performed: {len(all_results)}")
        print("  Great work, Paddu! 🧬")
    else:
        print("\n  No sequences entered.")
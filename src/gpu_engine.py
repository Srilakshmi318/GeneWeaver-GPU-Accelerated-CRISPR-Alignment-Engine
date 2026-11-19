
import numpy as np
from numba import cuda
import time


@cuda.jit
def crispr_kernel(genome, guide, max_mismatches,
                  positions, strands, mismatch_counts, penalties):

    idx = cuda.grid(1)

    n = genome.size
    m = guide.size
    total = n - m - 2

    if idx >= total:
        return

    # ---------------------------------
    # FORWARD STRAND: TARGET + NGG
    # ---------------------------------

    pam_start = idx + m

    if (genome[pam_start + 1] == 71 and
            genome[pam_start + 2] == 71):

        mismatches = 0
        penalty = 0

        for j in range(m):

            if genome[idx + j] != guide[j]:

                mismatches += 1

                if j >= m - 10:
                    penalty += 2
                else:
                    penalty += 1

        if mismatches <= max_mismatches:

            positions[idx] = idx
            strands[idx] = 1
            mismatch_counts[idx] = mismatches
            penalties[idx] = penalty

    # ---------------------------------
    # REVERSE STRAND: CCN + TARGET
    # ---------------------------------

    if (genome[idx] == 67 and
            genome[idx + 1] == 67):

        mismatches = 0
        penalty = 0

        for j in range(m):

            # Reverse complement comparison
            base = genome[idx + 3 + m - 1 - j]

            if base == 65:
                base = 84
            elif base == 84:
                base = 65
            elif base == 67:
                base = 71
            elif base == 71:
                base = 67

            if base != guide[j]:

                mismatches += 1

                if j >= m - 10:
                    penalty += 2
                else:
                    penalty += 1

        if mismatches <= max_mismatches:

            output_index = total + idx

            positions[output_index] = idx + 3
            strands[output_index] = -1
            mismatch_counts[output_index] = mismatches
            penalties[output_index] = penalty


def scan_crispr_gpu(genome, guide, max_mismatches=3):

    genome = genome.upper()
    guide = guide.upper()

    if not genome or not guide:
        raise ValueError("Genome and guide cannot be empty")

    if set(genome + guide) - set("ACGT"):
        raise ValueError("Only A, C, G and T are supported")

    n = len(genome)
    m = len(guide)

    if n < m + 3:
        return [], 0.0

    total = n - m - 2

    start_time = time.perf_counter()

    g = cuda.to_device(
        np.frombuffer(genome.encode("ascii"), dtype=np.uint8).copy()
    )

    q = cuda.to_device(
        np.frombuffer(guide.encode("ascii"), dtype=np.uint8).copy()
    )

    d_positions = cuda.to_device(
        np.full(2 * total, -1, dtype=np.int32)
    )

    d_strands = cuda.to_device(
        np.zeros(2 * total, dtype=np.int32)
    )

    d_mismatches = cuda.to_device(
        np.full(2 * total, -1, dtype=np.int32)
    )

    d_penalties = cuda.to_device(
        np.full(2 * total, -1, dtype=np.int32)
    )

    threads = 128
    blocks = (total + threads - 1) // threads

    crispr_kernel[blocks, threads](
        g, q, max_mismatches,
        d_positions, d_strands,
        d_mismatches, d_penalties
    )

    cuda.synchronize()

    positions = d_positions.copy_to_host()
    strands = d_strands.copy_to_host()
    mismatches = d_mismatches.copy_to_host()
    penalties = d_penalties.copy_to_host()

    elapsed = time.perf_counter() - start_time

    results = []

    for i in range(2 * total):

        if positions[i] >= 0:

            strand = "+" if strands[i] == 1 else "-"

            if strand == "+":
                target = genome[positions[i]:positions[i] + m]
                pam = genome[positions[i] + m:positions[i] + m + 3]

            else:
                target = "".join(
                    {"A": "T", "T": "A", "C": "G", "G": "C"}[b]
                    for b in reversed(
                        genome[positions[i]:positions[i] + m]
                    )
                )

                pam_reference = genome[positions[i] - 3:positions[i]]
                pam = "".join(
                    {"A": "T", "T": "A", "C": "G", "G": "C"}[b]
                    for b in reversed(pam_reference)
                )

            results.append({
                "position": int(positions[i]),
                "strand": strand,
                "target": target,
                "pam": pam,
                "mismatches": int(mismatches[i]),
                "penalty": int(penalties[i])
            })

    return results, elapsed


# ---------------------------------
# TEST GENEWEAVER CUDA ENGINE
# ---------------------------------

test_genome = "ATCGATCGATCGGATCCGATCGATCGATCG"
test_guide = "ATCGAT"

results, execution_time = scan_crispr_gpu(
    test_genome,
    test_guide,
    max_mismatches=2
)

print("\nGENEWEAVER - INTEGRATED CRISPR GPU ENGINE")
print("=" * 55)

print("GPU:", cuda.get_current_device().name)
print("Guide:", test_guide)

for result in results:

    print("\nPosition:", result["position"])
    print("Strand:", result["strand"])
    print("Target:", result["target"])
    print("PAM:", result["pam"])
    print("Mismatches:", result["mismatches"])
    print("Penalty:", result["penalty"])

print("\nTotal CRISPR candidates:", len(results))
print("Execution time:", execution_time, "seconds")
print("\nGPU CRISPR scanning completed!")

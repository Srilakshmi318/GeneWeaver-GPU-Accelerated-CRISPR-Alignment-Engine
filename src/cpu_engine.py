
import time


def find_matches(genome, guide, max_mismatches=2):

    genome = genome.upper()
    guide = guide.upper()

    if not genome or not guide:
        raise ValueError("Genome and guide cannot be empty")

    if max_mismatches < 0:
        raise ValueError("Maximum mismatches cannot be negative")

    if len(guide) > len(genome):
        return [], 0.0

    results = []

    start_time = time.perf_counter()

    for i in range(len(genome) - len(guide) + 1):

        target = genome[i:i + len(guide)]

        mismatches = [
            position
            for position, (a, b) in enumerate(zip(target, guide))
            if a != b
        ]

        if len(mismatches) <= max_mismatches:

            results.append({
                "position": i,
                "target": target,
                "mismatches": len(mismatches),
                "mismatch_positions": mismatches
            })

    elapsed_time = time.perf_counter() - start_time

    return results, elapsed_time


if __name__ == "__main__":

    genome = "ATCGATCGATCGGATCCGATCGATCGATCG"
    guide = "ATCGAT"

    matches, elapsed_time = find_matches(
        genome,
        guide,
        max_mismatches=1
    )

    print("\nGENEWEAVER - CPU ALIGNMENT ENGINE")
    print("=" * 45)

    print("Guide DNA:", guide)
    print("Maximum mismatches: 1")

    print("\nMatching locations:")

    for match in matches:
        print("\nPosition:", match["position"])
        print("Target:", match["target"])
        print("Mismatches:", match["mismatches"])
        print("Mismatch positions:", match["mismatch_positions"])

    print("\nTotal matches:", len(matches))
    print("Execution time:", elapsed_time, "seconds")
    print("\nCPU alignment completed!")


import time


def find_matches(genome, guide, max_mismatches=2):
    """
    Search DNA for guide sequences allowing mismatches.
    Returns matching positions and mismatch counts.
    """

    results = []
    guide_length = len(guide)

    start_time = time.perf_counter()

    for i in range(len(genome) - guide_length + 1):
        target = genome[i:i + guide_length]

        mismatches = sum(
            1 for a, b in zip(target, guide) if a != b
        )

        if mismatches <= max_mismatches:
            results.append({
                "position": i,
                "target": target,
                "mismatches": mismatches
            })

    elapsed_time = time.perf_counter() - start_time

    return results, elapsed_time


if __name__ == "__main__":

    genome = "ATCGATCGATCGGATCCGATCGATCGATCG"
    guide = "ATCGAT"

    matches, elapsed = find_matches(
        genome, guide, max_mismatches=1
    )

    print("\nGENEWEAVER - CPU ALIGNMENT ENGINE")
    print("-" * 40)

    for match in matches:
        print("Position:", match["position"])
        print("Target:", match["target"])
        print("Mismatches:", match["mismatches"])
        print("-" * 30)

    print("Total matches:", len(matches))
    print("Execution time:", round(elapsed, 6), "seconds")


def reverse_complement(sequence):
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    return "".join(
        complement[base] for base in reversed(sequence.upper())
    )


def has_ngg_pam(pam):
    pam = pam.upper()
    return len(pam) == 3 and pam[0] in "ACGT" and pam[1:] == "GG"


def score_target(guide, target):
    if len(guide) != len(target):
        raise ValueError("Guide and target lengths must match")

    mismatch_positions = []
    penalty = 0

    for position, (a, b) in enumerate(zip(guide, target)):
        if a != b:
            mismatch_positions.append(position)

            if position >= len(guide) - 10:
                penalty += 2
            else:
                penalty += 1

    return {
        "mismatches": len(mismatch_positions),
        "mismatch_positions": mismatch_positions,
        "penalty": penalty
    }


def find_pam_sites(genome, guide, max_mismatches=3):

    genome = genome.upper()
    guide = guide.upper()

    if not guide:
        raise ValueError("Guide cannot be empty")

    results = []
    guide_length = len(guide)

    for start in range(len(genome) - guide_length - 2):

        # Forward strand: target followed by NGG
        forward_target = genome[start:start + guide_length]
        forward_pam = genome[
            start + guide_length:start + guide_length + 3
        ]

        if has_ngg_pam(forward_pam):

            score = score_target(guide, forward_target)

            if score["mismatches"] <= max_mismatches:
                results.append({
                    "position": start,
                    "strand": "+",
                    "target": forward_target,
                    "pam": forward_pam,
                    **score
                })

        # Reverse strand: CCN appears before the target
        reverse_pam_reference = genome[start:start + 3]

        if (
            len(reverse_pam_reference) == 3
            and reverse_pam_reference[0:2] == "CC"
        ):

            reverse_target_reference = genome[
                start + 3:start + 3 + guide_length
            ]

            if len(reverse_target_reference) != guide_length:
                continue

            reverse_target = reverse_complement(
                reverse_target_reference
            )

            reverse_pam = reverse_complement(
                reverse_pam_reference
            )

            score = score_target(guide, reverse_target)

            if score["mismatches"] <= max_mismatches:
                results.append({
                    "position": start + 3,
                    "strand": "-",
                    "target": reverse_target,
                    "pam": reverse_pam,
                    **score
                })

    return results


if __name__ == "__main__":

    genome = "ATCGATCGATCGGATCCGATCGATCGATCG"
    guide = "ATCGAT"

    matches = find_pam_sites(genome, guide, max_mismatches=2)

    print("\nGENEWEAVER - TWO-STRAND PAM SCANNER")
    print("=" * 45)

    for match in matches:
        print("\nPosition:", match["position"])
        print("Strand:", match["strand"])
        print("Target:", match["target"])
        print("PAM:", match["pam"])
        print("Mismatches:", match["mismatches"])
        print("Mismatch positions:", match["mismatch_positions"])
        print("Penalty:", match["penalty"])

    print("\nTotal PAM-compatible matches:", len(matches))

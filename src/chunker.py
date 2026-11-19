
def split_sequence(sequence, chunk_size=10, overlap=3):
    """
    Split a DNA sequence into overlapping chunks.
    """

    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")

    if overlap < 0 or overlap >= chunk_size:
        raise ValueError("overlap must be between 0 and chunk_size - 1")

    chunks = []
    step = chunk_size - overlap

    for start in range(0, len(sequence), step):
        end = min(start + chunk_size, len(sequence))

        chunks.append({
            "start": start,
            "end": end,
            "sequence": sequence[start:end]
        })

        if end == len(sequence):
            break

    return chunks


if __name__ == "__main__":
    sample_sequence = "ATCGATCGATCGGATCCGATCGATCGATCG"

    chunks = split_sequence(
        sample_sequence,
        chunk_size=10,
        overlap=3
    )

    print("\nGENEWEAVER - DNA CHUNKER")
    print("=" * 40)

    for number, chunk in enumerate(chunks, start=1):
        print(f"\nChunk {number}")
        print("Start:", chunk["start"])
        print("End:", chunk["end"])
        print("DNA:", chunk["sequence"])

    print("\nTotal chunks:", len(chunks))

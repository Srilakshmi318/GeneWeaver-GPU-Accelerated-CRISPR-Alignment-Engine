
from pathlib import Path

VALID_BASES = set("ACGTRYSWKMBDHVN")


def parse_fasta(file_path):

    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(
            f"FASTA file not found: {file_path}"
        )

    records = []
    current_id = None
    current_sequence = []

    with open(file_path, "r") as file:

        for line in file:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if current_id is not None:
                    records.append(
                        create_record(current_id, current_sequence)
                    )

                current_id = line[1:].strip()

                if not current_id:
                    raise ValueError("Empty FASTA header")

                current_sequence = []

            else:

                if current_id is None:
                    raise ValueError(
                        "Sequence found before FASTA header"
                    )

                current_sequence.append(line.upper())

    if current_id is not None:
        records.append(
            create_record(current_id, current_sequence)
        )

    if not records:
        raise ValueError("No valid DNA sequences found")

    return records


def create_record(record_id, sequence_parts):

    sequence = "".join(sequence_parts)

    if not sequence:
        raise ValueError(f"Empty sequence: {record_id}")

    invalid = set(sequence) - VALID_BASES

    if invalid:
        raise ValueError(
            f"Invalid bases in {record_id}: {sorted(invalid)}"
        )

    return {
        "id": record_id,
        "sequence": sequence,
        "length": len(sequence)
    }


def display_genome_info(records):

    print("\n" + "=" * 45)
    print("       GENEWEAVER - GENOMIC PARSER")
    print("=" * 45)

    total_bases = 0

    for record in records:

        print("\nChromosome:", record["id"])
        print("Sequence length:", record["length"])
        print("DNA sequence:", record["sequence"])

        total_bases += record["length"]

    print("\nTotal chromosomes:", len(records))
    print("Total DNA bases:", total_bases)

    print("\nFASTA parsing completed successfully!")
    print("=" * 45)


if __name__ == "__main__":

    project_root = Path(__file__).resolve().parent.parent

    fasta_file = (
        project_root / "data" / "input" / "sample_genome.fasta"
    )

    try:

        records = parse_fasta(fasta_file)
        display_genome_info(records)

    except (FileNotFoundError, ValueError) as error:

        print("ERROR:", error)

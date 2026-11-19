
from dask import delayed, compute
import time


def process_chunk(chunk, guide, max_mismatches=2):
    from src.cpu_engine import find_matches

    results, elapsed = find_matches(
        chunk["sequence"],
        guide,
        max_mismatches
    )

    for result in results:
        result["position"] += chunk["start"]

    return results


def scan_distributed(chunks, guide, max_mismatches=2):
    start_time = time.perf_counter()

    tasks = [
        delayed(process_chunk)(chunk, guide, max_mismatches)
        for chunk in chunks
    ]

    chunk_results = compute(*tasks)

    results = []

    for chunk_result in chunk_results:
        results.extend(chunk_result)

    elapsed = time.perf_counter() - start_time

    return results, elapsed

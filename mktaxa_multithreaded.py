#! /usr/bin/env python
import os
import sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from pandas import DataFrame
from qiime2 import Artifact
from qiime2.plugins import feature_classifier


def run_blast(pid, reference_reads_path, reference_taxonomy_path, repseqs_path):
    # Each worker process loads its own artifacts.
    reference_reads = Artifact.load(reference_reads_path)
    reference_taxonomy = Artifact.load(reference_taxonomy_path)
    repseqs = Artifact.load(repseqs_path)

    result = feature_classifier.pipelines.classify_consensus_blast(
        query=repseqs,
        reference_reads=reference_reads,
        reference_taxonomy=reference_taxonomy,
        maxaccepts=1,
        perc_identity=float(pid),
        query_cov=0.4,
        num_threads=1,  # BLAST will be single-threaded in -subject mode
    ).classification.view(DataFrame)

    return pid, result


def main():
    if len(sys.argv) != 4:
        print("Usage: ./mktaxa.py ref_reads ref_taxa rep_seqs")
        sys.exit(1)

    reference_reads_path = sys.argv[1]
    reference_taxonomy_path = sys.argv[2]
    repseqs_path = sys.argv[3]

    confidence = np.linspace(1.0, 0.70, 80)

    # Start with 2. Increase to 3 or 4 only if RAM permits.
    workers = 2

    results = {}

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                run_blast,
                pid,
                reference_reads_path,
                reference_taxonomy_path,
                repseqs_path,
            ): pid
            for pid in confidence
        }

        for future in as_completed(futures):
            pid = futures[future]

            try:
                pid, result = future.result()
                results[pid] = result
                print(f"Finished: {pid:.4f}")
            except Exception as e:
                print(f"Error at {pid:.4f}: {e}")
                results[pid] = None

    # Restore original PID ordering.
    taxonomy = [results.get(pid) for pid in confidence]

    if taxonomy[0] is None:
        raise RuntimeError("The first BLAST run failed.")

    combined = taxonomy[0].copy()
    combined["Confidence"] = combined["Consensus"].astype(float)

    for conf, result in zip(confidence[1:], taxonomy[1:]):
        if result is None:
            print(f"Skipping confidence {conf:.4f} due to missing result.")
            continue

        result = result.copy()
        result["Confidence"] = result["Consensus"].astype(float) * conf

        unassigned = combined["Taxon"] == "Unassigned"
        combined.loc[unassigned, "Taxon"] = result.loc[unassigned, "Taxon"]
        combined.loc[unassigned, "Confidence"] = result.loc[
            unassigned, "Confidence"
        ]

    Artifact.import_data(
        "FeatureData[Taxonomy]", combined
    ).save("superblast_taxonomy")


if __name__ == "__main__":
    main()
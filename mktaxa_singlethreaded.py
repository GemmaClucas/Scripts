#! /usr/bin/env python
import sys
from pandas import DataFrame
import numpy as np

if len(sys.argv) != 4:
    print("Usage: ./mktaxa.py ref_reads ref_taxa rep_seqs")
    sys.exit()

from qiime2 import Artifact
from qiime2.plugins import feature_classifier

reference_reads = Artifact.load(sys.argv[1])
reference_taxonomy = Artifact.load(sys.argv[2])
repseqs = Artifact.load(sys.argv[3])

def run_blast(pid):
    result = feature_classifier.pipelines.classify_consensus_blast(
        query=repseqs,
        reference_reads=reference_reads,
        reference_taxonomy=reference_taxonomy,
        maxaccepts=1,
        perc_identity=pid,
        query_cov=0.4
    ).classification.view(DataFrame)
    print(f"Result at confidence {pid}:")
    print(result.head())
    print(result.columns)
    return result

confidence = np.linspace(1.0, 0.70, 80)

taxonomy = []
for pid in confidence:
    try:
        print(f"Processing confidence level: {pid}")
        result = run_blast(pid)
        taxonomy.append(result)
    except Exception as e:
        print(f"Error processing confidence {pid}: {e}")
        taxonomy.append(None)

combined = taxonomy[0]
combined["Confidence"] = combined["Consensus"].astype(float)  # Initialize Confidence

for i in range(1, len(confidence)):
    conf = confidence[i]
    result = taxonomy[i]
    if result is None:
        print(f"Skipping confidence {conf} due to missing result.")
        continue

    result["Confidence"] = result["Consensus"].astype(float) * conf
    for ind in combined.index:
        if combined.loc[ind, "Taxon"] == "Unassigned":
            combined.loc[ind, "Taxon"] = result.Taxon[ind]
            combined.loc[ind, "Confidence"] = result.Confidence[ind]

try:
    Artifact.import_data("FeatureData[Taxonomy]", combined).save("superblast_taxonomy")
except Exception as e:
    print(f"Failed to save taxonomy artifact: {e}")

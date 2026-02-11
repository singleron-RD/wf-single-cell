"""Process results for reporting for a single sample."""
import json
from pathlib import Path

from workflow_glue.models.common import SampleType
from workflow_glue.models.custom import (
    SCResultsContents, SCSample, UmapResults)
from workflow_glue.sc_util import StatusAggregator
from workflow_glue.util import get_named_logger, wf_parser


def parse_metadata(meta_path):
    """Get metadata info."""
    with open(meta_path) as fh:
        meta = json.load(fh)
        sample_meta = {
            "n_seqs": meta["n_seqs"],
            "type": meta["type"],
            "alias": meta["alias"],
            "barcode": meta["barcode"]
        }
    return sample_meta


def main(args):
    """Run the entry point."""
    logger = get_named_logger("process_results")

    sample_meta = parse_metadata(args.metadata)
    alias = sample_meta['alias']
    status = StatusAggregator.load(args.statuses)

    results = SCResultsContents(
        umap=UmapResults(
            gene_status=status['gene_umap_status'],
            gene_umap_file_names=status.get('gene_umap_replicate_path'),
            transcript_status=status['transcript_umap_status'],
            transcript_umap_file_names=status.get('transcript_umap_replicate_path')
        )
    )
    sample = SCSample(
        alias=alias,
        barcode=sample_meta['barcode'],
        sample_type=SampleType(sample_meta['type']),
        sample_pass=None,  # init None, validation/QC wfs will fill
        sample_checks=[],  # init empty list, validation/QC wfs will fill
        results=results,
        config=None
    )
    sample.to_json(args.results_json)
    logger.info(f"Results written to {args.results_json}")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("process_results")
    parser.add_argument("results_json", help="Workflow results output JSON")
    parser.add_argument("metadata", type=Path, help="Sample metadata JSON.")
    parser.add_argument("statuses", type=Path, help="Dir of JSON status files.")

    return parser

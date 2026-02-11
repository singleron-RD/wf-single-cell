"""Aggregate per-sample results to get WorkflowResult model for reporting."""
from enum import Enum
import json
from pathlib import Path

import dacite
from workflow_glue.models.custom import (
    SCSample, SCWorkflowResult
)

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("combine_results")

    samples = []
    for sample_input in args.inputs.iterdir():
        with open(sample_input) as fh:
            data = json.load(fh)
            sample_obj = dacite.from_dict(
                data_class=SCSample,
                data=data,
                config=dacite.Config(cast=[Enum])
            )
            samples.append(sample_obj)
    workflow = SCWorkflowResult(samples=samples)
    workflow.to_json(args.results_output)
    logger.info(f"Results written to {args.results_output}")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("combine_results")
    parser.add_argument(
        "inputs", type=Path, help="Workflow results output JSONs directory")
    parser.add_argument(
        "results_output", type=Path, help="Combined results output JSON")

    return parser

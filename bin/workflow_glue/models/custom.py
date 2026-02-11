"""Extended models for the workflow."""
from dataclasses import dataclass, field

from workflow_glue.models.common import (
    ResultsContents, Sample, WorkflowBaseModel, WorkflowResult)


@dataclass
class UmapResults(WorkflowBaseModel):
    """A place to store UMAP results."""

    gene_status: str = field(
        metadata={
            "title": "Gene UMAP status",
            "description": "'OK' on successful UMAP generation or an error message"
        })

    transcript_status: str = field(
        metadata={
            "title": "Transcript UMAP status",
            "description": "'OK' on successful UMAP generation or an error message"
        })

    gene_umap_file_names: list[str] | None = field(
        metadata={
            "title": "Gene UMAP file names",
            "description": (
                "File names of successfully generated gene expression UMAP "
                "matrices.")
        })

    transcript_umap_file_names: list[str] | None = field(
        metadata={
            "title": "Transcript UMAP file names",
            "description": (
                "File names of successfully generated transcript expression UMAP",
                "matrices.")
        })


@dataclass
class SCResultsContents(ResultsContents):
    """Results for a sample."""

    umap: UmapResults = field(
        metadata={
            "title": "UMAP status",
            "description": "Status of UMAP generation"
        })


@dataclass
class SCSample(Sample):
    """A sample sheet entry and its corresponding checks and related results."""

    results: SCResultsContents = field(
        default=None, metadata={
            "title": "Sample results",
            "description": "Further specific workflow results for this sample"})


@dataclass
class SCWorkflowResult(WorkflowResult):
    """
    Definition for results that will be returned by this workflow.

    This structure will be passed through by Gizmo speaking clients as
    WorkflowInstance.results.
    """

    samples: list[SCSample] = field(
        metadata={
            "title": "Samples",
            "description": "Samples in this workflow instance"})

    def __post_init__(self):
        """Determine overall status for the workflow given the sample results."""
        self.workflow_pass = all(
            sample.sample_pass for sample in self.samples)

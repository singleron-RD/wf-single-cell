"""Test process matrix.py."""

import os

import numpy as np
import pandas as pd
import pytest
from workflow_glue.process_matrix import ExpressionMatrix, make_umaps
from workflow_glue.sc_util import StatusRecorder


"""Pytest fixtures and tests for process_matrix.

Each case describes:
  n_genes x n_cells: shape of the expression matrix
  [count_low, count_high]: bounds for random uniform value generation
  fail: whether UMAP generation is expected to fail and return a status error
"""
CASES = [
    # Passing cases
    dict(n_genes=2, n_cells=8, count_low=0, count_high=100, fail=False),
    dict(n_genes=2, n_cells=4, count_low=0, count_high=100, fail=False),
    dict(n_genes=1, n_cells=4, count_low=0, count_high=100, fail=False),
    # 0 replicates does not make sense, but should be handled
    dict(n_genes=1, n_cells=4, count_low=0, count_high=100, reps=0, fail=False),
    # A matrix of all zeros should still pass
    dict(n_genes=2, n_cells=8, count_low=0, count_high=0, fail=False),

    # Failing cases, too few cells
    dict(n_genes=2, n_cells=2, count_low=0, count_high=100, fail=True),
    dict(n_genes=400, n_cells=3, count_low=0, count_high=100, fail=True),
    dict(n_genes=5000, n_cells=3, count_low=10, count_high=100, fail=True),
    dict(n_genes=400, n_cells=2, count_low=0, count_high=100, fail=True),
]


def _case_id(c):
    """Build a compact id from case dict fields.

    Format: {genes}g_{cells}c_{low}-{high}
    Example: 2g_4c_0-100
    """
    return f"{c['n_genes']}genes_{c['n_cells']}cells_{c['count_low']}-{c['count_high']}"


TEST_CASES = [pytest.param(c, id=_case_id(c)) for c in CASES]


@pytest.mark.parametrize("params", TEST_CASES)
def test_make_umaps(tmp_path, params):
    """Test the UMAP generation with varying matrix sizes.

    This test just checks that the UMAPs are generated without error and that the
    correct status is returned due to low cell count.
    """
    rng = np.random.RandomState(2)

    cells = [f"cell{i+1}" for i in range(params["n_cells"])]
    features = [f"feature{j+1}" for j in range(params["n_genes"])]

    mat = rng.uniform(
        params["count_low"], params["count_high"],
        size=(params["n_genes"], params["n_cells"]))

    mtx = ExpressionMatrix(
        matrix=mat,
        features=np.array(features, dtype=bytes),
        cells=np.array(cells, dtype=bytes),
        sparse=False)
    umap_tsv = tmp_path / "umap.tsv"
    # desired principal components, will be set to  min(pcn, *matrix._matrix.shape)
    # in UMAP code
    pcn = 100

    replicates = 3
    if "reps" in params:
        replicates = params["reps"]

    max_umap_cells = 100

    os.environ["NUMBA_CACHE_DIR"] = str(tmp_path / "numba_cache")

    status_recorder = StatusRecorder(sample='TEST', path=tmp_path / "status.json")

    make_umaps(
        mtx,
        pcn,
        replicates,
        max_umap_cells,
        n_neighbors=15,
        min_dist=0.1,
        umap_tsv=umap_tsv,
        feature='gene',
        status_logger=status_recorder
    )

    returned_status = status_recorder.data['gene_umap_status']

    if params['fail'] is True:
        assert returned_status == (
            "gene UMAP cannot be generated "
            " with matrix containing fewer than 4 cells.")
    else:
        assert returned_status == 'OK'
        # Check the dimensions of the UMAPs

        if replicates > 0:
            umap_paths = status_recorder.data['gene_umap_replicate_path']
            for path in umap_paths:
                umap = pd.read_csv(
                    path, sep="\t", index_col=0)
                assert umap.shape == (params["n_cells"], 2)

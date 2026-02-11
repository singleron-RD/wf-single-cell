"""Pytests conftest."""
import os
import pathlib

import pytest


@pytest.fixture(autouse=True)
def ensure_cwd():
    """Ensure the current working directory is set to the test directory."""
    try:
        os.getcwd()
    except FileNotFoundError:
        # fallback into the test directory
        os.chdir(str(pathlib.Path(__file__).parent))
    yield

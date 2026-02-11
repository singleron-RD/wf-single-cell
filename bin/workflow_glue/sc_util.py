"""Common code to be used across workflow scripts."""
import collections
from dataclasses import dataclass, field
import json

kit_adapters = {
    '3prime': {
        'adapter1': 'CTACACGACGCTCTTCCGATCT',
        'adapter2': 'ATGTACTCTGCGTTGATACCACTGCTT'
    },
    'multiome': {
        'adapter1': 'CTACACGACGCTCTTCCGATCT',
        'adapter2': 'ATGTACTCTGCGTTGATACCACTGCTT'
    },
    'visium': {
        'adapter1': 'CTACACGACGCTCTTCCGATCT',
        'adapter2': 'ATGTACTCTGCGTTGATACCACTGCTT'
    },
    '5prime': {
        'adapter1': 'CTACACGACGCTCTTCCGATCT',
        'adapter2': 'GTACTCTGCGTTGATACCACTGCTT'
    }
}

revcomp_map = str.maketrans("ACGTacgt", "TGCAtgca")


def rev_cmp(seq):
    """Reverse complement a DNA sequence."""
    return seq[::-1].translate(revcomp_map)


class StatsSummary(collections.Counter):
    """Summary dictionary for storing."""

    fields = {}  # subclasses should fill this in

    def __init__(self, *args, **kwargs):
        """Count some numbers."""
        self.update(*args, **kwargs)

    @classmethod
    def from_pandas(cls, df):
        """Create an instance from a pandas dataframe."""
        raise NotImplementedError("This method has not been implemented.")

    def to_dict(self):
        """Create dictionary with explicit zeroes."""
        return {k: self[k] for k in self}

    @classmethod
    def from_json(cls, fname):
        """Create and instance from a JSON file."""
        with open(fname, "r") as fh:
            data = json.load(fh)
        return cls(data)

    def to_json(self, fname):
        """Save to JSON."""
        with open(fname, "w") as fh:
            json.dump(self.to_dict(), fh, indent=4)


@dataclass
class StatusRecorder:
    """Save process status information to a JSON file."""

    sample: str  # noqa NT001
    path: str # noqa NT001
    data: dict = field(default_factory=dict)  # noqa NT001
    # Determine whether a status key is expected ot be a single value of several values
    schema = {
        "single_values": [
            "transcript_umap_status", "gene_umap_status"],
        "multiple_values": [
            "transcript_umap_replicate_path", "gene_umap_replicate_path"]
    }

    def add(self, key, value):
        """Set a status key to a value."""
        if key in self.schema['single_values']:
            self.data[key] = value
        elif key in self.schema['multiple_values']:
            self.data.setdefault(key, []).append(value)
        else:
            raise KeyError(f"Unknown status key: {key}")

    def write_json(self):
        """Write status to file."""
        with open(self.path, "w") as fh:
            json.dump({self.sample: self.data}, fh, indent=2)


class StatusAggregator(dict):
    """Aggregate workflow status information from multiple JSON files.

    Each JSON file comes from a different step in the workflow and all
    derive from a single sample.
    """

    @classmethod
    def load(cls, path):
        """Load status data from a directory of JSON files."""
        data = {}
        for p in path.iterdir():
            with open(p) as fh:
                for _, kv in json.load(fh).items():
                    data.update(kv)
        return cls(data)

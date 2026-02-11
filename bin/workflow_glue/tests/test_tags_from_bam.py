"""Test tag_bam.py"."""
from unittest.mock import Mock

import pandas as pd
import pysam
from workflow_glue.tags_from_bam import main


def create_test_bam(dir_):
    """Create a BAM file with CB and other tags."""
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 1575, 'SN': 'chr1'}]
    }
    #  Visium HD barcodes
    barcodes = [
        's_002_um_00100_00100-1',   # BC 1 2um bin 00025,00025
        's_002_um_00103_00103-1',   # BC 1 2um bin 00025,00025
        's_002_um_00104_00104-1',   # BC 2 2um bin 00026,00026
    ]
    unsorted_bam = dir_ / 'unsorted.bam'
    sorted_bam = dir_ / 'sorted.bam'

    with pysam.AlignmentFile(unsorted_bam, "wb", header=header) as outf:
        for i, barcode in enumerate(barcodes):
            a = pysam.AlignedSegment()
            a.query_name = f"read{i}"
            a.query_sequence = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 0
            a.mapping_quality = 60
            a.cigarstring = "M36"
            a.set_tag("CB", barcode)
            a.set_tag('CY', '????')
            a.set_tag('UY', '????')
            a.set_tag('CR', 'AAAT')
            outf.write(a)
    pysam.sort("-o", str(sorted_bam), str(unsorted_bam))
    pysam.index(str(sorted_bam))
    return sorted_bam


def test_main(tmp_path):
    """Test the tags files and barcode counts."""
    bam = create_test_bam(tmp_path)
    tags_out = tmp_path / 'test_tags.tsv'
    barcode_couts_out = tmp_path / 'test_barcode_counts.tsv'

    args = Mock()
    args.bam_in = bam
    args.tags_out = tags_out
    args.barcode_counts_out = barcode_couts_out
    args.tags = 'CR CB CY UR UB UY'.split()
    args.chrom = 'chr1'

    main(args)

    bc_counts_results = pd.read_csv(
        barcode_couts_out, sep='\t', index_col=0)  # barcode index
    # Barcode counts are binned to 8um for use in the report
    # 10x visium HD tags have a '-1' suffix. The resulting binned tags do not,
    # The naming of the tags is not important as they are intermediate data used
    # for genreating rank plots.
    assert len(bc_counts_results) == 2
    assert bc_counts_results.loc['s_008_um_00025_00025', 'count'] == 2
    assert bc_counts_results.loc['s_008_um_00026_00026', 'count'] == 1

    # Barcodes in tags file are un-binned, and should be unchanged
    tags_results = pd.read_csv(tags_out, sep='\t', index_col=0)  # read_id index
    assert tags_results.loc['read0', 'CB'] == 's_002_um_00100_00100-1'
    assert tags_results.loc['read1', 'CB'] == 's_002_um_00103_00103-1'
    assert tags_results.loc['read2', 'CB'] == 's_002_um_00104_00104-1'

from argweavers.io import read_bed, subset_bed


def test_read_bed(bedfile):
    df = read_bed(bedfile)
    assert df.shape == (1828, 5)


def test_subset_bed(bedfile):
    df = read_bed(bedfile)
    subset = subset_bed(df, "chr:1000-2000")
    assert subset.shape == (24, 5)

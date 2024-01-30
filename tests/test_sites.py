from argweaver import read_sites


def test_read_sites():
    sites = read_sites("examples/sim1/sim1.sites")
    assert sites.chrom == "chr"
    assert sites.start == 1
    assert sites.end == 100000
    assert sites.data.shape == (170, 9)

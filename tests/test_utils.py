import pytest

from argweavers.utils import parse_newick, parse_nhx, parse_region


def test_parse_region():
    assert parse_region("chr:1000-2000") == ("chr", 1000, 2000)
    with pytest.raises(ValueError):
        parse_region("not valid")
    with pytest.raises(ValueError):
        parse_region("chr:2000-1000")


def test_parse_nhx():
    assert parse_nhx("&&NHX:age=0.1") == {"age": 0.1}


def test_parse_newick():
    tree = parse_newick(
        "((6:1002.805899[&&NHX:age=0.000000],4:1002.805899[&&NHX:age=0.000000])12:17041.819563[&&NHX:age=1002.805899],(2:18044.625462[&&NHX:age=0.000000],((5:5363.846221[&&NHX:age=0.000000],(7:122.586947[&&NHX:age=0.000000],1:122.586947[&&NHX:age=0.000000])14:5241.259274[&&NHX:age=122.586947])10:6697.962294[&&NHX:age=5363.846221],(3:1545.314509[&&NHX:age=0.000000],0:1545.314509[&&NHX:age=0.000000])11:10516.494006[&&NHX:age=1545.314509])9:5982.816948[&&NHX:age=12061.808515])13:0.000000[&&NHX:age=18044.625462])8[&&NHX:age=18044.625462];"
    )
    assert tree

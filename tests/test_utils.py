import pytest

from argweavers.utils import parse_region


def test_parse_region():
    assert parse_region("chr:1000-2000") == ("chr", 1000, 2000)
    with pytest.raises(ValueError):
        parse_region("not valid")
    with pytest.raises(ValueError):
        parse_region("chr:2000-1000")

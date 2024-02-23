import pytest

from argweavers.bin import require_executable


def test_require_executable():
    require_executable("python3")
    with pytest.raises(FileNotFoundError):
        require_executable("__NOT_FOUND__")

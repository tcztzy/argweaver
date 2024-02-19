import pathlib
import shutil

import pytest

from argweaver.scripts.smc2bed_all import main as smc2bed_all

examples_dir = pathlib.Path(__file__).parent.parent / "examples"

smc_files = list((examples_dir / "sim1" / "sim1.sample").glob("out.*.smc.gz"))
log_file = examples_dir / "sim1" / "sim1.sample" / "out.log"


@pytest.fixture
def sim1_sample(tmp_path):
    for smc_file in smc_files:
        shutil.copy(smc_file, tmp_path / smc_file.name)
    shutil.copy(log_file, tmp_path / log_file.name)
    return tmp_path


@pytest.fixture
def bedfile(sim1_sample):
    smc2bed_all([str(sim1_sample / "out")])
    return sim1_sample / "out.bed.gz"

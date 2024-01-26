import hashlib
import pathlib
import shutil

import pytest

from argweaver.bin import smc2bed
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


def sha256sum(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def test_smc2bed(sim1_sample):
    p = smc2bed(
        ["--sample", "0", sim1_sample / "out.0.smc.gz"],
        capture_output=True,
        return_process=True,
    )
    sha256 = "8fbf2326c4ed37e7232c77b882e41964d75d617031cd81021d842c9cbe33eaf3"
    assert hashlib.sha256(p.stdout).hexdigest() == sha256


def test_smc2bed_all(sim1_sample):
    sha256 = "577951a303869a2ab3aaf99e500ff29ce5613c27ac1a7cbef7dbfa63524136c3"
    smc2bed_all([str(sim1_sample / "out")])
    assert (sim1_sample / "out.bed.gz").exists() and (
        sim1_sample / "out.bed.gz.tbi"
    ).exists()
    print((sim1_sample / "out.bed.gz").stat().st_size)
    assert sha256sum(sim1_sample / "out.bed.gz") == sha256

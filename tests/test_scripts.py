import hashlib

from argweaver.scripts.smc2bed import main as smc2bed


def sha256sum(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def test_smc2bed(sim1_sample):
    out = smc2bed(
        ["--sample", "0", sim1_sample / "out.0.smc.gz"],
        capture_output=True,
    )
    sha256 = "8fbf2326c4ed37e7232c77b882e41964d75d617031cd81021d842c9cbe33eaf3"
    assert hashlib.sha256(out).hexdigest() == sha256


def test_smc2bed_all(bedfile):
    sha256 = "577951a303869a2ab3aaf99e500ff29ce5613c27ac1a7cbef7dbfa63524136c3"
    assert sha256sum(bedfile) == sha256

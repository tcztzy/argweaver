import hashlib

from argweavers.scripts.smc2bed import main as smc2bed


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
    sha256 = "94012cc3ecadddaf19f823ec361a72ac9b8073fe33dc75c1f731498c789d7023"
    assert sha256sum(bedfile) == sha256

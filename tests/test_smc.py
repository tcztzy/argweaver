from argweavers.smc import SMC


def test_smc(sim1_sample):
    SMC.from_path(sim1_sample / "out.0.smc.gz")

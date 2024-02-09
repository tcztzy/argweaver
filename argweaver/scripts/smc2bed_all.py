import argparse
import os
import sys
from itertools import count

from argweaver.bin import require_executable
from argweaver.scripts.smc2bed import main as smc2bed


class MetavarFormatter(argparse.HelpFormatter):
    def _get_default_metavar_for_optional(self, action):
        return f"<{action.dest}>"

    def _get_default_metavar_for_positional(self, action):
        return f"<{action.dest}>"


def main(argv=None):
    sort_bed = require_executable("sort-bed", "It is included in `bedops`")
    bgzip = require_executable("bgzip", "It is included in `samtools`")
    tabix = require_executable("tabix", "It is included in `samtools`")
    parser = argparse.ArgumentParser(
        description="Convert all *.smc files from a single arg-sample run into "
        "a bed.gz file which can be parsed by arg-summarize",
        formatter_class=MetavarFormatter,
    )
    parser.add_argument(
        "-s",
        dest="startnum",
        type=int,
        default=0,
        help="start with this MCMC rep of arg-sample (default: 0)",
    )
    parser.add_argument(
        "-e",
        dest="endnum",
        type=int,
        default=-1,
        help="end with this MCMC rep (default: stop when no file exists)",
    )
    parser.add_argument(
        "-i",
        dest="interval",
        type=int,
        default=0,
        help="sampling interval (will be auto-detected from output if not "
        "specified)",
    )
    parser.add_argument(
        "-r",
        dest="region",
        help="region (in format START-END, 1-based, inclusive) to pass to "
        "smc2bed (default: run on all coordinates)",
    )
    parser.add_argument(
        "baseout",
        help="base name of arg-sample output (specified with arg-sample -o). "
        "Creates a file called <baseout>.bed.gz.",
    )
    args = parser.parse_args(argv)

    baseout = args.baseout
    startnum = args.startnum
    interval = args.interval
    endnum = args.endnum
    region = args.region

    startfile = f"{baseout}.{startnum}.smc.gz"
    if not os.path.exists(startfile):
        sys.stderr.write(f"ERROR: {startfile} does not exist\n")
        sys.exit(1)

    if interval == 0:
        for interval in range(1, 1001):
            nextnum = startnum + interval
            if os.path.exists(f"{baseout}.{nextnum}.smc.gz"):
                break
        else:
            sys.stderr.write(
                "ERROR detecting sampling interval; try specifying with -i\n"
            )
            sys.exit(1)
    sys.stderr.write(f"starting at rep startnum={startnum}\n")
    sys.stderr.write(f"using sampling interval={interval}\n")

    num = startnum
    regionarg = ["--region", region] if region else []
    out = b""
    for num in count(startnum, interval):
        file = f"{baseout}.{num}.smc.gz"
        if (endnum != -1 and num > endnum) or not os.path.exists(file):
            num -= interval
            sys.stderr.write(f"ended at sample={num}\n")
            break
        sys.stderr.write(f"{num} {file}\n")
        print(os.stat(file).st_size)
        out += smc2bed(
            ["--sample", str(num), *regionarg, file],
            capture_output=True,
        )
        print(len(out))
        num += interval
    out = sort_bed(["-"], input=out, capture_output=True).stdout
    with open(f"{baseout}.bed.gz", "wb") as f:
        bgzip(["-"], input=out, stdout=f)
    tabix(["-p", "bed", f"{baseout}.bed.gz"])

    sys.stderr.write(f"wrote and indexed {baseout}.bed.gz\n")


if __name__ == "__main__":
    main()

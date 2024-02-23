import io
import sys
from contextlib import redirect_stdout

from argweavers import smc2bed


def main(args=None, capture_output=False):
    if args is not None:
        args = [str(arg) for arg in args]
        args = ["smc2bed"] + args
    if capture_output:
        buf = io.TextIOWrapper(io.BytesIO())
    else:
        buf = sys.stdout
    with redirect_stdout(buf):
        smc2bed(args)
    if capture_output:
        buf.seek(0)
        return buf.buffer.read()

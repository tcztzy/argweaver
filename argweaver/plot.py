from rpy2.robjects.packages import importr

argweaver = importr("argweaver")


def plot_tree(bedfile):
    argweaver.plotTreesFromBed(bedfile)

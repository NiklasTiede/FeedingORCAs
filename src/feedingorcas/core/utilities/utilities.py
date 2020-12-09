
from matplotlib.text import Annotation
from mpl_toolkits.mplot3d.proj3d import proj_transform
from rdkit import Chem


class Annotation3D(Annotation):
    """ Annotate the point xyz with text s """
    def __init__(self, s, xyz, *args, **kwargs):
        Annotation.__init__(self, s, xy=(0, 0), *args, **kwargs)
        self._verts3d = xyz

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.xy = (xs, ys)
        Annotation.draw(self, renderer)


def annotate3d(ax, s, *args, **kwargs):
    """ add annotation text s to to Axes3d ax """

    tag = Annotation3D(s, *args, **kwargs)
    ax.add_artist(tag)


def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts(
        "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


# This timer function converts elapsed seconds into a more under understandable format

# # copy this:
# from src.my_timer import time_converter
# import time
# start_time1 = time.perf_counter()
#
# print('code to be executed')
#
# print(f"\nScript execution finished after {time_converter(round(time.perf_counter()-start_time1, 3))}")


def time_converter(seconds):
    """Converts time (in seconds) in more understandable format.
    Example: time_conversion(131623.456) --> Script execution
    finished after 1 day, 12 hours and 33 minutes."""
    total_seconds = seconds
    days = seconds // (24 * 3600)
    seconds = seconds % (24 * 3600)
    hours = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    microseconds = (seconds % 1) * 1000
    days, hours, minutes, seconds, microseconds = (
        int(days),
        int(hours),
        int(minutes),
        int(seconds),
        int(microseconds),
    )
    if total_seconds < 1:
        return f"{microseconds} microsecond(s)"
    elif 1 <= total_seconds < 60:
        return f"{seconds} second(s) and {microseconds} microsecond(s)"
    elif 60 <= total_seconds < 3600:
        return f"{minutes} minute(s), {seconds} second(s) and {microseconds} microsecond(s)"
    elif 3600 <= total_seconds < 86400:
        return f"{hours} hour(s), {minutes} minute(s) and {seconds} second(s)"
    elif 86400 <= total_seconds:
        return f"{days} day(s), {hours} hour(s) and {minutes} minute(s)"

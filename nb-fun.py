import csv
import os
from IPython.display import display, Image
from sage.matrix.constructor import Matrix
from sage.misc.temporary_file import tmp_dir, graphics_filename


# TODO: find better way than defining a global
dirname = tmp_dir()


def nb_show(graphics_object,ext='png'):
    name = os.path.join(dirname,
                sage.misc.temporary_file.graphics_filename(ext))
    graphics_object.save(name)
    display(Image(name))

def matrix_of_csv(fname):
    with open(fname, 'rb') as csvfile:
        mreader = csv.reader(csvfile, delimiter=',', quotechar="'")
        return Matrix([[int(x) for x in row] for row in mreader])

def mtx_to_dict(m, keys=None):
    d = {}
    if keys == None or len(keys) != m.nrows():
        keys = range(m.nrows())

    for i in range(m.nrows()):
        d[keys[i]] = []
        for j in range(m.ncols()):
            if m[i][j] != 0:
                d[keys[i]].append(keys[j])

    return(d)

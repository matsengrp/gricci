import csv
from IPython.display import display, Image
from sage.matrix.constructor import MatrixFactory

# TODO: find better way
dirname=sage.misc.temporary_file.tmp_dir()

def nb_show(graphics_object,ext='png'):
    name = os.path.join(dirname,
                sage.misc.temporary_file.graphics_filename(ext))
    graphics_object.save(name)
    display(Image(name))

#def arr_arr_of_csv(fname):
def matrix_of_csv(fname):
    with open(fname, 'rb') as csvfile:
        mreader = csv.reader(csvfile, delimiter=',', quotechar="'")
        return Matrix([[int(x) for x in row] for row in mreader])

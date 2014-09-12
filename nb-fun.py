from IPython.display import display, Image

# TODO: find better way
dirname=sage.misc.temporary_file.tmp_dir()

def nb_show(graphics_object,ext='png'):
    name = os.path.join(dirname,
                sage.misc.temporary_file.graphics_filename(ext))
    graphics_object.save(name)
    display(Image(name))


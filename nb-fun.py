from IPython.display import display, Image

def random_file_path():
    return "_tmp/"+os.urandom(16).encode('hex')

def nb_show(graphics_object,suffix=".png"):
    name = random_file_path()+suffix
    graphics_object.save(name)
    display(Image(name))


# vim:ft=python

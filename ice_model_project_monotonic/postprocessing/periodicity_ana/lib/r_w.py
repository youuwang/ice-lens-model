def read_file(fname):
    try:
        fobj = open(fname, 'r')
        content = fobj.read()
        fobj.close()
        del fobj # release file handle
        return content
    except IOError as e:
        print(str(e))
        raise RuntimeError("Error reading/opening file '{}'".format(fname))

def write_file(fname, contents):
    try:
        fobj = open(fname, 'w')
        fobj.write(contents)
        fobj.close()
        del fobj # release file handle
    except IOError as e:
        print(str(e))
        raise RuntimeError("Error opening file '{}' for writing".format(fname))
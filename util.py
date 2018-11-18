from itertools import starmap


def format_table(names, align, formats, table, sep='   '):
    header = ['# ' + names[0]] + list(names[1:])
    cells = [header] + [
        list(starmap(format, zip(items, formats)))
        for items in table
    ]
    widths = list(map(max, zip(*[
        map(len, items) for items in cells
    ])))
    return '\n'.join([
        sep.join(starmap(adjust, zip(items, align, widths)))
        for items in cells
    ])


_adjust = {
    'l': str.ljust,
    'r': str.rjust,
}


def adjust(s, a, w):
    return _adjust[a](s, w)


def format_strengths(data):
    return ''.join([
        '{} = {!r};\n'.format(k, v)
        for k, v in data.items()
    ])

import re
import operator
import itertools
import collections
import numpy as np
import pandas as pd

def multimode(data):
    counts = collections.Counter(iter(data)).most_common()
    _, mode_items = next(itertools.groupby(counts, key=operator.itemgetter(1)), (0, []))
    return list(map(operator.itemgetter(0), mode_items))

def parse_genomic_location(loc):
    '''
    Parses an eccDNA genomic location string into coordinates.
    return: (chromosome, start, end)
    '''
    pile = re.compile(r':|-|\(|\)')
    tmp = pile.split(loc)
    if not tmp[3]:
        tmp[3] = '-'
    tmp = list(filter(None, tmp))
    tmp[1:3] = list(map(int, tmp[1:3]))
    return tuple(tmp)

def fragment_priority_adjuster(string):
    '''
     Adjusts fragment order to prioritize specific fragments while maintaining the original sequence of others.
    '''
    tmp = string.split('|')
    indices = np.arange(len(tmp)) + 1
    indices = indices.astype(str)
    minIndex = np.argmin(tmp)
    return ('|'.join((*tmp[minIndex:], *tmp[: minIndex])), string, '|'.join((*indices[minIndex:], *indices[:minIndex])))
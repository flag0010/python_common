#modified to work with python3
##############################
import operator
import random
from collections import defaultdict
from functools import reduce
from bisect import bisect_right

###Stats and math functions
def weighted_sampler(pop_dict, k=1):
    """randomly sample a dictionary's keys based on weights stored as values example:
       m = {'a':3, 'b':2, 'c':5}
       samps = weighted_sampler(m, k=1000)
       #samps should be a ~ 300, b ~ 200, and c ~ 500
       >>> samps.count('a')
       304
       >>> samps.count('b')
       211
       >>> samps.count('c')
       485
   of course, being a random sampler your results will vary"""
    vals = list(pop_dict.keys())
    weights = [pop_dict[i] for i in vals]
    return random.choices(vals, weights = weights, k=k)

def choose(n,k):
    '''implements binomial coefficient function
       see: https://en.wikipedia.org/wiki/Binomial_coefficient 
       performance not tested on really large values'''
    return reduce(lambda a,b: a*(n-b)/(b+1),range(k),1)

def sampler(pop, size, replacement=False):
    '''a quick re-implementation of the python random sampler that
       allows for sampling with or without replacement (diff pythno
       builtin perform these two functions)'''
    if replacement:
        #return [random.choice(pop) for i in range(size)] #old way
        return random.choices(pop, k=size) #new way since they added choices, prob faster
    else:
        return random.sample(pop, size)          

def rank(x):
    '''returns the sample rank of the elements in a list.
       if ties, first observed gets lower rank, then second,
       and so on'''
    out={}
    idx=0
    for i in x:
        out[idx] = i
        idx+=1
    p1 =  (j[0] for j in sorted(sort_dict_by_val(out), key=lambda s: s[1]))
    p2 = list(range(len(x)))
    idx=0
    for i in p1:
        p2[i] = idx
        idx+=1
    return p2

def order(x):
    '''returns the sample indeces that would return the list in sorted order
       ie: 
       x = (4,3,406,5)
       sorted(x) == [x[i] for i in order(x)]'''
    out={}
    idx=0
    for i in x:
        out[idx] = i
        idx+=1
    p1 =  [j[0] for j in sorted(sort_dict_by_val(out), key=lambda s: s[1])]
    return p1

###Useful functions for bioinformatics
###NOTE: biopython offers more robust versions, but sometimes you just want something basic

def revcom (s):
    '''returns the reverse complement of a DNA sequence string
       only accepts ACGT, upper or lowercase'''
    trans = str.maketrans('atcgATCG', 'tagcTAGC')
    rv_s = s[::-1] #strange python string reversal, it works!
    rv_comp_s = rv_s.translate(trans)
    return rv_comp_s

def get_fasta(file_name):
    '''read a properly formated fasta and return a dict
       with key=readname and value=sequence
       reads the whole file in'''
    d = [i.strip() for i in open(file_name,'r')]
    out={}
    for i in d:
        if i.startswith('>'):
            curr_seq = i[1:] 
            out[curr_seq] = []
        else:
            out[curr_seq].append(i)
    for i in out:
        out[i] = ''.join(out[i])
    return out

def get_fasta_buffer(file_name):
    '''An efficient fasta reader that is buffered and therefore
       useful for big fasta files.  It returns each fasta one by 
       as a tuple -> (name, sequence). '''
    file_iter = open(file_name)
    current_seq = [] # a dummy, needed to get through the 1st read only
    for line in file_iter:
        if not line.startswith('>'):
            current_seq.append(line.strip())
        else:
            if len(current_seq) != 0:
                yield (current_name, ''.join(current_seq))
            current_name = line[1:].strip()
            current_seq = []
    yield (current_name, ''.join(current_seq))

#Misc
def get_file(filename, splitchar = 'NA', buffered = False):
    if not buffered:
        if splitchar == 'NA':
            return [i.strip().split() for i in open(filename)]
        else: return [i.strip().split(splitchar) for i in open(filename)]
    else:
        if splitchar == 'NA':
            return (i.strip().split() for i in open(filename))
        else: return (i.strip().split(splitchar) for i in open(filename))

def sort_dict_by_val(aDict):
    '''returns a list of tuples sorted by the dict values'''
    return sorted(iter(aDict.items()), key=lambda k_v: (k_v[1],k_v[0]))

def pairwise(li):  
    '''a convienience function that produces all pairwise comparisons from a list'''
    for i in range(len(li)):
        j = i+1
        while j < len(li):
            yield (li[i], li[j])
            j += 1

def count_all(xlist, proportions=False):
    '''Count all the items in a list, return a dict
       with the item as key and counts as value.
       If proportions are set to True, the values 
       are the proportions not counts'''
    out =  defaultdict(int)
    for i in xlist: out[i]+=1
    if proportions:
        out2 = {}
        tot_sz = float(sum(out.values()))
        for i in out: out2[i] = out[i] / tot_sz
        return out2
    else: return out

class rank_list():
    '''maintain ranked list of predetermined size.
       useful when you want sample say the top 5
       things from a long running piece of code
       i.e.
       >>>r = rank_list(maxlen=3)
       >>>r.add(10); r.add(3); r.add(7)
       >>>r.x
        [3, 7, 10]
       >>>r.add(1)
       >>>r.x
        [1, 3, 7]'''
    def __init__(self, maxlen=5):
        self.x = []
        self.maxlen = maxlen
    def add(self, i):
        lx = len(self.x)
        if lx == 0: self.x.append(i)
        else:
            b = bisect_right(self.x, i)
            self.x.insert(b, i)
            self.x = self.x[:self.maxlen]

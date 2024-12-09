import argparse
import random

def load_binding(in_filename):
    infile = open(in_filename,'r')
    line = infile.readline()
    #creating dictionaries that map zf names to their sequences, codons, and scores
    zf_seq = {}
    zf_codon = {}
    zf_score = {}
    linenum = 1
    while line:
        data = line.split('\t')
        if len(data) < 4:
            print('Missing data for ZF on line {}; skipping'.format(linenum))
        else:
            name = data[0]
            seq = data[1]
            codon = data[2]
            score = data[3]
            if name in zf_seq:
                #accounting for circumstance where a single zf can bind to multiple codons
                #just internally counting those as two separate zfs with different codons and scores
                name = name + '\t' + codon

            zf_seq[name] = seq
            zf_codon[name] = codon
            zf_score[name] = float(score)
        linenum += 1
        line = infile.readline()
        
    infile.close()
    return zf_seq, zf_codon, zf_score

def load_binding_no_score(in_filename):
    infile = open(in_filename,'r')
    #creating dictionaries that map zf names to their sequences and codons
    zf_seq = {}
    zf_codon = {}
    linenum = 1
    line = infile.readline()
    while line:
        data = line.split('\t')
        if len(data) < 3:
            print('Missing data for ZF on line {}; skipping'.format(linenum))
        else:
            name = data[0]
            seq = data[1]
            codon = data[2]
            if name in zf_seq:
                #accounting for circumstance where a single zf can bind to multiple codons
                #just internally counting those as two separate zfs with different codons 
                name = name + '\t' + codon

            zf_seq[name] = seq
            zf_codon[name] = codon
        linenum += 1
        line = infile.readline()
        
    infile.close()
    return zf_seq, zf_codon

def load_transitions(transition_file, zf_seq):
    infile = open(transition_file,'r') #formatted with each line as "zf1_name\tzf2_name"
    line = infile.readline()
    zf_transitions = {}
    while line:
        skip = 0
        data = line.split('\t')
        zf1 = data[0]
        zf2 = data[1][:-1]
        if not zf1 in zf_seq.keys():
           #print('ZF {} in transitions but not binding data file; skipping'.format(zf1))
           skip = 1
        if not zf2 in zf_seq.keys():
           #print('ZF {} in transitions but not binding data file; skipping'.format(zf2)) 
           skip = 1
        if not skip:
            if zf1 in zf_transitions:
                zf_transitions[zf1].append(zf2)
            else:
                zf_transitions[zf1] = [zf2]
        line = infile.readline()
    infile.close()
    return zf_transitions
    
def invert_zf_codon(zf_codon):
    codon_zfs = {}
    for zf in zf_codon.keys():
        codon = zf_codon[zf]
        if codon in codon_zfs.keys():
            codon_zfs[codon].append(zf)
        else:
            codon_zfs[codon]=[zf]
    return codon_zfs

def split_seq(sequence):
    #splits a DNA sequence into a list of codons, truncating at the 3' end if not a multiple of 3 bases long
    #and returning the codons from 3' to 5' (the same order as the binding ZFs will be in from N to C terminal)
    #eg split_seq('ATTGACGGG') = ['GGG','GAC','ATT']
    if len(sequence) < 3:
        return []
    if len(sequence) % 3 != 0:
        print('Sequence {} does not have length that is a multiple of three; truncating to {}'.format(sequence, sequence[:-(len(sequence) % 3)]))
        sequence = sequence[:-(len(sequence) % 3)]
    temp_seq =sequence
    output = []
    while temp_seq:
        output.append(temp_seq[-3:])
        temp_seq = temp_seq[:-3]
    return output

class zf_tree:
    #object represents a set of zinc fingers that bind to the same DNA sequence. Each node represents a given zinc finger domain,
    #while the children are themselves trees which are the zfs that can directly follow the parent tree's ZF
    def __init__(self, zf):
        self.zf = zf
        self.children = []
    def add_child(self, child):
        self.children.append(child)
    def return_arrays(self):
        if self.children:
            output = []
            for child in self.children:
                for array in child.return_arrays():
                    output.append([self.zf] + array)
            return output
        else:
            return [[self.zf]]
            
def seq_tree(sequence, missing_codons,codon_zfs,zf_transitions,zf_possible_per_codon,codon_transitions):
    #with seq_tree_recursive, creates a list of trees that represents all zf arrays that can bind to the input sequence. 
    #note that this only looks at the positive strand sequence; should run on reverse complement also (since the ZF can bind in either orientation)
    codon_list = split_seq(sequence)
    for i in range(len(codon_list)):
        if codon_list[i] in missing_codons: #return empty list without starting 
            return None                     #if you know there are codons no zfs bind to in the sequence
        if i < len(codon_list) - 1:
            if not (codon_list[i+1] in codon_transitions[codon_list[i]]): #or if a codon that cannot follow another appears in the sequence
                return None
    trees = []
    for zf_1 in codon_zfs[codon_list[0]]: #for every zf that can bind the first codon
        base_tree = []
        if zf_1 in zf_transitions.keys(): #and which can be followed by any other zfs
            for zf_2 in zf_possible_per_codon[zf_1+codon_list[1]]: #then for every zf that can follow the first zf and bind the second codon
                tree = seq_tree_recursive(codon_list[1:],zf_2, missing_codons,codon_zfs,zf_transitions,zf_possible_per_codon) #create the tree with that zf at the root
                if tree: #if there's a nonempty return, then some zf tree can bind the sequence, and should be added to the tree
                    if not base_tree:
                        base_tree = zf_tree(zf_1)
                    base_tree.add_child(tree)
        if base_tree:
            trees.append(base_tree)
    return trees


def seq_tree_recursive(codon_list,zf_1, missing_codons,codon_zfs,zf_transitions,zf_possible_per_codon,skip_list = {}):
    if len(codon_list) == 1: #if only one codon, return the finger you're starting with since you don't need to add anything else
        return zf_tree(zf_1)
    # if we are on ZF_A, which can be followed by either ZF_B1 or ZF_B2
    # but both ZF_B1 and ZF_B2 can be followed by ZF_C
    # then the tree beyond ZF_C will be the same in both cases and we don't need to generate it twice
    # so for every such ZF_C we add it and the child with it as the root to next_skip
    # and pass that to every subsequent recursive call (in this step) as skip_list
    next_skip = {} 
    if zf_1 in zf_transitions.keys(): #if no zfs can follow the one you're starting with then this is skipped and None is returned
        base_tree = []
        for zf_2 in zf_possible_per_codon[zf_1 + codon_list[1]]: #for every zf that can follow the current zf and bind the next codon
            if not zf_2 in skip_list.keys(): # if we haven't seen it already in this position from the previous step
                tree = seq_tree_recursive(codon_list[1:],zf_2, missing_codons,codon_zfs,zf_transitions,zf_possible_per_codon,next_skip) #create the tree with that zf at the root
                if tree: 
                    for child in tree.children: #note that we've seen all the possible ZFs that can follow this ZF and their subsequent children
                        next_skip[child.zf] = child #add them so as not to redo that work
            else: # if we have seen this ZF in the next position
                tree = skip_list[zf_2] #then we know the child that has to follow it
            if tree: #if there's a nonempty return, then some zf tree can bind the rest of sequence, and should be added to the tree
                if not base_tree:
                    base_tree = zf_tree(zf_1)
                base_tree.add_child(tree)
        if base_tree:
            return base_tree         
            
def num_arrays(tree_list):
    #provides the number of arrays (paths from starting root to final child) in a list of trees from seq_tree
    output = 0
    for tree in tree_list:
        output = output + len(tree.return_arrays())
    return(output)

def grade_array(array,zf_score):
    #provides the score of a zf array (a list of zf names) as the sum of the individual zf scores
    score = 0
    for zf in array:
        score = score + zf_score[zf]
    return score

def ranked_arrays(tree_list,zf_score):
    #returns a list of all the arrays in trees produced by seq_tree ranked by score from highest to lowest
    all_arrays = return_list_arrays(tree_list)
    scored_arrays = [[array,grade_array(array,zf_score)] for array in all_arrays]
    output = sorted(scored_arrays, key = lambda x: x[1],reverse=True)
    return output     
    
def return_list_arrays(tree_list):
    all_arrays = []
    for tree in tree_list:
        all_arrays += tree.return_arrays()
    return all_arrays
    
def write_array(array_list):
    #returns names of ZFs that constitute the array
    output = ''
    for zf in array_list:
        if '\t' in zf: #if zf name was modified due to having multiple affinities
            zf = zf[:-4]    #truncate it, removing trinucleotide affinity and codon
        output = output + zf + '-'
    return output[:-1]
    
def array_sequence(array,zf_seq):
    #given a zf array as a list of zf names, returns the amino acid sequence for that array
    output = zf_seq[array[0]]
    for zf in array[1:]:
        output = output + 'TGERP' + zf_seq[zf]
    return output   

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--target', help='targeted DNA sequence', type=str, required=True)
    parser.add_argument('-o', '--output', help='output file', type=str, required=True)
    parser.add_argument('-b', '--binding', help='ZF binding data, either Z[ifRC], D[eepZF] or custom filename', type=str, required=True)
    parser.add_argument('-m', '--mhc', help='MHC-accepted ZF pairs, either M[ARIA], N[etMHC2Pan] or custom filename', type=str, required=True)
    parser.add_argument('-s', '--score', help='use ZF score or not, default to 1', type=int, required=False, default=1)
    parser.add_argument('-n', help='number of ZF arrays to return, or 0 if all should be returned. if -s is 1 the top scoring arrays'
                                   'will be returned, otherwise a random selection will be, default is 0', type=int, required=False,default=0)
                                  
    
    args = parser.parse_args()
    arguments = vars(args)
    
    if arguments['binding']=='Z':
        binding_file = 'ZifRC_data.txt'
    elif arguments['binding']=='D':
        binding_file = 'DeepZF_data.txt'
    else:
        binding_file = arguments['binding']
    score_use = arguments['score']
    if score_use:
        zf_seq, zf_codon, zf_score = load_binding(binding_file)
    else:
        zf_seq, zf_codon = load_binding_no_score(binding_file)
    
    #inverting zf -> codon dictionary to get all zfs that bind a given codon
    codon_zfs = invert_zf_codon(zf_codon)
                
                
    #creating a list of all possible trinucleotides (64)
    possible_codons = set()
    for n1 in ['A', 'G', 'C','T']:
        for n2 in ['A', 'G', 'C','T']:
            for n3 in ['A', 'G', 'C','T']:
                possible_codons.add(n1+n2+n3)
                
    #creating list of nucleotides no zf in the set binds to
    missing_codons = possible_codons - codon_zfs.keys()
    
    #creating dictionary that maps ZF names to a list of all ZFs that can follow it
    if arguments['mhc']=='M':
        transition_file = 'MARIA_transitions.txt'
    elif arguments['mhc']=='N':
        transition_file = 'NetMHCII_transitions.txt'
    else:
        transition_file = arguments['mhc']
    zf_transitions = load_transitions(transition_file, zf_seq)
    
    #creating dictionary whose keys are concatenations of ZF names and codons (eg "ZNF250 finger 3ATT")
    #whose entries are lists of ZFs that are able to follow the key's ZF and bind the key's codon
    zf_possible_per_codon = {}
    for zf1 in zf_transitions.keys():
        for codon in possible_codons:
            name = zf1 + codon
            zf_list = []
            for zf2 in zf_transitions[zf1]:
                if zf_codon[zf2] == codon:
                    zf_list.append(zf2)
            zf_possible_per_codon[name] = zf_list
    
    codon_transitions = {}
    for zf1 in zf_transitions.keys():
        codon2 = zf_codon[zf1]
        for zf2 in zf_transitions[zf1]:
            codon1 = zf_codon[zf2]
            if codon1 in codon_transitions.keys():
                if not (codon2 in codon_transitions[codon1]):
                    codon_transitions[codon1].append(codon2)
            else:
                codon_transitions[codon1] = [codon2]


    for codon in possible_codons:
        if not codon in codon_transitions.keys():
            codon_transitions[codon] = []
    target = arguments['target']
    tree_list = seq_tree(target, missing_codons,codon_zfs,zf_transitions,zf_possible_per_codon,codon_transitions)
        
    if not tree_list:
        print('No ZFs can target sequence {}'.format(target))
    else:
        outfile = open(arguments['output'],'w')
        if score_use:
            sorted_arrays = ranked_arrays(tree_list,zf_score)
            if arguments['n'] == 0:
                max_arrays = len(sorted_arrays)
            else:
                max_arrays = arguments['n']
            for i in range(max_arrays):
                array = sorted_arrays[i]
                array_list = array[0]
                score = array[1]
                zf_names = write_array(array_list)
                array_seq = array_sequence(array_list,zf_seq)
                outfile.write('>{}|{}\n{}\n'.format(zf_names,score,array_seq))
        else:
            all_arrays = return_list_arrays(tree_list)
            if arguments['n'] == 0:
                max_arrays = len(all_arrays)
            else:
                max_arrays = arguments['n']
            
            random.shuffle(all_arrays) #presumably if you aren't using scores you want a representative sample of the ZF arrays
            for i in range(max_arrays): #so you want arrays which are dissimilar to one another, hence the shuffling
                array_list=all_arrays[i]
                zf_names = write_array(array_list)
                array_seq = array_sequence(array_list,zf_seq)
                outfile.write('>{}\n{}\n'.format(zf_names,array_seq))
        outfile.close()
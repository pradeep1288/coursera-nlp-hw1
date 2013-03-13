#! /usr/bin/python

__author__="Pradeep Nayak <pradeep1288 at gmail dot com"
__date__ ="March 8, 2013"

import sys
from collections import defaultdict
import math
    

def simple_conll_corpus_iterator(corpus_file):
    """
    Get an iterator object over the corpus file. The elements of the
    iterator contain (word, ne_tag) tuples. Blank lines, indicating
    sentence boundaries return (None, None).
    """
    l = corpus_file.readline()
    while l:
        line = l.strip()
        if line: # Nonempty line
            # Extract information from line.
            # Each line has the format
            # word pos_tag phrase_tag ne_tag
            fields = line.split(" ")
            ne_tag = fields[-1]
            #phrase_tag = fields[-2] #Unused
            #pos_tag = fields[-3] #Unused
            word = " ".join(fields[:-1])
            yield word, ne_tag
        else: # Empty line
            yield (None, None)                        
        l = corpus_file.readline()

def sentence_iterator(corpus_iterator):
    """
    Return an iterator object that yields one sentence at a time.
    Sentences are represented as lists of (word, ne_tag) tuples.
    """
    current_sentence = [] #Buffer for the current sentence
    for l in corpus_iterator:        
            if l==(None, None):
                if current_sentence:  #Reached the end of a sentence
                    yield current_sentence
                    current_sentence = [] #Reset buffer
                else: # Got empty input stream
                    sys.stderr.write("WARNING: Got empty input file/stream.\n")
                    raise StopIteration
            else:
                current_sentence.append(l) #Add token to the buffer

    if current_sentence: # If the last line was blank, we're done
        yield current_sentence  #Otherwise when there is no more token
                                # in the stream return the last sentence.

class viterbi(object):
    """viterbi algorithm tagger class"""
    def __init__(self):
        super(viterbi, self).__init__()
        self.n = 3
        self.emission_counts = defaultdict(int)
        self.ngram_counts = defaultdict(int)
        self.word_map = defaultdict(int)
        self.all_states = set()

    
    def build_word_map(self, corpus_file):
        """
        build word map along with their frequencies of occurence
        """
        for line in corpus_file:
            word = line.split(" ")
            if self.word_map.has_key(word[0]):
                self.word_map[word[0]] += 1
            else:
                self.word_map[word[0]] = 1

    def read_counts(self, corpusfile):
        for line in corpusfile:
            parts = line.strip().split(" ")
            count = float(parts[0])
            if parts[1] == "WORDTAG":
                ne_tag = parts[2]
                word = parts[3]
                self.emission_counts[(word, ne_tag)] = count
                self.all_states.add(ne_tag)
            elif parts[1].endswith("GRAM"):
                n = int(parts[1].replace("-GRAM",""))
                ngram = tuple(parts[2:])
                self.ngram_counts[n-1][ngram] = count

    def compute_emission(self, word, ne_tag):
        if self.word_map.has_key(word):
            return self.emission_counts[(word, ne_tag)]/self.ngram_counts[0].get((ne_tag,))
        else:
            return self.emission_counts[("RARE", ne_tag)]/self.ngram_counts[0].get((ne_tag,))


class Tagger(object):
    """
    My generic tagger class
    """
    def __init__(self):
        self.word_map = defaultdict(int)
        self.n = 3
        self.emission_counts = defaultdict(int)
        self.ngram_counts = [defaultdict(int) for i in xrange(self.n)]
        self.all_states = set()

    def build_word_map(self, corpus_file):
        """
        build word map along with their frequencies of occurence
        """
        for line in corpus_file:
            word = line.split(" ")
            if self.word_map.has_key(word[0]):
                self.word_map[word[0]] += 1
            else:
                self.word_map[word[0]] = 1

    def replace_low_freq_words(self, corpus_file, corpus_file_out):
        """
        replace the words with frequencies less than 5 with rare
        """
        for line in corpus_file:
            word = line.split(" ")
            if self.word_map[word[0]] < 5:
                corpus_file_out.write("RARE " + word[1])  
            else:
                corpus_file_out.write(line)
            
    def read_counts(self, corpusfile):
        for line in corpusfile:
            parts = line.strip().split(" ")
            count = float(parts[0])
            if parts[1] == "WORDTAG":
                ne_tag = parts[2]
                word = parts[3]
                self.emission_counts[(word, ne_tag)] = count
                self.all_states.add(ne_tag)
            elif parts[1].endswith("GRAM"):
                n = int(parts[1].replace("-GRAM",""))
                ngram = tuple(parts[2:])
                self.ngram_counts[n-1][ngram] = count

    def compute_emission(self, word, ne_tag):
        if self.word_map.has_key(word):
            return self.emission_counts[(word, ne_tag)]/self.ngram_counts[0].get((ne_tag,))
        else:
            return self.emission_counts[("RARE", ne_tag)]/self.ngram_counts[0].get((ne_tag,))

    def max_emission(self, word):
        if (self.compute_emission(word , "I-GENE") >= self.compute_emission(word , "O")):
            return "I-GENE"
        else:
            return "O"

    def tag_words(self, corpus_file_in, corpus_file_out):
        for line in corpus_file_in:
            if line.strip() == "":
                corpus_file_out.write(line)
            else:
                tag = self.max_emission(line.strip())
                corpus_file_out.write(line.strip() + " " + tag + "\n")

def usage():
    print """
    python count_freqs.py [input_file] > [output_file]
        Read in a gene tagged training input file and produce counts.
    """

if __name__ == "__main__":

    if len(sys.argv)!=2: # Expect exactly one argument: the training data file
        usage()
        sys.exit(2)

    try:
        input = file(sys.argv[1],"r")
        output = file("test.out", "w")
    except IOError:
        sys.stderr.write("ERROR: Cannot read inputfile %s.\n" % arg)
        sys.exit(1)
    tagger = Tagger()
    tagger.build_word_map(input)
    input.seek(0)
    tagger.replace_low_freq_words(input, output)
    input.close()
    output.close()
    gene_count = file("gene.count", "r")
    tagger.read_counts(gene_count)
    gene_count.close()
    
    #Read test data
    gene_input_file = file("gene.test", "r")
    gene_out_file = file("gene_test.p1.out", "w")
    tagger.tag_words(gene_input_file, gene_out_file)
    gene_input_file.close()
    gene_out_file.close()

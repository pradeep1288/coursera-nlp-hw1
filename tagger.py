#! /usr/bin/python

__author__="Pradeep Nayak <pradeep1288 at gmail dot com"
__date__ ="March 8, 2013"

import sys
from collections import defaultdict
import math


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
    gene_input_file = file("gene.dev", "r")
    gene_out_file = file("gene.out", "w")
    tagger.tag_words(gene_input_file, gene_out_file)
    gene_input_file.close()
    gene_out_file.close()
    
    print tagger.max_emission("sent")


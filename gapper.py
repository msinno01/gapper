import os
import random
import numpy as np

class Gapper:
    '''

    A class for producing randomised alignments from a starting alignment.
    Primarily developed for decoy modelling purposes.

    '''
    
    def __init__(self, alignment, no_output=300, out_dir='alignments'):
        self.read_alignment(alignment)
        self.no_output = no_output
        self.out_dir = out_dir

    def read_alignment(self, alignment):
        with open(alignment,'r') as f:
            g = f.read().splitlines()
            f.close()

        self.title1 = g[1]
        self.label1 = g[2]

        ali_seq1 = []
        ali_seq2 = []

        x = 3
        while len(g[x]) == 75 or g[x].endswith('*'):
            ali_seq1.append(g[x])
            x+=1

        self.title2 = g[x+2]
        self.label2 = g[x+3]

        x+=4
        try:
            while len(g[x]) == 75 or g[x].endswith('*'):
                # print(g[x])
                ali_seq2.append(g[x])
                x+=1
        except:
            pass

        self.ali_seq1 = ''.join(ali_seq1)
        self.ali_seq2 = ''.join(ali_seq2)

    def get_gap_no(self,gap):
        if len(gap) < 3:
            gap_nos = [1,2,3,4]
        else:
            gap_nos = [1,2]

        return random.choice(gap_nos)

    def get_new_ali(self,alignment,gap,gap_no):
        new_ali = list(alignment[:])
        for gn in np.arange(gap_no):
            rand_ind = random.randint(0,len(new_ali)-1)
            new_ali.insert(rand_ind, gap)
        return ''.join(new_ali)

    def make_alignment_lines(self, alignment):
        alignment_line = list(alignment)
        x = 0
        c = 0
        for letter in alignment_line:
            x+=1
            c+=1
            if x == 75:
                alignment_line.insert(c, '\n')
                x = 0
        return ''.join(alignment_line)

    def write_new_ali(self, new_ali1, new_ali2, model_no):
        try:
            os.mkdir(self.out_dir)
        except:
            pass
        with open('%s/model%s.ali' % (self.out_dir, model_no), 'w+') as f:
            f.write('\n%s\n%s\n%s\n\n%s\n%s\n%s' % (self.title1, self.label1, new_ali1, self.title2, self.label2, new_ali2))
            f.close()

    def gapify(self):
        gaps = ['-','--','---','----','-----']

        for x in np.arange(self.no_output):
            gap1 = random.choice(gaps)
            gap2 = random.choice(gaps)

            gap_no1 = self.get_gap_no(gap1)
            gap_no2 = self.get_gap_no(gap2)

            new_ali1 = self.get_new_ali(self.ali_seq1, gap1, gap_no1)
            new_ali2 = self.get_new_ali(self.ali_seq2, gap2, gap_no2)

            new_ali_line1 = self.make_alignment_lines(new_ali1)
            new_ali_line2 = self.make_alignment_lines(new_ali2)

            self.write_new_ali(new_ali_line1, new_ali_line2, x)



gapper = Gapper('native.ali', out_dir='alignments_gap')
gapper.gapify()

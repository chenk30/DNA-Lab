from __future__ import annotations
import re

def find_all(strand: str, sequence: str):
    return [(match.start(), match.end()) for match in re.finditer(sequence, strand)]


class StrandHelper:
    inverses = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', ' ': ' '}

    @staticmethod
    def validate(sequence: str, opposite: str = None):
        sequence = sequence.upper()
        for c in sequence:
            if c not in "ACTG ":
                raise Exception(f"Invalid char in sequence: {c}")
        
        if opposite:
            if len(sequence) != len(opposite):
                raise Exception(f"Invalid opposite: len({sequence}) != len({opposite})")
            for c1, c2 in zip(sequence, opposite):
                if c1 != " " and c2 != " " and StrandHelper.inverses[c1] != c2:
                    raise Exception(f"Invalid opposite: {c1} != {c2} in {sequence}, {opposite}")

    @staticmethod
    def get_inverse_str(sequence: str):
        result = []
        for c in sequence:
            result.append(StrandHelper.inverses[c])
        return "".join(result)

# -------------- = strand
# -----          = opposite
class DNA:
    main_strand: str # a strand of DNA starting from 5' to 3'
    opposite:    str # a strand of DNA starting from 3' to 5'
    single:      bool

    def __init__(self, strand: str, opposite: str = None, create_double_helix=False):
        self.main_strand = strand
        self.opposite = " " * len(self.main_strand)

        if create_double_helix:
            self.opposite = StrandHelper.get_inverse_str(self.main_strand)
        elif opposite:
            self.opposite = opposite

        # validate opposite
        StrandHelper.validate(self.main_strand, self.opposite)

        # set if single strand
        self.single = self.opposite == " " * len(self.main_strand)

        # get rid of unneeded spaces
        if self.single:
            self.main_strand = self.main_strand.strip()
        else:
            # trim spaces if there are any
            space_end_index = 0
            space_start_index = len(self.main_strand)
            for i in range(len(self.main_strand)):
                if self.main_strand[i] != " " or self.opposite[i] != " ":
                    space_end_index = i
                    break
            for i in range(len(self.main_strand) - 1, -1, -1):
                if self.main_strand[i] != " " or self.opposite[i] != " ":
                    space_start_index = i + 1
                    break

            self.opposite = self.opposite[space_end_index:space_start_index]
            self.main_strand = self.main_strand[space_end_index:space_start_index]

    def __contains__(self, sequence: str):
        return sequence in self.main_strand or sequence in self.opposite[::-1]

    def __len__(self):
        return len(self.main_strand)

    def __repr__(self):
        if self.single:
            return self.main_strand
        return "\n"+self.main_strand + "\n" + self.opposite

    def will_be_double_helix(self, other: DNA):
        return self.main_strand.get_sequence() == other.main_strand.get_inverse().get_sequence()

    def separate(self):
        # the opposite becomes the main strand, so we flip it's order (3->5 to 5->3)
        separated_strands = [DNA(strand[::-1]) for strand in self.opposite.split()]
        self.main_strand = self.main_strand.strip()
        self.opposite = " " * len(self.main_strand)
        self.single = True
        return separated_strands

    def cleave(self, sequence: str):
        new_dnas = []
        main_strand_start = 0
        opposite_strand_start = 0
        for match_start, match_end in find_all(self.opposite, sequence):
            # we cut "up" at the start and "down" at the end
            # ------------ with seq "AG" will be cut to ----   ---   -----
            #    -AGTAGT                                   -AG   TAG   T

            # split into two new DNAs
            main_strand_end = match_start
            opposite_strand_end = match_end

            new_main_strand = self.main_strand[main_strand_start: main_strand_end]
            new_opposite = self.opposite[opposite_strand_start:opposite_strand_end]
            if len(new_main_strand) > 0:
                if new_opposite[main_strand_end - 1] == " ":
                    new_dnas.append(DNA(new_main_strand))
                    new_dnas.append(DNA(new_opposite[main_strand_end - 1:][::-1]))
                else:
                    if len(new_opposite) > len(new_main_strand):
                        padding = len(new_opposite) - len(new_main_strand)
                        new_main_strand = new_main_strand + " " * padding
                        new_dnas.append(DNA(new_main_strand, new_opposite))
                    else:
                        padding = len(new_main_strand) - len(new_opposite)
                        new_opposite = new_opposite + " " * padding
                        new_dnas.append(DNA(new_main_strand, new_opposite))
            elif len(new_opposite) > 0:
                new_dnas.append(DNA(new_opposite[::-1]))

            main_strand_start = main_strand_end
            opposite_strand_start = opposite_strand_end

        # add dna at the end, if there is one
        if main_strand_start != 0 and main_strand_start != len(self.main_strand):
            new_main_strand = self.main_strand[main_strand_start:len(self.main_strand)]
            new_opposite = len(sequence) * " " + self.opposite[opposite_strand_start:len(self.opposite)]
            new_dnas.append(DNA(new_main_strand, new_opposite))

        # change current dna to one of the cuts, so we won't need to remove it from the tube
        if len(new_dnas) > 0:
            curr_dna = new_dnas.pop()
            self.main_strand = curr_dna.main_strand
            self.opposite = curr_dna.opposite
            self.single == self.opposite.count(" ") == len(self.opposite)

        return new_dnas

    def get_combine_index(self, other):
        other_main = other.main_strand
        other_opposite = other.opposite
        inverse_main = StrandHelper.get_inverse_str(other.main_strand)
        
        # Check main strand against inverse
        for match_start, match_end in find_all(self.main_strand, inverse_main):
            opposite_in_match = self.opposite[match_start:match_end]
            if opposite_in_match.count(" ") == len(opposite_in_match):
                return match_start
                
        # Check sticky ends
        if not other.single:
            # Look for matching regions where opposite strands are empty
            for match_start, match_end in find_all(self.main_strand, other_main):
                if (self.opposite[match_start:match_end].count(" ") == len(self.opposite[match_start:match_end]) and
                    other_opposite.count(" ") == len(other_opposite)):
                    return match_start
                    
        return -1

    def add_opposite(self, start_index: int, strand: str):
        self.opposite = (
            self.opposite[:start_index] + 
            strand + 
            self.opposite[start_index + len(strand):]
        )
        if len(strand) > 0:
            self.single = False

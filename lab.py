from __future__ import annotations
from itertools import combinations
from copy import deepcopy
from random import random
from DNA import *

class Tube:
    def __init__(self, tube=None, filter_fail_percent=0.01, pcr_fail_percent=0.01):
        if not tube:
            self.tube = {}
        else:
            self.tube = tube
        self.filter_fail_percent = filter_fail_percent
        self.pcr_fail_percent = pcr_fail_percent

    def duplicate(self):
        # simluates duplicating the tube to another tube (for example, by amplyfing all strands and pouring half the tube to another one)
        return Tube(deepcopy(self.tube))

    def filter(self, sequence: str):
        # simulate using magnetic beads to filter out strands not containing the given sequence
        new_tube = {}
        for strand in self.tube:
            if sequence in strand and random() >= self.filter_fail_percent:
                new_tube[strand] = self.tube[strand]

        self.tube = new_tube

    def length_sort(self):
        # simulates jel electrophoresis to sort the DNA strands by length
        return [(len(dna), dna) for dna in sorted(self.tube, key=lambda dna: len(dna))]

    def amplify(self, sequence: str, iterations):
        # simluates the PCR process which amplifies all strands containing the given sequence
        # each iterations multiplies the amount of the strand by 2
        for dna in self.tube:
            # 1 percent fail chance
            if sequence in dna and random() >= self.pcr_fail_percent:
                self.tube[dna] *= 2**iterations

    def separate(self):
        # simluate heating the tube to 95 degress so the DNA separates
        to_add = []
        for dna in self.tube:
            amount = self.tube[dna]
            new_strands = dna.separate()
            to_add += [(new_strands, amount)]

        for strands, amount in to_add:
            for strand in strands:
                if strand in self.tube:
                    self.tube[strand] += amount
                else:
                    self.tube[strand] = amount

    def cleave(self, sequence: str):
        # simluate using restriction enzymes to cleave the DNA strands at matching sequences
        to_add = []
        for dna in self.tube:
            amount = self.tube[dna]
            new_strands = dna.cleave(sequence)
            to_add += [(new_strands, amount)]

        for new_strands, amount in to_add:
            for dna in new_strands:
                self.tube[dna] = amount

    def cool(self):
        # simluate cooling the tube to 55 degrees so the DNA strands can combine where possible
        to_remove = []
        to_add = []

        # go over all pairs in tube
        for dna1, dna2 in combinations(self.tube, 2):
            if dna1 not in to_remove and dna2 not in to_remove:
                if dna2.single and len(dna2) <= len(dna1):
                    single_dna = dna2
                    main_dna = dna1
                elif dna1.single and len(dna1) <= len(dna2):
                    single_dna = dna1
                    main_dna = dna2
                else:
                    single_dna = dna1
                    main_dna = dna2

                combine_index = main_dna.get_combine_index(single_dna)
                if combine_index != -1:
                    if self.tube[single_dna] <= self.tube[main_dna]:
                        combined_amount = self.tube[single_dna]
                        self.tube[single_dna] = 0
                    else: # self.tube[dna2] > self.tube[dna1]
                        combined_amount = self.tube[main_dna]
                        self.tube[single_dna] -= self.tube[main_dna]

                    self.tube[main_dna] -= combined_amount

                    combined_dna = deepcopy(main_dna)
                    combined_dna.add_opposite(combine_index, single_dna.main_strand)
                    to_add.append((combined_dna, combined_amount))

                if self.tube[main_dna] == 0:
                    to_remove.append(main_dna)
                if self.tube[single_dna] == 0:
                    to_remove.append(single_dna)

        for dna in to_remove:
            self.tube.pop(dna)

        for dna, amount in to_add:
            self.tube[dna] = amount

    def add(self, other):
        # simluate adding the contents of another tube to this tube
        for dna, amount in other.tube.items():
            if dna in self.tube:
                self.tube[dna] += amount
            else:
                self.tube[dna] = amount

    def __repr__(self):
        return str(self)

    def __str__(self):
        # return a string representation of the tube
        return str(self.tube)
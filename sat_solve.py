from __future__ import annotations
from itertools import product
from lab import *
from random import choices
from pprint import pprint

class SatLiteral:
    def __init__(self, s: str):
        self.s = s
        self.neg = s[0] == "!"
        self.name = s[int(self.neg):]

    def __repr__(self):
        return self.s

    def is_sat(self, val: bool):
        return self.neg ^ val

    def get_inverted(self):
        if self.neg:
            return SatLiteral(self.name)
        return SatLiteral(f"!{self.name}")

class SatClause:
    def __init__(self, clause: str):
        assert(clause[0] == '(' and clause[-1] == ")")
        self.clause_str = clause
        self.literals = [SatLiteral(s) for s in clause[1:-1].split('V')]

    def __repr__(self):
        return self.clause_str

    def is_sat(self, values:list[bool]):
        assert len(values) == len(self.literals)
        return any([literal.is_sat(value) for literal, value in zip(self.literals, values)])

class SAT:
    def __init__(self, cnf: str):
        self.clauses = [SatClause(clause) for clause in cnf.split('^')]
        self.nLiteral = len(set([literal.name for clause in self.clauses for literal in clause.literals]))

    def convert_to_graph(self):
        vertices = set()
        edges = set()

        curr_v_i = 1
        for clause in self.clauses:
            for literal in clause.literals:
                if literal.name not in vertices:
                    left_v=f"a{curr_v_i}"
                    right_v=f"a{curr_v_i+1}"
                    literal_v=literal.name
                    literal_neg_v="!"+literal.name
                    vertices.add(left_v)
                    vertices.add(right_v)
                    vertices.add(literal_v)
                    vertices.add(literal_neg_v)

                    edges.add((left_v, literal_v))
                    edges.add((left_v, literal_neg_v))
                    edges.add((literal_v, right_v))
                    edges.add((literal_neg_v, right_v))

                    curr_v_i += 1

        return Graph(vertices, edges)


class Graph:
    def __init__(self, v: set[str], e: set[tuple[str, str]]):
        self.v = v
        self.e = e

    def _get_dna(self, existing, size):
        found = False
        while not found:
            seq = "".join(choices("AGTC", k=size))
            if seq not in existing and StrandHelper.get_inverse_str(seq) not in existing:
                return seq

    def convert_to_dna(self):
        v_to_dna = {}
        dna_to_v = {}
        dna_len = 2 * len(self.v)
        for v in self.v:
            seq = self._get_dna(v_to_dna, dna_len)
            v_to_dna[v] = seq
            dna_to_v[seq] = v

        e_to_dna = {}
        dna_to_e = {}
        for e in self.e:
            left_v, right_v = e
            dna = StrandHelper.get_inverse_str(v_to_dna[left_v][dna_len // 2:]) + StrandHelper.get_inverse_str(v_to_dna[right_v][:dna_len // 2])
            e_to_dna[e] = dna
            dna_to_e[dna] = e

        return v_to_dna, dna_to_v, e_to_dna, dna_to_e


class SatSolver:
    def __init__(self, cnf:str):
        self.sat = SAT(cnf)
        self.graph = self.sat.convert_to_graph()
        self.nLiteral = self.sat.nLiteral
        self.v_to_dna, self.dna_to_v, self.e_to_dna, self.dna_to_e = self.graph.convert_to_dna()

    def convert_literal_values_to_dna(self, literal_values: list[bool]):
        dna = ""
        dna_opposite = " " * (len(self.v_to_dna["a1"]) // 2)
        for i in range(1, len(literal_values) + 1):
            a_start = f"a{i}"
            a_end = f"a{i+1}"
            chosen_literal = f"x{i}" if literal_values[i-1] else f"!x{i}"
            dna += self.v_to_dna[f"a{i}"] + self.v_to_dna[chosen_literal]
            dna_opposite += self.e_to_dna[(a_start, chosen_literal)] + self.e_to_dna[(chosen_literal, a_end)]
        dna += self.v_to_dna[f"a{self.nLiteral + 1}"]
        dna_opposite += " " * (len(self.v_to_dna["a1"]) // 2)
        return DNA(dna, dna_opposite)

    def convert_literal_values_and_clause_to_strands(self, literal_values: list[bool], clause: SatClause):
        strands = []
        for i, literal in enumerate(clause.literals):
            if literal_values[i]: # assignment is True
                literal_choice = literal.name
            else: # assignment is False
                literal_choice = "!" + literal.name
            strands.append(self.v_to_dna[literal_choice])
        return strands

    def convert_dna_to_literal_values(self, dna: DNA):
        literal_values = []
        for i in range(1, self.nLiteral + 1):
            literal_values.append(self.v_to_dna[f"x{i}"] in dna)
        return literal_values

    def solve(self):
        d = {}
        # add all possible values for the literals as DNA sequences
        for literal_values in product([False, True], repeat=self.nLiteral):
            dna = self.convert_literal_values_to_dna(literal_values)
            d[dna] = 1

        tube = Tube(d)
        print("All possible values:")
        pprint(tube)

        # filter per clause - save only the DNA sequences that satisfy the clause
        for i, clause in enumerate(self.sat.clauses):
            filtered_tubes = []
            for literal_values in product([False, True], repeat=3):
                if clause.is_sat(literal_values):
                    strands = self.convert_literal_values_and_clause_to_strands(literal_values, clause)
                    new_tube = tube.duplicate()
                    for strand in strands:
                        new_tube.filter(strand)
                    filtered_tubes.append(new_tube)
            new_tube = Tube()
            for tube in filtered_tubes:
                new_tube.add(tube)
            tube = new_tube
            print(f"All possible values satisfying clauses until clause {i+1}:")
            pprint(tube)

        possible_solutions = tube.length_sort()
        if len(possible_solutions) > 0:
            print(f"{len(possible_solutions)} satisfying solutions:")
            pprint([dna for _, dna in possible_solutions])
            pprint([self.convert_dna_to_literal_values(dna) for _, dna in possible_solutions])
        else:
            print("No solutions")

    def solve2(self):
        d = {}
        aEnd = f"a{(len(self.v_to_dna) // 2) - 1}"
        dv = {DNA(seq): 128 for v, seq in self.v_to_dna.items() if v != "a1" and v != aEnd}
        seq = self.v_to_dna["a1"]
        dv.update({DNA(seq, StrandHelper.get_inverse_str(seq[:len(seq)//2]) + " " * (len(seq) // 2)): 128})
        seq = self.v_to_dna[aEnd]
        dv.update({DNA(seq, " " * (len(seq) // 2) + StrandHelper.get_inverse_str(seq[len(seq)//2:])): 128})
        de = {DNA(seq): 64  for seq in self.e_to_dna.values()}
        d.update(dv)
        d.update(de)

        tube = Tube(d)
        print(tube)

        tube.cool()
        tube.cool()
        tube.cool()
        tube.cool()

        print(tube)

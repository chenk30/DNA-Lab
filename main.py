from lab import *
from sat_solve import SatSolver

def test():
    # test space trimming
    d = DNA("   ATCG   ")
    print("Main:", d.main_strand)

    # test space trimming with opposite
    d = DNA("   ATCG   ", "   TAGC   ")
    print("Main:", d.main_strand)
    print("Opp:", d.opposite)

    d = DNA("   ATC    ", "     GC   ")
    print("Main:", d.main_strand)
    print("Opp:", d.opposite)

    # test combine
    d = {DNA("TACGA"): 1, DNA("ATGC"): 2, DNA("ATACG"): 1}
    t = Tube(d)
    t.cool()
    print("Cool result", t)

    # test cleave
    t.cleave("TG")
    print("Cleave result", t)

    # test separation
    t.separate()
    print("Separate result", t)

def main():
    sat ="(x1Vx2Vx3)^(x1Vx2V!x3)^(x1V!x2Vx3)^(x1V!x2V!x3)^(!x1Vx2Vx3)^(!x1Vx2V!x3)^(!x1V!x2Vx3)"
    #sat ="(x1Vx2Vx3)^(!x1V!x2Vx3)^(x1V!x2V!x3)"
    #sat = "(x1Vx2Vx3)"
    solver = SatSolver(sat)
    solver.solve()

if __name__ == '__main__':
    main()

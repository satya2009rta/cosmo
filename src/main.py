from distutils.cmd import Command
import sys
import os
import importlib
import subprocess
import re
import time

MAIN_DIR = os.path.dirname(os.path.abspath(__file__))

def ltlsynt(APins, LTL):
    print("----------------")
    print(f'Calling ltlsynt for formula \n {LTL} \n with inputs \n {APins}\n')
    start = time.time()

    out = subprocess.run(["ltlsynt",
                            f"--ins={','.join(APins)}",
                            f"-f {LTL}",
                            "--decompose=no",
                            "--print-game-hoa=t"],
                            capture_output=True,
                            text=True)
    if out.returncode != 0:
        print(f'Error: ltlsynt errored\n{out.stderr}', file=sys.stderr)
        sys.exit(out.returncode)
    print(f"Done in {time.time() - start} s")
    game = out.stdout
    return game

def solveltl(LTL):
    start = time.time()

    out = subprocess.run(["ltlsynt",
                            "--ins=jl50",
                            f"-f {LTL}",
                            "--decompose=no"],
                            capture_output=True,
                            text=True)
    # print(f'-----\nSynthesizing strategy for formula \n {LTL}\n')
    if out.returncode != 0 and out.returncode != 1:
        print(f'Error: negotiate errored\n{out.stderr}', file=sys.stderr)
        sys.exit(out.returncode)
    elif out.returncode != 0:
        print ("UNREALIZABLE with cooperation!")
    print(f"Cooperative synthesis done in {time.time() - start} s")
    strategy = out.stdout
    return strategy


def negotiate(game0, game1):
    print("----------------")
    print("Negotiating between players...\n")
    start = time.time()

    command = MAIN_DIR + "/../build/negotiate.o"
    args = [f"{game0}", f"{game1}"]
    out = subprocess.run([command],
                          input=''.join(args),
                          capture_output=True,
                          text=True)
    if out.returncode != 0 and out.returncode != 1:
        print(f'Error: negotiate errored\n{out.stderr}', file=sys.stderr)
        sys.exit(out.returncode)
    elif out.returncode == 0:
        print("UNREALIZABLE!")
    elif out.returncode == 1:
        print("#SUCCESS : initial vertex is winning")
    print(f"\nDone in {time.time() - start} s")
    print("----------------\n")
    contract = out.stdout
    return contract


def main(LTL0, LTL1, APins):
    
    start = time.time()

    game0 = ltlsynt(APins, LTL0)
    game1 = ltlsynt(APins, LTL1)

    contract = negotiate(game0,game1)

    print(contract)

    print(f"COMPLETED in {time.time() - start} s\n")

    strategy = solveltl(f"( {LTL0} ) & ( {LTL1} )")

    print(game0)
    print(game1)


def dummyLTLgen(k,i):
    if k == 1:
        return f"( X m{i} )"
    return f"( ( X ! m{i} ) & X {dummyLTLgen(k-1,i)} )"



if __name__ == '__main__':
    rc = 1
    try:
        if len(sys.argv) > 1:
            k1 = int (sys.argv[1])
            k2 = int (sys.argv[2])
            APins = ["go1","m1"]
            LTL1 = f" ( G ! ( go1 & go2 ) ) & G F go1 & m1 & G ( m1 => {dummyLTLgen(k1,1)})"
            LTL2 = f" ( G ! ( go1 & go2 ) ) & G F go2 & m2 & G ( m2 => {dummyLTLgen(k2,2)})"
        else:
            LTL1 = input("\nSpecification of Player 1 (who plays first) -- in LTL:\n")
            LTL2 = input("\nSpecification of Player 2 (who plays second) -- in LTL:\n")
            APins = input("\n1st varibales of the plays: Output Atomic Propositions of Player 1 (who plays first) -- separated by space:\n").split()
        main(LTL1, LTL2, APins)
        rc = 0
    except Exception as e:
        print('Error: %s' % e, file=sys.stderr)
    sys.exit(rc)
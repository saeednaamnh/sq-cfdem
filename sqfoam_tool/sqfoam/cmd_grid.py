import argparse
from pathlib import Path
from config import SHAPES
from cmd_new import cmd_new

ST_VALUES=[1000,100,30,10,3,1]

def cmd_grid(a):
    if not a.full:
        print("use --full to generate all 24 cases"); return
    made=0
    for shape in SHAPES:
        for St in ST_VALUES:
            ns=argparse.Namespace(shape=shape,St=float(St),name=None,
                                  steps=5000000,out=a.out)
            try: cmd_new(ns); made+=1
            except SystemExit: pass
    print(f"\ngrid: attempted {len(SHAPES)*len(ST_VALUES)} cases into {a.out}/")

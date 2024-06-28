from pathlib import Path
from argparse import ArgumentParser

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from arip3.ARIP_main import load_atom_models, parse_atom_model
from arip3.typing import DTYPE


def vis_xyz(args):
    points = np.loadtxt(args.f, dtype=DTYPE)
    print('points.shape:', points.shape)
    p0 = points[0]
    print('|p0|:', np.linalg.norm(p0))

    ax = plt.axes(projection='3d')
    ax.scatter3D(points[:, 0], points[:, 1], points[:, 2], s=1)
    plt.savefig(f'{args.f.stem}.png', dpi=300)


def vis_pdb(args):
    for idx, a_model in load_atom_models('', args.f):
        plt.figure()
        df_tmp = parse_atom_model(a_model, 0)
    
        df = pd.concat(df_tmp, axis=0)
        ax = plt.axes(projection='3d')
        ax.scatter3D(df['x'], df['y'], df['z'], s=1)
        if idx == -1: plt.savefig(f'{args.f.stem}.png', dpi=300)
        else:         plt.savefig(f'{args.f.stem}_MODEL_{idx}.png', dpi=300)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('f', type=Path, help='*.pdb or *.xyz file')
    args = parser.parse_args()

    fp: Path = args.f
    try: 
        if    fp.suffix == '.xyz': vis_xyz(args)
        else: vis_pdb(args)
    except:
        print(f'unsupported file: {fp}')

from pathlib import Path
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
import os

from arip3 import arip_main

def process_file(in_fp, o, ref, dx, t, d, p, a, w, r, s, z, disable_print):
    try:
        if not disable_print:
            print(f'>> process: {in_fp}, ref: {ref}')
        arip_main(in_fp, o, ref, dx, t, d, p, a, w, r, s, z, disable_print)
    except SystemExit:
        print(f'Error processing file: {in_fp}')
    except:
        print(f'The file {in_fp} cannot be processed')

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-ref', default=Path(__file__).parent / 'data/s5000.xyz', type=Path, help='*.xyz coord file, sphere grid points')
    parser.add_argument('input', type=Path, help='input *.pdb file or folder path')
    parser.add_argument('-o', default=Path(__file__).parent / 'out', type=Path, help='output path')
    parser.add_argument('-dx', default=0.2, type=float, help='cubic grid distance')
    parser.add_argument('-e', action='store_true', help='enhanced precision, use 15092 dots and 0.1 interval as surface and volume')
    parser.add_argument('-a', action='store_true', help='calculate solvent accessible surface area')
    parser.add_argument('-c', nargs='*', default=None, type=float, help='surface and volume lower cutoff, two values')
    parser.add_argument('-d', nargs='?', default=None, const=0, type=float, help='closest distance between two atoms used to determine contact')
    parser.add_argument('-t', default=os.cpu_count(), type=int, help='number of threads')
    parser.add_argument('-p', action='store_true', help='analyze only the interactions mediated by water molecules of hydrophilic atoms')
    parser.add_argument('-w', action='store_true', help='use atomic overlapping weighted algorithm for volume')
    parser.add_argument('-r', action='store_true', help='generate .csv file for every residue')
    parser.add_argument('-s', action='store_true', help='only calculate surface area')
    parser.add_argument('-z', action='store_true', help='use .gz compressed format save result')
    args = parser.parse_args()

    # Enhanced precision mode
    if args.e: 
        args.ref = Path(__file__).parent / 'data/s15092.xyz'
        args.dx = 0.1

    input: Path = args.input
    assert input.is_file() or input.is_dir(), f'{input} is not a valid file or directory'
    
    # Optional threshold, no default
    threshold = args.c

    if args.c is not None:
        # Use the default threshold
        if len(args.c) == 0:
            threshold = [0.5, 0.2]
        # The user provided two values, use the user's values
        elif len(args.c) == 2:            
            threshold = list(map(float, args.c))
        # The user did not provide the correct values
        else:            
            print('Lower cutoff must be TWO values like 1.0 0.5, or leave it blank to use the default values 0.5 0.2')
            exit(1)
    
    # The closest distance for determining whether two atoms are in contact
    dist = args.d
    
    # The number of threads should not exceed the number of CPUs. If the user enters a larger number, the number of threads is determined by the number of CPUs.
    num_threads = min(args.t, os.cpu_count())
    
    # Only analyze the interactions mediated by water molecules of hydrophilic atoms
    polar = True if args.p else False
    if polar:
        args.ref = Path(__file__).parent / 'data/s15092.xyz'
        args.dx = 0.05

    # Calculate SASA
    accessible = True if args.a else False
    
    # Use atomic overlapping weighted algorithm
    weighted = True if args.w else False
    
    # Save .csv for every residue
    each = True if args.r else False
    
    # Calculate surface area only
    only = True if args.s else False
    
    # Save the result in .gz compressed format
    compress = True if args.z else False
    
    # Don't print details for dir
    in_fps = [input] if input.is_file() else list(input.iterdir())
    disable_print = False if input.is_file() else True
    
    # Create output directory
    output_dir = args.o.resolve()
    print(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        executor.map(process_file, in_fps, [args.o]*len(in_fps), [args.ref]*len(in_fps),
                     [args.dx]*len(in_fps), [threshold]*len(in_fps), [dist]*len(in_fps),
                     [polar]*len(in_fps), [accessible]*len(in_fps), [weighted]*len(in_fps),
                     [each]*len(in_fps), [only]*len(in_fps), [compress]*len(in_fps),
                     [disable_print]*len(in_fps))
        
        
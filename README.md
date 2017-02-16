# vneighbor.cpp - `g++ -o vneighbor vneighbor.cpp` to compile the code

## Usage

`vneighbor FileName1 FileName2 min-distance max-distance`
-FileName1:   bulk supercell
-FileName2: defect supercell

e.g.,

- `vneighbor POSCAR1 POSCAR2 1.2 3.0` (for 1st NNs of fcc, a=3.54)
- `vneighbor POSCAR1 POSCAR2 3.0 3.9` (for 2nd NNs of fcc, a=3.54)

The input file must be in the following POSCAR format (no other lines like select relaxation):

    Headline of POSCAR
    3.5
    1 0 0
    0 1 0
    0 0 1
    Ni Fe
    2 2
    D
    .0 .0 .0
    .0 .5 .5
    .5 .0 .5
    .5 .5 .0

The code then outputs the neighbors of the vacancy within the given shell. The output file is `vneighbor_min-max`.

## Output

The output file can be divided into four parts:

1. Atomic index and species of the vacancy in bulk.
2. The number of bondings for the vacancy with any element.
3. The index of neighbor atoms within that shell for the vacancy. Index is the sequence in POSCAR.
4. The distance between the neighbor atoms and the atom used to be in the bulk;
   The distance between the neighbor atoms and the vacancy.

## Tricks in practice

If you don't know the distance of 1st shell, just make a guess, like 1.0-4.0, the code will print all distances. Then you can see the distance of 1st, or even 2nd, neighbors.

If you want to study only a certain type or types of atom in a crystal, delete the other types of atom from your POSCAR. Then use this code to calculate the neighborhood.

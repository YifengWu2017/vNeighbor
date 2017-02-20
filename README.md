# vneighbor.cpp
- `g++ -o vneighbor vneighbor.cpp` to compile the code

## Usage

`vneighbor FileName1 FileName2 min-distance max-distance`
- FileName1:   bulk supercell
- FileName2: defect supercell

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

NAME
    FES-naive

DESC
    a naive implementation of fast exhaustive search over GF(2)

EXAMPLE
    > cmake
    > make
    > ./FastExhaustiveSearch_naive
            field: GF(2)
    		number of variables: 20
    		number of equations: 40
    		seed: 2
    		reading equations...
    		brute forcing...
    		found valid solution: 30468
    		solution valid
    		solution: [0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0]

 REF
    [1] Cahrles Bouill, et al. Fast Exhaustive Search for Polynomial Systems in F2
            https://eprint.iacr.org/2010/313.pdf
    [2] Cahrles Bouill, et al. Fast Exhaustive Search for Quadratic Systems in F2 on FPGAs
            https://www.win.tue.nl/~tchou/papers/forcemp-fpga.pdf
    [3] GPU-version of libFES
            http://polycephaly.org/projects/forcemq/data/MQForceGPU-20180724.tgz
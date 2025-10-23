# Incremental_Driver: Constitutive Model Calibrator and Tester in Modern Fortran

## Compiling

A `fpm.toml` file is provided for compiling incremental-driver with the [Fortran Package Manager (fpm)](https://github.com/fortran-lang/fpm). To install fpm, I recommend installing it with conda. You can install conda [here](https://www.anaconda.com/docs/getting-started/miniconda/install) 

To build the program with fpm use
```
fpm build
```

This will build the program in debug mode. Which is likely what you want so that the debugger can step into your umat.

To run the unit tests:

```
fpm test
```

To generate the documentation using [ford](https://github.com/Fortran-FOSS-Programmers/ford), run: ```ford ford.md```


The latest API documentation can be found [here (Not Made yet)](). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

## License

The finterp source code and related files and documentation are distributed under a permissive free software [license](https://github.com/CriticalSoilModels/Incremental_Driver/LICENSE) (BSD-style).

# Note
Currently a large part (almost all) of the incremental driver software was taken directly from: [Niemunis SoilModels](https://soilmodels.com/idriver/)

Some slight modifications have been made so that it can compile using as a .f90 file and to
move some of the functions to another file to make it easier to see the main program.

A big thank you to Andrzej Niemunis for all the work done to write this code.

Please cite Andrzej Niemunis if this code is used in your research.
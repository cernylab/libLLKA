Classify and write CIF example
===

This is an example of a program that uses libLLKA to load a structure from a mmCIF file, performs NtC assigment and prints the result to standard output. The results is a new mmCIF file with additional categories that contain the results of the NtC assignment.

Build and use
---
This program is built automatically if libLLKA is built with `BUILD_EXAMPLES` option set to `ON`. To run the program, execute

```
./classify_and_write_cif <input_file.cif>
```

Note that the [assignment parametrization files](../../README.md#NtC-assignment-parametrization) must be present in the working directory (the directory from which the program is launched) for the assignemnt to work. Otherwise, the program will fail to load the parametrization files.

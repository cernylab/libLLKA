Examples
===

libLLKA currently provides the following examples to demonstrate the library's capabilities and API:

- [mirror_cif](mirror_cif)

  A simple code that loads and parses a mmCIF file and prints its contents back to standard output. This example demonstrates how to load and read mmCIF files using libLLKA.

- [simple_NtC_assignment](simple_NtC_assignment)

  An example program that reads a structure from a mmCIF file, breaks the structure down to __dinucleotide steps__ and performs NtC assignment for each __dinucleotide step__. The results are written as a simple list to the standard output.

- [classify_and_write_cif](classify_and_write_cif)

  A program that performs the same operation an `simple_NtC_assignment` but outputs the result as mmCIF-formatted text that contains additional categories with the results of the NtC assignment. This example program can be readily used to annotate a structure with NtC assignment data and store the output in a commonly used data format.

- [gui_assigner](gui_assigner)

  A program that performs NtC assigment but also comes with a simple graphical user interface to display and explore the results.


Building examples
---
All examples are built automatically if the `BUILD_EXAMPLES` option is set to `ON` during configuration. This is the default. You can also pass `-DBUILD_EXAMPLES=ON` parameter on the command line to make sure that the option is enabled.

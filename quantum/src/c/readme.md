=== Compilation ===
Copy the machine-specific Makefile local_MACHINE.mk to local.mk and edit if necessary. "make all" will build the code, "make doc" will build the documentation in the subdirectory doc.

=== Running ===
To run, type

./driver.x

All parameters are read from the file parameters.in, so edit this accordingly. The order of parameters matters (and the labels are actually ignored). You might want to copy the file parameters_template.in to parameters.in to make a start.
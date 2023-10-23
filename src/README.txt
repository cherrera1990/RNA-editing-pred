The files in this folder are the C++ source code for all the programs (except linearfold) used to create the initial RNA-editing datasets
for both the Random Foret and Neural Networks approaches. For details on the overall pipeline and order of execution, conssult the supplementary
methods in Zawisza et al., 2023. For details about what each program does, consult the comment at the start of the source code of each program.

Before execution, the programs need to be compiled to create an executable file. You can use any C++ compiler, but an example in linux would be
using the g++ compiler:
g++ program.cpp -o program.x

Then, in order to run the program, if the executable file is in thw working directory, the command should be:
./program.x [arguments]
If you are running the program from a different directory, you can run it with the path to the file instead.

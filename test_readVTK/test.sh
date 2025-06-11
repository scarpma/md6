make
cc -c test_readVTK/test_readVTK.c -o test_readVTK/test_readVTK.o
cc exec/write_to_vtk.o test_readVTK/test_readVTK.o -o test_readVTK/test_readVTK
./test_readVTK/test_readVTK

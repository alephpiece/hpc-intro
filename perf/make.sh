CC=g++

$CC mm.cxx -o mm.exe
$CC -O3 -fopt-info-vec mm.cxx -o vec_mm.exe
$CC -fopenmp par_mm.cxx -o par_mm.exe
$CC -fopenmp -O3 -fopt-info-vec par_mm.cxx -o par_vec_mm.exe

CC=hipcc
$CC simt_mm.cxx -o simt_mm

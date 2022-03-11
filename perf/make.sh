CC=g++

echo "Compiling mm.cxx"
$CC -march=native mm.cxx -o mm.exe

echo "Compiling mm.cxx with vectorization"
$CC -march=native -O3 -fopt-info-vec mm.cxx -o vec_mm.exe

echo "Compiling par_mm.cxx with multithreading"
$CC -march=native -fopenmp par_mm.cxx -o par_mm.exe

echo "Compiling par_mm.cxx with multithreading & vectorization"
$CC -march=native -fopenmp -O3 -fopt-info-vec par_mm.cxx -o par_vec_mm.exe

echo "Compiling simt_mm.cxx with SIMT (HIP)"
hipcc simt_mm.cxx -o simt_mm

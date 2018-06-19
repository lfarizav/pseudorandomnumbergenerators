# pseudorandomnumbergenerators
You can use Box-Muller and Ziggurat algorithms optimized using SIMD instructions.

# Prerequisites
1. Intall the Math Kernel Library (MKL) from Intel: source l_mkl_2018.3.222/./install.sh
2. Run the script compilervars: source /opt/intel/bin/./compilervars.sh intel64
# Compilation Box-Muller
gcc -o boxmuller boxmuller.c gaussdouble.c display.c -lm -mavx2 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
# Compilation Ziggurat
gcc -o ziggurat ziggurat.c gaussdouble.c display.c -lm -mavx2 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread

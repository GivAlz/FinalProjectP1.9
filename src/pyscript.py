import ctypes

# gcc -shared -Wl,-install_name,test1.so -o test1.so -fPIC test1.c

cprog = CDLL('./ljmd.so')
cprog.test_few_atoms()

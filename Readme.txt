sbox10.cpp : 
// This code perfoms the steepest descent-like iterative search algorithm in [1] within the space of
// rotation symmetric S-boxes (RSSBs), 2-RSSBs, 5-RSSBs, and the concatenation of RSSBs in dimension 10.
// The files "mersenne.cpp" and "randomc.h" are used to generate random S-box.
// The file "H10.txt" contains the Walsh-Hadamard matrix.
// The file "ANF10.txt" is used to compute algebraic degree.
// Result are written to the file "S10_P.txt".

lat_ddt_10.cpp :
// This code finds difference distribution table (DDT) and linear approximation table (LAT) 
// of some rotation symmetric S-boxes and their concatenations (in dimension 10), which are 
// obtained by the steepest descent like iterative search algorithm [1].
// Result are written to the files "DDT.txt" and "LAT.txt".


DDT_Si.txt and LAT_Si.txt (i=1,2,3,4): The whole DTT and LAT tables of S1, S2, S3, and S4 obtained in [1]. 
DDT_S5.txt and LAT_S5 : DDT and LAT tables of the AES S-box.

//S1 (rotation symmetric) with nonlinearity 456, absolute indicator 192, algebraic degree 9, and differential uniformity 10 
//S2 (rotation symmetric) with nonlinearity 454, absolute indicator 184, algebraic degree 9, and differential uniformity 8 
//S3 (concatenation) with nonlinearity 456, absolute indicator 192, algebraic degree 9, and differential uniformity 12 
//S4 (concatenation) with nonlinearity 454, absolute indicator 184, algebraic degree 9, and differential uniformity 8 
//S5 (AES S-box)

[1] Selçuk Kavut. Cryptographically Optimized Large S-boxes in Some Subspaces (in Turkish). Submitted to EMO Bilimsel Dergi.

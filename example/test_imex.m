%%
clear
gen_mex('imexTest.cpp');

%%
N = 10;
A = randi(100, N, 'int32');
B = eye(N);

%%
C = -1j*imexTest(A, 1j*B);
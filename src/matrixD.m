function K = matrixD(Z,ZD,ZDD,N,d)
K = zeros(N);
%todo: vectorize this 
f = 0.25/pi;
for aa = 1:N
     bb = (1:N ~= aa);
     K(aa,aa) =     0.5 + f*imag(ZDD(aa)/ZD(aa));
     K(aa,bb) =         + f*imag(ZD(aa)*cot( 0.5*(Z(aa)-Z(bb))) );   
     K(aa,1:N)  = K(aa,1:N) - f*imag(ZD(aa)*cot( 0.5*(Z(aa) - conj(Z) + 2i*d) ));
end
end
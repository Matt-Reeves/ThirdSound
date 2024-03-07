function w = velocity(a,ad,Z,ZD,ZDD,N,d)

% TODO: Check if this can be optimized. Roberts mentions calculating cot(Z) 
%and using difference formulae for cotangents. 

w = zeros(1,N);
nall = 1:N;
i_on_4pi = 1i*0.25/pi;
for jj = nall
    offd = (1:N ~= jj);
    
    w(jj) = -i_on_4pi*(  sum(a(offd).*cot( 0.5*(Z(jj)  -      Z(offd)        ))) ...
                        -sum(a(nall).*cot( 0.5*(Z(jj)  - conj(Z(nall)) + 2i*d))) ); 
end

%N.B.: Typo in Mercer and Roberts [Physics of Fluids 11, 1051 (1999)] -- 
%factor of 1i/4/pi is missing on second two terms in the paper.
w = w + (0.5*a - i_on_4pi*( ZDD.*a./ZD - 2*ad ))./ZD;
end
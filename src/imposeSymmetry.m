function [x,y,phi] = imposeSymmetry(x,y,phi,N)

check = size(x);
if check(1) > check(2) 
    error('Error: x must be row vector')
end

check = size(y);
if check(1) > check(2) 
    error('Error: y must be row vector')
end

check = size(phi);
if check(1) > check(2) 
    error('Error: phi must be row vector')
end
 
x(1) = 0; 
x(N/2+1) = pi; 
x(N/2+2:N) = 2*pi - x(N/2:-1:2);
y(N/2+2:N) = fliplr(y(2:N/2));
phi(N/2+2:N) = fliplr(phi(2:N/2)); 

end
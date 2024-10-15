function L = DiffM(m)
L = zeros(m,m);
for i = 1:m-1
    L(i,i) = 1;
    L(i,i+1) = -1;
end
L(m,1) = -1;
L(m,m) = 1;


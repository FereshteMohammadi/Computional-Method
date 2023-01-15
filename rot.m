function A = rot(fi)
sfi = sin(fi);
cfi = cos(fi);
A = [cfi, -sfi
    sfi, cfi];
end
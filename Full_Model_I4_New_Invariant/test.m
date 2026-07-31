% test function

function out = test(a,b,c)

if (~exist('s1','var'))
    s1 = 0;
end

out = a+b;

if out>0
    out = c;
end

end
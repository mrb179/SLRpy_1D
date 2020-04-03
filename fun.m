function [ w ] = fun( z )
%  fun: function definition
%
% INPUTS
%
%  z     : function argument
%
% OUTPUTS
%
%  w     : function value
%

w = py.pyfunc.f(real(z), imag(z));
disp(w);
disp(class(w));
%disp(I(1));
%w = complex(I(1), I(2));
%disp(w);

end






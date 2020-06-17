

% rectangular domain definition z=x+jy x\in[xb,xe] , y\in[yb,ye]
xb = 3.0;     % real part range begin 
xe = 4.0 ;     % real part range end
yb = -0.1;     % imag part range begin
ye = 0;      % imag part range end 
r = 0.01;     % initial mesh step

NewNodesCoord = rect_dom(xb,xe,yb,ye,r); % initial mesh generation

Tol = 1e-4; % accuracy (candidate region size)

visual=0; % mesh visualization:  0 - turned off,   1  - only last iteration,   2  - all iterations

ItMax=100; % max number of iterations

NodesMax=500000; % max number of nodes

SkinnyTriangle=3; % sinny triangle definition


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Lattice Parameters	%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


RNP  = 25.000000;	% nanoparticle radius [nm]

Nind = 1.000000;	% background refractive index




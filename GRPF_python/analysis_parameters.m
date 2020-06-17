

% rectangular domain definition z=x+jy x\in[xb,xe] , y\in[yb,ye]
xb = 2.4635*0.9;     % real part range begin 
xe = 2.4635*1.1 ;     % real part range end
yb = -0.1;     % imag part range begin
ye = 0;      % imag part range end 
r = 0.01;     % initial mesh step

NewNodesCoord = rect_dom(xb,xe,yb,ye,r); % initial mesh generation

Tol = 1e-5; % accuracy (candidate region size)

visual=0; % mesh visualization:  0 - turned off,   1  - only last iteration,   2  - all iterations

ItMax=100; % max number of iterations

NodesMax=500000; % max number of nodes

SkinnyTriangle=3; % sinny triangle definition


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Lattice Parameters	%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_parallel = 9.519978e-03; % in-plane wave vector [nm^-1]

RNP  = 20.000000;	% nanoparticle radius [nm]

Nind = 1.500000;	% background refractive index

d    = 330.000000;	% unit cell length [nm]

t    = 165.000000;	% spacing between NPs within the unit cell [nm]







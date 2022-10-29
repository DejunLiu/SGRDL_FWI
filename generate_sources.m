function b = generate_sources(dom,sources)
% GENERATE_SOURCES   Create sources corresponding to physical locations.
%   B = GENERATE_SOURCES(DOM,SOURCES) returns a matrix B whose jth column
%   is a point source specified by the jth column of SOURCES on the grid
%   associated with the DOMAIN object. SOURCES should ideally be generated
%   by SOURCES_AND_RECEIVERS.
%
% References:
%A. Cosse,  S. D. Shank, and L. Demanet.  A short note on rank-2 relaxation for waveform inversion, 
%Seg Technical Program Expanded Abstracts, 2015, 1344-1350.


ns = size(sources,2);		% number of sources
b = zeros(dom.N,ns);		% initialize right-hand sides
% Create discrete deltas at source location
for j = 1:ns
	b(:,j) = dom.pt_src(sources(1,j),sources(2,j));
end
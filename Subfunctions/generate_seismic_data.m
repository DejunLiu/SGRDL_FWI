function d = generate_seismic_data(dom,sources,receivers,all_freqs,m_true,sigma)
% GENERATE_SEISMIC_DATA   Creating seismic data for given parameters.
%   D = GENERATE_SEISMIC_DATA(DOM,SOURCES,RECEIVERS,ALL_FREQS,M_TRUE,SIGMA)
%   returns a tensor D such that D(:,i,j) is noisy data sampled at
%   receivers at locations specified by RECEIVERS. DOM is an instance of
%   the domian object.
%
%   In the above, index i denotes the ith frequency in the vector ALL_FREQS
%   and index j denotes the jth source at the locations specified by
%   SOURCES.  SOURCES and RECEIVERS can easily be generated by.
%   SOURCES_AND_RECIEVERS. M_TRUE is the true background squared slowness
%   and SIGMA is the noise level.
%
%   See also DOMAIN, SOURCES_AND_RECEIVERS.

nf = length(all_freqs);			% number of frequencies
ns = size(sources,2);			% number of sources
nr = size(receivers,2);			% number of receivers
d = zeros(nr,nf,ns);			% seismic data
r_ind = dom.loc2ind(receivers);	% flat indices of receivers

% Generate true data
fprintf('Generating data... frequency: ')
rand('seed',0);% set random seed
parfor i = 1:nf
	fprintf('%d, ',i)
	% Invert Helmholtz operator corresponding to true model
	A_true = invertA(helmholtz_2d(m_true,all_freqs(i),dom));
	for j = 1:ns
		% Apply to jth source
		u_true = A_true.apply(dom.pt_src(sources(1,j),sources(2,j))); %#ok<*PFBNS>
		% Gather seismic data
		d(:,i,j) = u_true(r_ind)+sigma.*(rand(nr,1)+1i*rand(nr,1));
	end
end

function set_path();
%function [F,p]=sippi_set_path();
[p]=fileparts(which('set_path.m'));
if (isempty(p)|strcmp(p,'.'))
    p=pwd;
end

% adding toolboxes shipped with PQN
i=0;
i=i+1;F{i}=p;
i=i+1;F{i}=[p,filesep,'omp2'];
i=i+1;F{i}=[p,filesep,'tool'];

for i=1:length(F);
	
    if exist(F{i},'dir')
        try
            addpath(F{i});
            stat='OK';
        catch
            stat='FAILED';
        end
        disp(sprintf('trying to add path %s [%s]',F{i},stat));
    else
        stat='NONEXISTENT';
        disp(sprintf('DIRECTORY DOES NOT EXIST : %s ',F{i}))
    end
    
end
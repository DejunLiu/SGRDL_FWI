
function set_path();
%function [F,p]=sippi_set_path();
[p]=fileparts(which('set_path.m'));
if (isempty(p)|strcmp(p,'.'))
    p=pwd;
end

% adding toolboxes shipped with PQN
i=0;
i=i+1;F{i}=p;
%i=i+1;F{i}=[p,filesep,'TSDL',filesep,'omp2',filesep,'private'];
i=i+1;F{i}=[p,filesep,'omp2'];
i=i+1;F{i}=[p,filesep,'ksvd'];
i=i+1;F{i}=[p,filesep,'Test_Images'];
i=i+1;F{i}=[p,filesep,'Subfunctions'];
i=i+1;F{i}=[p,filesep,'tool'];
i=i+1;F{i}=[p,filesep,'marm_data'];
i=i+1;F{i}=[p,filesep,'BG_data'];
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
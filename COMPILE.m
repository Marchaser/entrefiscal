function COMPILE

% mex -g update_dist.cpp -DUSE_OPENMP COMPFLAGS="$COMPFLAGS /Qopenmp /Qansi_alias"
% if ~isequal(exist('update_dist','file'),3)
mex update_dist.cpp  -DUSE_OPENMP "C:\Program Files (x86)\Intel\Composer XE\mkl\lib\intel64\mkl_intel_lp64.lib" "C:\Program Files (x86)\Intel\Composer XE\mkl\lib\intel64\mkl_core.lib" "C:\Program Files (x86)\Intel\Composer XE\mkl\lib\intel64\mkl_intel_thread.lib" -I"C:\Program Files (x86)\Intel\Composer XE\mkl\include" -I"C:\Program Files (x86)\Intel\Composer XE 2013 SP1\compiler\include" COMPFLAGS="$COMPFLAGS /Qopenmp /Qansi_alias"
mex trans_occ.cpp  -DUSE_OPENMP "C:\Program Files (x86)\Intel\Composer XE\mkl\lib\intel64\mkl_intel_lp64.lib" "C:\Program Files (x86)\Intel\Composer XE\mkl\lib\intel64\mkl_core.lib" "C:\Program Files (x86)\Intel\Composer XE\mkl\lib\intel64\mkl_intel_thread.lib" -I"C:\Program Files (x86)\Intel\Composer XE\mkl\include" -I"C:\Program Files (x86)\Intel\Composer XE 2013 SP1\compiler\include" COMPFLAGS="$COMPFLAGS /Qopenmp /Qansi_alias"
%{
mex update_dist.cpp -DUSE_OPENMP COMPFLAGS="$COMPFLAGS /Qopenmp /Qansi_alias"
mex trans_occ.cpp
make_dpopt('cutoff');
%}
% mex update_dist.cpp -g COMPFLAGS="$COMPFLAGS /Qopenmp /Qansi_alias"
% end

if ~isequal(exist('dpopt_worker','file'),3)
% make_flag.donlp2 = 'on';
% make_flag.snopt = 'on';
make_flag.MaxData = 30;
make_flag.MaxDim = 10;
make_dpopt('worker', make_flag);
end

if ~isequal(exist('dpopt_manager','file'),3)
% make_flag.donlp2 = 'on';
% make_flag.snopt = 'on';
make_flag.MaxData = 30;
make_flag.MaxDim = 10;
make_dpopt('manager', make_flag);
end
%}

% if ~isequal(exist('dpopt_manager_c','file'),3)
% make_flag.donlp2 = 'on';
% make_flag.MaxData = 30;
% make_flag.MaxDim = 10;
% make_dpopt('manager_c', make_flag);
% end

end
function gen_mex(mexfile, debug)

validateattributes(mexfile, {'char'}, {'row'});

%# defines
defs = ' -largeArrayDims -DMEX ';

if nargin == 2 && ~isempty(debug)
    validateattributes(debug, {'numeric', 'logical'}, {'scalar'});
    if debug, defs = [defs '-g -DDEBUG']; end
end

%# paths
opts.defs = { ...
    };

opts.incpaths = { ...
    ... 'UNCOMMENT AND SET PATH TO EIGEN', ...
    '../include' ...
    };

opts.libpaths = { ...
    };

opts.libs = { ...
    };

[pathstr, ~, ~] = fileparts(mfilename('fullpath'));
cmd = ['mex ' pathstr filesep mexfile ' -outdir ' pathstr defs];

defs = strcat(' -D', opts.defs);
incp = strcat(' -I', opts.incpaths);
libp = strcat(' -L', opts.libpaths);
libs = strcat(' -l', opts.libs);
cmd = [cmd, defs{:}, incp{:}, libp{:}, libs{:}];

clc
clear mex
eval(cmd)

end
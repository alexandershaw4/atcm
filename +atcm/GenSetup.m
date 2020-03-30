fprintf(['BEWARE: spawning will overwrite existing' ...
    ' Setup.m\nin current directory, if it exists!\n\n']);

OK = input('Spawn new Setup script? (y/n)','s');

if strcmp(OK,'y')

    [p,f,e] = fileparts( mfilename('fullpath') );

    orig_cfg = [p '/Setup.m'];

    curr_dir = evalinContext('pwd');

    new_file = [curr_dir '/Setup.m']; 

    str = ['cp ' orig_cfg ' ' new_file];

    unix( str );

    edit(new_file);
end
function addDirectories
% addDirectories - Set the path and add required functions to the path.
    %
    % addDirectories is a handy script which: 
    % 1. It adds the necessary directories, such as rawData to the
    % MATLAB path for the current user.
    % 2. It runs latex_interpreter function once, which helps make plot
    % text like LaTeX.
    % 
    % Usage:
    % addDirectories
    % 
    % Example:
    % addDirectories
    % 
    % See also: cd, addpath
    
    % Check the current user's username
   
    
    % Add relevant directories to the MATLAB path
    addpath('rawData\')
    latex_interpreter
end
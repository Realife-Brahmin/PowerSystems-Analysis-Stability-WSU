%% EE 521 Power System Analysis
% Preamble and Control Inputs

start = tic;
clearVars = true;
localUsername = getenv('USERNAME');
listOfUsernames = {'aryan', 'Aryan Ritwajeet Jha'};
if ismember(localUsername,  listOfUsernames)
    % For user 'aryan', change the current directory to the specified path
    cd(strcat("C:", filesep, "Users", filesep, localUsername, filesep, "Documents", filesep, ...
        "documents_general", filesep, "PowerSystems-Analysis-Stability-WSU", filesep) )
    addpath(genpath('functions\'))
else
    fprintf("Are you not me? Might want to add the folder to the path or add folder to the workspace.\n");
end
addDirectories;

clearVariables(clearVars);

% systemName = "PhD-QE-CaseA-3"
% systemName = "ieee14"
% systemName = "anamika3A"
systemName = "anamika3B"
powerFlowMethod =  "NRPF" 
% powerFlowMethod =  "Decoupled NRPF"
% powerFlowMethod =  "Fast Decoupled NRPF"
doContinuationPowerFlow = false;
displayCPFResults = true;
plotCPFPlots = true;
CPF_Bus = 4;
useSparseDSA = false;
includeOPFScenarios = false;
showOPFFormulae = false;
showOPFValues = false;
numIterations = 50; %I don't wait for the system to converge, 
printPowerFlowConvergenceMessages = false;
% neither do I care if the system converges earlier.
toleranceLimit = 1e-3; %mean of absolute values of 
% corrections should be less than this for convergence to be achieved.
displayRawData = false;
displayYBus = true;
displayTables = true; %show busData, branchData, ybus, 
% basically data structures which are not the final output.
printJacobians = true  ; %Print Jacobians during NRPF iterations? Does not work if displayTables is off.
displayLUFactors = true;
printMismatches = true; %Print Mismatches during NRPF iterations? Does not work if displayTables is off.
printCorrections = true;
disableTaps = false; %Disable Tap-changers when commputing YBus?
showPlots = false;
% showPlots = true;
displayResults = true;
reducedBranchColumnsCDFReading = true; 
saveBusDataAndBranchData = false;
showImages = true; %might add iteration specific images later. 
verboseCDFReading = false; %Will give a verbose output when reading CDF files.
MVAb = 100; %Currently the same for all systems in database.
% 

folder_rawData = "rawData/"; %location of CDF .txt file for the system
file_rawData = strcat(folder_rawData, systemName, "cdf.txt"); %Exact location of CDF .txt file for the system
folder_processedData = "processedData/";
% Should configure it to be read from the CDF file later.
latex_interpreter %for LaTeX typesetting in plots

if contains(systemName, 'ieee11')
    systemName1 = 'ieee11';
else
    systemName1 = systemName;
end
% Read CDF file and store the data in neat MATLAB |tables: busData| and |branchData|.

[busData, branchData, N, numBranch] = ...
    readCDF(file_rawData, reducedBranchColumnsCDFReading, verboseCDFReading);
if displayTables && displayRawData
    displayRawDataAsTables(busData, branchData, N, numBranch);
end

if saveBusDataAndBranchData
    fileType = '.csv';
    filenameBusData = strcat(folder_processedData, systemName1, "/busData", fileType);
    writetable(busData, filenameBusData);
    filenameBranchData = strcat(folder_processedData, systemName1, "/branchData", fileType);
    writetable(branchData, filenameBranchData);
end
% Extract $Y_{Bus}$,  Adjacency List $E$ from the |branchData table.|

if useSparseDSA && ~doContinuationPowerFlow && ~includeOPFScenarios
    [nnzYBus, NYBus] = makeSparseYBus(busData, branchData, displayTables, displayYBus);
else
    [ybus, BMatrix, ~, ~, ~, E] = ybusGenerator(busData, branchData, 'verbose', true);
    % [ybusByType, BMatrixByType, ~, ~, ~, ENew] = ybusGenerator(busData, branchData, 'disableTaps', false, 'sortBy', "busTypes", "verbose", true, 'saveTables', true, 'saveLocation', folder_processedData, 'systemName', systemName1);
    % [ybusByType, BMatrixByType, ~, ~, ~, ENew] = ybusGenerator(busData, branchData, 'disableTaps', false, "verbose", true, 'saveTables', true, 'saveLocation', folder_processedData, 'systemName', systemName1);

end
% Run Newton Raphson Power Flow and obtain a steady state snapshot of the system variables $P_i, Q_i, V_i, \delta_i$ $\forall$ buses $i \in [1, N], i \in \mathbb{N}$

[PSpecified, QSpecified, V, delta, ...
    listOfPQBuses, listOfPVBuses, nPQ, nPV, ...
    listOfNonSlackBuses] = initializeVectors(busData);

if useSparseDSA && ~doContinuationPowerFlow && ~includeOPFScenarios   
    doTheSparseThing(PSpecified, QSpecified, V, delta, nnzYBus, NYBus, busData); 
else
    if contains(systemName, 'caseThree')
        resultsFromCaseTwo = load("processedData\ieee11-caseTwoResults");
        resultsFromCaseTwo = resultsFromCaseTwo.resultTable;
        wiggle = 0.35; %minimum 0.35 value for good result
        V = V*(1-wiggle) + wiggle*resultsFromCaseTwo.V;
        delta = delta*(1-wiggle) + wiggle*resultsFromCaseTwo.delta;
    end
    
    if doContinuationPowerFlow && ~includeOPFScenarios
        desiredOutput = 'both';
        [V_CPF, delta_CPF, lambda_CPF, iter_CPF, sigma_CPF] =...
        continuationPowerFlow(busData, ybus, PSpecified, QSpecified, V, delta, BMatrix, E, ...
        nPQ, nPV, listOfPQBuses, listOfNonSlackBuses, powerFlowMethod, desiredOutput, ...
        CPF_Bus, displayCPFResults, plotCPFPlots);
    else
        [P, Q, V, delta] = solveForPowerFlow(PSpecified, QSpecified, V, delta, ...
            ybus, BMatrix, E, nPQ, nPV, ...
            listOfPQBuses, listOfNonSlackBuses, ...
            numIterations, toleranceLimit, powerFlowMethod, displayTables, ...
            printJacobians, printMismatches, printCorrections, displayLUFactors, ...
            printPowerFlowConvergenceMessages);
        resultTable = displayPowerFlowResults(P, Q, V, delta, displayResults);
        
        if contains(systemName, 'caseTwo') && strcmp(powerFlowMethod, 'NRPF')
            save("processedData\ieee11-caseTwoResults", "resultTable");
        end
    end
end
% Compare obtained snapshot values of $V_i$ and $\delta_i$ against the ones given in the CDF file.

if ~doContinuationPowerFlow
    plotPowerFlowResults(showPlots, V, busData, systemName, powerFlowMethod, delta);
end
% Economic Dispatch and Optimal Power Flow Calculations:

if includeOPFScenarios && strcmp(systemName, 'ieee14') && ~useSparseDSA
    runOPFScenarios(busData, P, Q, V, delta, N, ybus, BMatrix, E, nPQ, nPV, listOfPQBuses, listOfNonSlackBuses, powerFlowMethod, MVAb, showOPFFormulae, showOPFValues, printPowerFlowConvergenceMessages)
end
toc(start)
%% Have a nice day! 
% 
% 
% In case you encounter a Java Heap Memory error, delete the above gif, or go 
% to |Preferences->General->Java Heap Memory| and increase the allocated size.
%%
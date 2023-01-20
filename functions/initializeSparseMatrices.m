function [LSparse, USparse] = initializeSparseMatrices(JSparse)
    [LSparse, USparse] = factorizeLU(JSparse);
end
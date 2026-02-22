function dataDict = csvToDictionary(csvFile, fromIdx, toIdx)
%CSVTODICTIONARY Convert a CSV file (or a row range) into a MATLAB dictionary.
%   dataDict = csvToDictionary(csvFile) imports all data rows.
%   dataDict = csvToDictionary(csvFile, fromIdx, toIdx) imports only rows
%   from fromIdx to toIdx (1-based, relative to data rows after header).
%
%   - The first CSV row is used as keys.
%   - Each value is stored as a column vector (N x 1).
%
%   Note:
%   MATLAB dictionary objects typically require cell-wrapped values for
%   variable-size per-key arrays. This function returns a containers.Map
%   so each lookup returns the column vector directly.

    if nargin < 2 || isempty(fromIdx)
        fromIdx = 1;
    end
    if nargin < 3 || isempty(toIdx)
        toIdx = inf;
    end

    if ~(ischar(csvFile) || isstring(csvFile))
        error('csvToDictionary:InvalidPath', 'csvFile must be a character vector or string.');
    end
    if ~isfile(csvFile)
        error('csvToDictionary:FileNotFound', 'CSV file not found: %s', string(csvFile));
    end

    validateattributes(fromIdx, {'numeric'}, {'scalar', 'integer', '>=', 1}, mfilename, 'fromIdx');
    if ~isinf(toIdx)
        validateattributes(toIdx, {'numeric'}, {'scalar', 'integer', '>=', fromIdx}, mfilename, 'toIdx');
    end

    opts = detectImportOptions(csvFile, 'NumHeaderLines', 0);
    opts.VariableNamingRule = 'preserve';

    startRowInFile = fromIdx + 1;
    if isinf(toIdx)
        opts.DataLines = [startRowInFile inf];
    else
        endRowInFile = toIdx + 1;
        opts.DataLines = [startRowInFile endRowInFile];
    end

    T = readtable(csvFile, opts);
    keys = string(T.Properties.VariableNames);

    dataDict = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for col = 1:numel(keys)
        colData = T.(keys(col));
        dataDict(char(keys(col))) = colData(:);
    end
end

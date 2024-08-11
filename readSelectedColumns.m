function tbl = readSelectedColumns(filePath)
    % Create import options
    opts = detectImportOptions(filePath, 'VariableNamingRule', 'preserve');
    
    % Specify which columns to import (2nd and 4th)
    opts.SelectedVariableNames = opts.VariableNames([2, 3]);
    
    % Read the table with the specified options
    tbl = readtable(filePath, opts);
end

function tableDict = readCsvTables(fileName)
    % READCSVTABLES Reads a CSV file with multiple semicolon-delimited tables
    % Input: fileName - Path to the CSV file
    % Output: tableDict - Struct where each field is a table name containing a struct of variable-value arrays
    % Compatible with MATLAB R2018b
    % Tables are terminated by: 1) blank line, 2) line starting with **, 3) line starting with ;
    
    % Initialize the output dictionary
    tableDict = struct();
    
    % Open the file
    fid = fopen(fileName, 'r');
    if fid == -1
        error('Cannot open file: %s', fileName);
    end
    
    % Initialize variables
    currentTable = '';
    headers = {};
    lineNum = 0;
    dataRows = {};
    expectData = false;
    
    % Read the file line by line
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        lineNum = lineNum + 1;
        
        % Check for table termination conditions
        if isempty(line) || strncmp(line, '**', 2) || strncmp(line, ';', 1)
            % If we have a table with headers, process accumulated data
            if ~isempty(currentTable) && ~isempty(headers)
                if isempty(dataRows)
                    warning('Line %d: Table %s has no data rows.', lineNum, currentTable);
                    tableDict.(currentTable) = struct(); % Store empty struct
                else
                    % Initialize table struct
                    tableVars = struct();
                    for j = 1:length(headers)
                        varName = matlab.lang.makeValidName(headers{j});
                        % Initialize array for this variable
                        values = cell(length(dataRows), 1);
                        
                        % Process each data row for this variable
                        for r = 1:length(dataRows)
                            tokens = dataRows{r};
                            if j <= length(tokens)
                                value = tokens{j};
                                % Convert value to appropriate type
                                if strcmpi(value, 'true')
                                    values{r} = true;
                                elseif strcmpi(value, 'false')
                                    values{r} = false;
                                else
                                    % Try to convert to numeric
                                    numValue = str2double(value);
                                    if ~isnan(numValue)
                                        values{r} = numValue;
                                    else
                                        values{r} = value; % Store as string
                                    end
                                end
                            else
                                values{r} = NaN; % Handle missing values
                            end
                        end
                        
                        % Convert cell array to appropriate type
                        if all(cellfun(@islogical, values))
                            % Convert to logical array or scalar
                            if length(values) == 1
                                tableVars.(varName) = values{1}; % Scalar for single row
                            else
                                tableVars.(varName) = cell2mat(values); % Logical array
                            end
                        elseif all(cellfun(@isnumeric, values) & ~cellfun(@isnan, values))
                            % Convert to numeric array
                            tableVars.(varName) = cell2mat(values);
                        else
                            % Keep as cell array for mixed or string types
                            if length(values) == 1
                                tableVars.(varName) = values{1}; % Scalar for single row
                            else
                                tableVars.(varName) = values;
                            end
                        end
                    end
                    tableDict.(currentTable) = tableVars;
                end
            end
            
            % If the line starts with **, start a new table
            if strncmp(line, '**', 2)
                currentTable = matlab.lang.makeValidName(line(3:end)); % Remove ** and make valid field name
                headers = {};
                dataRows = {};
                expectData = false;
            else
                % Reset for blank line or comment
                currentTable = '';
                headers = {};
                dataRows = {};
                expectData = false;
            end
            continue;
        end
        
        % If we have a current table, process headers or data
        if ~isempty(currentTable)
            % Split the line by semicolon
            tokens = split(line, ';');
            tokens = strtrim(tokens); % Remove leading/trailing whitespace
            
            % If headers are not set, this is the header row
            if isempty(headers)
                headers = tokens;
                expectData = true; % Expect data rows next
            elseif expectData
                % Validate data row
                if length(tokens) ~= length(headers)
                    warning('Line %d: Number of values (%d) does not match headers (%d) in table %s. Skipping row.', ...
                        lineNum, length(tokens), length(headers), currentTable);
                    continue;
                end
                % Store data row
                dataRows{end+1} = tokens;
            end
        end
    end
    
    % Handle case where file ends with an open table
    if ~isempty(currentTable) && ~isempty(headers)
        if isempty(dataRows)
            warning('File ended with incomplete table %s: No data rows found.', currentTable);
            tableDict.(currentTable) = struct();
        else
            % Process accumulated data rows
            tableVars = struct();
            for j = 1:length(headers)
                varName = matlab.lang.makeValidName(headers{j});
                values = cell(length(dataRows), 1);
                for r = 1:length(dataRows)
                    tokens = dataRows{r};
                    if j <= length(tokens)
                        value = tokens{j};
                        if strcmpi(value, 'true')
                            values{r} = true;
                        elseif strcmpi(value, 'false')
                            values{r} = false;
                        else
                            numValue = str2double(value);
                            if ~isnan(numValue)
                                values{r} = numValue;
                            else
                                values{r} = value;
                            end
                        end
                    else
                        values{r} = NaN;
                    end
                end
                if all(cellfun(@islogical, values))
                    if length(values) == 1
                        tableVars.(varName) = values{1};
                    else
                        tableVars.(varName) = cell2mat(values);
                    end
                elseif all(cellfun(@isnumeric, values) & ~cellfun(@isnan, values))
                    tableVars.(varName) = cell2mat(values);
                else
                    if length(values) == 1
                        tableVars.(varName) = values{1};
                    else
                        tableVars.(varName) = values;
                    end
                end
            end
            tableDict.(currentTable) = tableVars;
        end
    end
    
    % Close the file
    fclose(fid);
    
    % Warn if no tables were processed
    if isempty(fieldnames(tableDict))
        warning('No valid tables were processed from %s.', fileName);
    end
end

function parsedMPR = parseBiologicMPR(filename)
%PARSEBIOLOGICMPR  Parse a Bio-Logic .mpr file that has ASCII "MODULE" markers.
%
% This code:
%   1) Checks for a "BIO-LOGIC MODULAR FILE" magic string.
%   2) Reads the entire remainder of the file as bytes.
%   3) Splits by ASCII "MODULE".
%   4) For each module chunk, parses the module header (short_name, long_name, etc.).
%   5) Dispatches to sub-functions: parseVMPSet, parseVMPData, parseVMPLog, parseVMPLoop, etc.
%
% NOTE:
%   - You must adapt offsets, data types, and endianness to match your actual .mpr structure.
%   - The parse*Module functions provided here are stubs. Insert real logic for each module type.
%
% Usage:
%   dataStruct = parseBiologicMPR('myFile.mpr')

    fid = fopen(filename, 'rb');
    if fid < 0
        error('Could not open file: %s', filename);
    end

    % 1) Check magic (roughly 64 bytes from typical exports, your version may differ)
    magicBytes = fread(fid, [1,28], '*char');
    magicStr = deblank(char(magicBytes));
    expectedMagic = 'BIO-LOGIC MODULAR FILE';
    if ~contains(magicStr, expectedMagic)
        fclose(fid);
        error('Not a recognized .mpr file. Expected magic: %s', expectedMagic);
    end

    % 2) Read the remainder of the file as raw bytes
    remainder = fread(fid, '*uint8')';
    fclose(fid);

    % Convert to char for naive splitting by "MODULE"
    textAll = char(remainder);
    parts = strsplit(textAll, 'MODULE');
    
    % The first chunk (parts{1}) is usually leftover text after the file magic block
    % or possibly empty. The subsequent chunks are each "modules."
    
    % Initialize output structure
    parsedMPR = struct();
    parsedMPR.Modules = {};
    %disp(length(parts))
    for iModule = 2:length(parts)
        rawModuleText = parts{iModule};
        rawModule = uint8(rawModuleText); % Convert back to bytes for binary reading
                
        % ===========================================================
        % Module Header Formats in Biologic VMP3 .mpr Files
        % ===========================================================
        % The module headers can follow two different formats:
        %
        % 1. **Newer Header Format (59 bytes total)**
        %    - Includes `max_length` and `newver` fields.
        %    - Used in newer versions where `newver` helps determine the file version.
        %    - Field layout:
        %      | Offset | Field        | Type    | Size (bytes) |
        %      |--------|-------------|---------|-------------|
        %      | 0      | short_name  | char[10]| 10          |
        %      | 10     | long_name   | char[25]| 25          |
        %      | 35     | max_length  | uint32  | 4           |
        %      | 39     | length      | uint32  | 4           |
        %      | 43     | oldver      | uint32  | 4           |
        %      | 47     | newver      | uint32  | 4           |
        %      | 51     | date        | char[8] | 8           |
        %      |--------|-------------|---------|-------------|
        %      | Total  |             |         | **59 bytes** |
        %
        % 2. **Older Header Format (51 bytes total)**
        %    - Does **not** include `max_length` or `newver`.
        %    - `oldver` alone is used to determine the version.
        %    - Field layout:
        %      | Offset | Field        | Type    | Size (bytes) |
        %      |--------|-------------|---------|-------------|
        %      | 0      | short_name  | char[10]| 10          |
        %      | 10     | long_name   | char[25]| 25          |
        %      | 35     | length      | uint32  | 4           |
        %      | 39     | oldver      | uint32  | 4           |
        %      | 43     | date        | char[8] | 8           |
        %      |--------|-------------|---------|-------------|
        %      | Total  |             |         | **51 bytes** |
                % Ensure we have enough bytes for header parsing


        HEADER_LENGTH_NEW = 59;
        HEADER_LENGTH_OLD = 51;


        if numel(rawModule) < min(HEADER_LENGTH_OLD, HEADER_LENGTH_NEW)
            warning('Module chunk %d is too short to parse header fully.', iModule-1);
            continue;
        end

        % Extract module header fields
        shortName  = char(rawModule(1:10));
        longName   = char(rawModule(11:35));
    
        lengthBytesOld = typecast(rawModule(36:39), 'uint32'); % Read length
        lengthBytesNew = typecast(rawModule(40:43), 'uint32'); % Read length

        if (length(rawModule) == lengthBytesOld + HEADER_LENGTH_OLD)
            debug_fprintf("Module matches old header type")
            lengthVal = double(lengthBytesOld);
            oldverBytes = typecast(rawModule(40:43), 'uint32'); % Read old version
            newverBytes = 0; % just set to zero if it doesn't exist
            version = double(oldverBytes);
            dateBytes   = char(rawModule(44:51)); % Read date (if applicable)
            % The actual module data portion is everything after these 55 bytes 
            moduleData = rawModule(52:end); 

        elseif (length(rawModule) == lengthBytesNew + HEADER_LENGTH_NEW)
            debug_fprintf("Module matches new header type!")
            lengthVal = double(lengthBytesNew);
            oldverBytes = typecast(rawModule(44:47), 'uint32'); % Read old version
            %disp(oldverBytes)

            newverBytes = typecast(rawModule(48:51), 'uint32'); % Read new version
            %disp(newverBytes)
            version = double(newverBytes);
            dateBytes   = char(rawModule(52:59)); % Read date (if applicable)
            moduleData = rawModule(60:end);             % The actual module data portion is everything after these 55 bytes 

        else
             disp("Error: Module does not match either old or new header type")
             error('Header module length does not match expected length')
        end

        % Determine minver based on version
        if version >= 10
            minver = "11.50";
        else
            minver = "10.40";
        end
    
        % Store in structure
        modHeader = struct();
        modHeader.short_name = strtrim(shortName);
        modHeader.long_name  = strtrim(longName);
        modHeader.length     = lengthVal;
        modHeader.oldver     = double(oldverBytes);
        modHeader.newver     = double(newverBytes);
        modHeader.version    = version;
        modHeader.date       = strtrim(dateBytes);
        modHeader.minver     = minver;
    
        debug_fprintf("Read '%s' with version '%d' ('%s')\n", modHeader.short_name, version, minver);


 
        %disp('First 32 bytes of moduleData (hex):'); 
        %disp(dec2hex(moduleData(1:32)));  % Print hex values for debugging 
        

        % 3) Dispatch based on short_name or long_name
        % We'll do simple string checks:
        switch modHeader.short_name
            case 'VMP Set'
                dataStruct = parseVMPSet(moduleData, modHeader);
            case 'VMP data'
                dataStruct = parseVMPData(moduleData, modHeader);
            case 'VMP LOG'
                %dataStruct = parseVMPLog(moduleData, modHeader);
            case 'VMP loop'
                %dataStruct = parseVMPLoop(moduleData, modHeader);
            otherwise
                % Unknown or not implemented
                warning('Module "%s" not recognized. Storing raw bytes.', modHeader.short_name);
                dataStruct = struct('rawModuleData', moduleData);
        end

        % Accumulate results
        parsedMPR.Modules{end+1} = struct( ...
            'Header', modHeader, ...
            'Data',   dataStruct ...
        );
    end
end

%% Below are STUBS for each module type. Fill these in with real offsets & logic.

function out = parseVMPSet(moduleData, modHeader)
%PARSEVMPSET  Parse the "VMP Set" (Settings) module data
%
%  This code mimics "process_settings" from Python:
%   - Byte 0: technique_id
%   - Then reads fields from a settings_dtypes map
%   - Potentially scans for param sets at offsets 0x0572, 0x1845, etc.
%
%  Inputs:
%    moduleData - the raw bytes of this module, after the 5-field header
%    modHeader  - struct with oldver/newver, date, etc.
%
%  Output:
%    out.Technique  - e.g. "GCPL", "CV", "CA", etc. (or raw ID if unknown)
%    out.Settings   - struct of fields like "comments", "electrode_area", etc.
%    out.Params     - array or struct of advanced technique parameters (if found)

    out = struct();
    % 1) Technique ID is typically at offset 0x0000
    if isempty(moduleData)
        warning('VMP Set moduleData is empty.');
        out.Technique  = 'Unknown';
        out.Settings   = struct();
        out.Params     = [];
        return;
    end

    techniqueByte = moduleData(1);  % 1-based indexing in MATLAB
    out.Technique = localLookupTechnique(techniqueByte);

    % 2) Read known settings fields (like offset 0x0007 => "comments", etc.)
    %    This map corresponds to your 'settings_dtypes' in Python.
    %    Key = decimal offset, Value = {dtypeString, fieldName}.
    settingsMap = containers.Map('KeyType','uint32','ValueType','any');
    % Example entries (adapt to your actual offsets from the Python code):
    settingsMap(uint32(7))    = {"pascal", "comments"};
    settingsMap(uint32(hex2dec('0107'))) = {"<f4", "active_material_mass"};
    settingsMap(uint32(hex2dec('010B'))) = {"<f4", "at_x"};
    settingsMap(uint32(hex2dec('010F'))) = {"<f4", "molecular_weight"};
    settingsMap(uint32(hex2dec('0113'))) = {"<f4", "atomic_weight"};
    settingsMap(uint32(hex2dec('0117'))) = {"<f4", "acquisition_start"};
    settingsMap(uint32(hex2dec('011B'))) = {"<u2", "e_transferred"};
    settingsMap(uint32(hex2dec('011E'))) = {"pascal", "electrode_material"};
    % etc...

    out.Settings = struct();
    keys = settingsMap.keys;
    for i = 1:numel(keys)
        offsetKey = keys{i};
        infoCell  = settingsMap(offsetKey);
        dtypeStr  = infoCell{1};
        fieldName = infoCell{2};

        % Ensure offset is within moduleData
        endPos = offsetKey + localByteSize(dtypeStr) - 1;  % approximate
        if endPos > numel(moduleData)
            % Not enough bytes
            continue;
        end

        rawVal = moduleData(offsetKey+1 : endPos); % 1-based indexing
        val = readBioLogicValue(rawVal, dtypeStr);
        out.Settings.(fieldName) = val;
    end

    % 3) Potentially parse technique parameters at offsets like 0x1845, 0x1846, ...
    out.Params = localParseTechniqueParams(moduleData, out.Technique, modHeader);

end

%% localLookupTechnique: Convert technique ID byte to a technique name
%% localLookupTechnique: Convert technique ID byte to a name
function techName = localLookupTechnique(byteVal)
    % Define technique ID-to-name mapping
    techniqueKeys = [...
        0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, ...
        0x0B, 0x0C, 0x0D, 0x0E, 0x0F, 0x10, 0x11, 0x12, 0x18 ...
    ];

    techniqueValues = {...
        'OCV',                  % 0x01 - Open Circuit Voltage
        'PEIS',                 % 0x02 - Potentio Electrochemical Impedance Spectroscopy
        'CA',                   % 0x03 - Chronoamperometry
        'GCPL',                 % 0x04 - Galvanostatic Cycling w/ Potential Limitation
        'CP',                   % 0x05 - Chronopotentiometry
        'CV',                   % 0x06 - Cyclic Voltammetry
        'DPV',                  % 0x07 - Differential Pulse Voltammetry
        'SWV',                  % 0x08 - Square Wave Voltammetry
        'NPV',                  % 0x09 - Normal Pulse Voltammetry
        'CV Staircase',         % 0x0A - Cyclic Voltammetry (Staircase)
        'DPA',                  % 0x0B - Differential Pulse Amperometry
        'SWV Amperometry',      % 0x0C - Square Wave Amperometry
        'NP Amperometry',       % 0x0D - Normal Pulse Amperometry
        'CPVS',                 % 0x0E - Chronopotentiometry with Voltage Step
        'GITT',                 % 0x0F - Galvanostatic Intermittent Titration Technique
        'PITT',                 % 0x10 - Potentiostatic Intermittent Titration Technique
        'MPT',                  % 0x11 - Multi-Potential Steps
        'IST',                  % 0x12 - Intermittent Signal Test
        'CA'                    % 0x18 - Chronoamperometry (Alternative)
    };

    % Create the mapping
    techniqueMap = containers.Map(techniqueKeys, techniqueValues);

    % Lookup or fallback
    if isKey(techniqueMap, byteVal)
        techName = techniqueMap(byteVal);
        %disp("technique name: ")
        %disp(byteVal)
        %disp(techName)
    else
        techName = sprintf('Unknown (0x%02X)', byteVal);
    end
end

%% localParseTechniqueParams: Read parameters based on known offsets
function paramStruct = localParseTechniqueParams(moduleData, techniqueName, modHeader)
    paramStruct = struct(); % Initialize output

    % Search for `ns` dynamically
    %base_offset_candidates = [hex2dec('1843'), hex2dec('1844'), hex2dec('1845'), hex2dec('1846')];
    ns_offset = hex2dec('1845');
    %ns_offset = hex2dec('1846');


    raw_ns = moduleData(ns_offset : ns_offset + 1);

    ns = typecast(raw_ns, 'uint16');
    swapped_ns = swapbytes(ns);
    ns = min(ns, swapped_ns); % Take the more reasonable value
    debug_fprintf('Final chosen ns_offset = 0x%X: ns=%d\n', ns_offset, ns);

    % Validate `ns`
    % if ns > 0 && ns < 100
    %     found_ns = true;
    % else 
    %     warning('invalid ns')
    % end



    %fprintf('Chosen ns_offset=0x%X: ns=%d\n', ns_offset, ns);

    % `n_params` is always **2 bytes after** `ns`
    np_offset = ns_offset + 2;
    raw_np = moduleData(np_offset : np_offset + 1);
    n_params = typecast(raw_np, 'uint16');
    swapped_np = swapbytes(n_params);
    
    % If swapped version makes more sense, use it
    if swapped_np > 0 && swapped_np < 100
        n_params = swapped_np;
    end

    debug_fprintf('Final chosen n_params offset = 0x%X: n_params=%d\n', np_offset, n_params);
    
    % Store results
    paramStruct.ns = ns;
    paramStruct.n_params = n_params;
end





%% readBioLogicValue: Interpret raw bytes given a dtype string
function val = readBioLogicValue(rawVal, dtypeStr)
    % Ensure rawVal is a column vector
    rawVal = rawVal(:);
    numBytes = numel(rawVal);

    % Handle empty input case
    if isempty(rawVal)
        val = [];
        return;
    end

    % Map Python dtype strings to MATLAB equivalents
    dtypeMap = containers.Map( ...
        {'<f4', '<f8', '<u2', '<u4', '<i2', '<i4'}, ...
        {'single', 'double', 'uint16', 'uint32', 'int16', 'int32'} ...
    );

    % Handle Pascal strings separately
    if strcmp(dtypeStr, 'pascal')
        if numBytes == 0
            val = '';
            return;
        end
        strLen = rawVal(1);
        if strLen > 0 && (1 + strLen) <= numBytes
            val = native2unicode(rawVal(2 : 1 + strLen), 'UTF-8');
        else
            val = '';
        end
        return;
    end

    % Default: treat as ASCII if dtype is unknown
    if ~isKey(dtypeMap, dtypeStr)
        val = strtrim(char(rawVal'));
        return;
    end

    % Get MATLAB type
    dtypeMatlab = dtypeMap(dtypeStr);

    % Determine required byte size
    bytesRequired = 0;
    switch dtypeMatlab
        case {'uint16', 'int16'}, bytesRequired = 2;
        case {'uint32', 'int32', 'single'}, bytesRequired = 4;
        case {'uint64', 'int64', 'double'}, bytesRequired = 8;
    end

    % Ensure rawVal has the correct length (pad with zeros if needed)
    padSize = mod(-numBytes, bytesRequired);
    if padSize > 0
        rawVal = [rawVal; zeros(padSize, 1, 'uint8')];
    end

    % Convert using typecast
    val = typecast(rawVal, dtypeMatlab);
end


%% localByteSize: approximate # of bytes needed to read one value of dtype
function n = localByteSize(dtypeStr)
    switch dtypeStr
        case '<f4'
            n = 4;
        case '<f8'
            n = 8;
        case '<u2'
            n = 2;
        case '<u4'
            n = 4;
        case 'pascal'
            % Hard to guess exact size; might read up to 256 bytes if the length byte is big
            % We'll assume at least 2 for now
            n = 16;  % a small guess: your file might store bigger
        otherwise
            n = 1;
    end
end


function out = parseVMPData(moduleData, modHeader)
%PARSEVMPDATA - Parses "VMP data" from Biologic MPR files row-by-row.
    % some code from Python started indices at 0 whereas matlab starts arrays at 1
    CONSTANT_ZERO_TO_ONE_ARRAY_INDEX = 1;
    
    out = struct();
    out.DataPoints  = [];
    out.Columns     = {};
    out.ColumnIDs   = [];

    % Step 1: Read header info
    %disp('moduleData: ')
    %disp(moduleData(1:10))
    
    nDataPoints = typecast(moduleData((0:3) + CONSTANT_ZERO_TO_ONE_ARRAY_INDEX), 'uint32');
    nColumns    = double(moduleData(4+CONSTANT_ZERO_TO_ONE_ARRAY_INDEX));
    debug_fprintf("nDataPoints: %d\n", nDataPoints)
    debug_fprintf("nColumns: %d\n", nColumns)


    if nColumns == 0 || nColumns > 100
        warning('[parseVMPData] Failed to detect valid nColumns.');
        return;
    end

    % Step 2: Extract column IDs
    colIDStart = 5;
    colIDEnd = colIDStart + (2 * nColumns) - 1;

    if colIDEnd > length(moduleData)
        warning('[parseVMPData] Not enough bytes for column IDs.');
        return;
    end

    %disp(modHeader)


    if ismember(modHeader.version, [10, 11])
        % Read big-endian (">u2" in NumPy)
        columnIDs = typecast(swapbytes(typecast(moduleData(6:6 + 2 * nColumns - 1), 'uint16')), 'uint16');
    elseif ismember(modHeader.version, [2, 3])
        % Read little-endian ("<u2" in NumPy)
        columnIDs = typecast(moduleData(6:6 + 2 * nColumns - 1), 'uint16');
    else 
        error('no version detected in modHeader')
    end


    %columnIDBytes = moduleData((colIDStart:colIDEnd)+CONSTANT_ZERO_TO_ONE_ARRAY_INDEX);
    %disp(columnIDBytes(1:10))

    %columnIDs = typecast(columnIDBytes, 'uint16');
    out.ColumnIDs = columnIDs;

    % Step 3: Get column names, data types, and bytes per column
    [columnNames, dataTypes, bytesPerCol, flagInfo] = lookupColumnTypes(columnIDs);

    bytesPerRow = sum(bytesPerCol);

    % Step 4: Determine data start offset

    if modHeader.version == 3
        dataStart = hex2dec('196') + CONSTANT_ZERO_TO_ONE_ARRAY_INDEX;
    elseif modHeader.version == 2
        dataStart = hex2dec('195') + CONSTANT_ZERO_TO_ONE_ARRAY_INDEX;
    elseif modHeader.version == 10 || modHeader.version == 11
        dataStart = hex2dec('3EF') + CONSTANT_ZERO_TO_ONE_ARRAY_INDEX;
    else
        disp(modHeader)
        error('no version detected'); % this probably indicates a problem, maybe not
        dataStart = hex2dec('195') + CONSTANT_ZERO_TO_ONE_ARRAY_INDEX;
    end

    debug_fprintf("Data offset: 0x%x\n", dataStart)

    bytesNeeded = nDataPoints * bytesPerRow;
    dataEnd = dataStart + bytesNeeded - 1;

    if dataEnd > length(moduleData)
        nDataPoints = floor((length(moduleData) - dataStart) / bytesPerRow);
    end

    if nDataPoints <= 0
        return;
    end

    % Step 5: Allocate Data Matrix
    dataMatrix = zeros(nDataPoints, nColumns);

    % Step 6: Read data row-by-row
    for rowIdx = 1:nDataPoints
        rowStart = dataStart + (rowIdx - 1) * bytesPerRow;
        rowBytes = moduleData(rowStart : rowStart + bytesPerRow - 1);
        nColumnsToIterate = numel(bytesPerCol);  % ✅ Corrected this line
    
        % Process each column in the row
        rowOffset = 1;
        for colIdx = 1:nColumnsToIterate
            colType = dataTypes{colIdx};
            numBytes = bytesPerCol(colIdx);
            colBytes = rowBytes(rowOffset:rowOffset + numBytes - 1);
    
            % Only print debugging info for first 10 rows
            if rowIdx <= 10
                debug_fprintf("Row %d, Column %d (%s, %d bytes): %s\n", rowIdx, colIdx, colType, numBytes, mat2str(colBytes));
            end
    
            % Convert based on type
            switch colType
                case 'single'
                    colData = typecast(uint8(colBytes), 'single');
                case 'uint8'
                    colData = typecast(uint8(colBytes), 'uint8');
                case 'uint16'
                    colData = typecast(uint8(colBytes), 'uint16');
                case 'int16'
                    colData = typecast(uint8(colBytes), 'int16');
                case 'uint32'
                    colData = typecast(uint8(colBytes), 'uint32');
                case 'int32'
                    colData = typecast(uint8(colBytes), 'int32');
                case 'double'
                    colData = typecast(uint8(colBytes), 'double');
                otherwise
                    warning("[parseVMPData] Unknown type for column %d: %s", colIdx, colType);
                    colData = nan;
            end
    
            % Store value in matrix
            dataMatrix(rowIdx, colIdx) = colData;
    
            % Advance row offset
            rowOffset = rowOffset + numBytes;
        end
    
        % Only print row summary for first 10 rows
        % if rowIdx <= 10
        %     fprintf("Row %d Extracted Data: %s\n", rowIdx, mat2str(dataMatrix(rowIdx, :)));
        % elseif rowIdx == 11
        %     fprintf("\n=== Debugging output limited to first 10 rows ===\n");
        % end
    end



    % Step 7: Handle Flags if present
    if isfield(flagInfo, 'Flags')
        flagIdx = find(strcmp(columnNames, 'Flags'));
        flagData = dataMatrix(:, flagIdx);
        unpackedFlags = zeros(nDataPoints, numel(fieldnames(flagInfo.Flags))); 

        flagNames = fieldnames(flagInfo.Flags);
        for i = 1:numel(flagNames)
            bitmask = flagInfo.Flags.(flagNames{i}).bitmask;
            unpackedFlags(:, i) = bitand(flagData, bitmask) > 0;
        end

        % Insert unpacked flag columns
        dataMatrix = [dataMatrix(:, 1:flagIdx-1), unpackedFlags, dataMatrix(:, flagIdx+1:end)];
        columnNames = [columnNames(1:flagIdx-1), flagNames', columnNames(flagIdx+1:end)];
    end

    fprintf("Processed %d rows of data\n", nDataPoints)
    % Step 8: Store final parsed data
    out.DataPoints = dataMatrix;
    out.Columns = columnNames;
end




function [columnNames, dataTypes, bytesPerCol, flagInfo] = lookupColumnTypes(columnIDs)
% LOOKUPCOLUMNTYPES - Maps Biologic column IDs to data type, name, and bitmask details.

    FLAGS_GROUP = [1,2,3,21,31,65];

    data_columns = containers.Map('KeyType', 'double', 'ValueType', 'any');
    flag_columns = containers.Map('KeyType', 'double', 'ValueType', 'any');

    % Define standard columns
    data_columns(4)   = {'double', 'time', 8};  % Time (s)
    data_columns(5)   = {'single', 'control_V_I', 4};  % Control Voltage/Current (V/mA)
    data_columns(6)   = {'single', 'Ewe', 4}; % Working Electrode Voltage (V)
    data_columns(7)   = {'double', 'dq', 8};  % Charge increment (mA·h)
    data_columns(8)   = {'single', 'I', 4};  % Current (mA)
    data_columns(9)   = {'single', 'Ece', 4}; % Counter Electrode Voltage (V)
    data_columns(11)  = {'double', '<I>', 8}; % Mean Current (mA)
    data_columns(13)  = {'double', '(Q-Qo)', 8}; % Charge Difference (mA·h)
    data_columns(16)  = {'single', 'Analog IN 1', 4}; % Analog Input 1 (V)
    data_columns(17)  = {'single', 'Analog IN 2', 4}; % Analog Input 2 (V)
    data_columns(19)  = {'single', 'control_V', 4}; % Control Voltage (V)
    data_columns(20)  = {'single', 'control_I', 4}; % Control Current (mA)
    data_columns(23)  = {'double', 'dQ', 8}; % Charge increment (mA·h)
    data_columns(24)  = {'double', 'cycle number', 8}; % Cycle Number
    data_columns(32)  = {'single', 'freq', 4}; % Frequency (Hz)
    data_columns(33)  = {'single', '|Ewe|', 4}; % Absolute Working Electrode Voltage (V)
    data_columns(34)  = {'single', '|I|', 4}; % Absolute Current (A)
    data_columns(35)  = {'single', 'Phase(Z)', 4}; % Phase Angle (deg)
    data_columns(36)  = {'single', '|Z|', 4}; % Impedance Magnitude (Ω)
    data_columns(37)  = {'single', 'Re(Z)', 4}; % Real Part of Impedance (Ω)
    data_columns(38)  = {'single', '-Im(Z)', 4}; % Imaginary Part of Impedance (Ω)
    data_columns(39)  = {'uint16', 'I_Range', 2}; % Current Range
    data_columns(69)  = {'single', 'R', 4}; % Resistance (Ω)
    data_columns(70)  = {'single', 'P', 4}; % Power (W)
    data_columns(74)  = {'double', '|Energy|', 8}; % Absolute Energy (W·h)
    data_columns(75)  = {'single', 'Analog OUT', 4}; % Analog Output (V)
    data_columns(76)  = {'single', '<I>', 4}; % Mean Current (mA)
    data_columns(77)  = {'single', '<Ewe>', 4}; % Mean Working Electrode Voltage (V)
    data_columns(78)  = {'single', 'Cs⁻²', 4}; % Capacitance Inverse (µF⁻²)
    data_columns(96)  = {'single', '|Ece|', 4}; % Absolute Counter Electrode Voltage (V)
    data_columns(98)  = {'single', 'Phase(Zce)', 4}; % Phase of Zce (deg)
    data_columns(99)  = {'single', '|Zce|', 4}; % Impedance of Zce (Ω)
    data_columns(100) = {'single', 'Re(Zce)', 4}; % Real Part of Zce (Ω)
    data_columns(101) = {'single', '-Im(Zce)', 4}; % Imaginary Part of Zce (Ω)
    data_columns(123) = {'double', 'Energy_charge', 8}; % Energy Charge (W·h)
    data_columns(124) = {'double', 'Energy_discharge', 8}; % Energy Discharge (W·h)
    data_columns(125) = {'double', 'Capacitance_charge', 8}; % Capacitance Charge (µF)
    data_columns(126) = {'double', 'Capacitance_discharge', 8}; % Capacitance Discharge (µF)
    data_columns(131) = {'uint16', 'Ns', 2}; % Number of Samples
    data_columns(163) = {'single', '|Estack|', 4}; % Absolute Stack Voltage (V)
    data_columns(168) = {'single', 'Rcmp', 4}; % Compensation Resistance (Ω)
    data_columns(169) = {'single', 'Cs', 4}; % Capacitance (µF)
    data_columns(172) = {'single', 'Cp', 4}; % Capacitance Parallel (µF)
    data_columns(173) = {'single', 'Cp⁻²', 4}; % Capacitance Inverse (µF⁻²)
    data_columns(174) = {'single', '<Ewe>', 4}; % Mean Working Electrode Voltage (V)
    data_columns(178) = {'single', '(Q-Qo)', 4}; % Charge Difference (C)
    data_columns(179) = {'single', 'dQ', 4}; % Charge Increment (C)
    data_columns(182) = {'double', 'step time', 8}; % Step Time (s)
    data_columns(185) = {'single', '<Ece>', 4}; % Mean Counter Electrode Voltage (V)
    data_columns(211) = {'double', 'Q_charge_or_discharge', 8}; % Charge or Discharge (C)
    data_columns(217) = {'single', 'THD Ewe', 4}; % Total Harmonic Distortion Ewe (%)
    data_columns(241) = {'single', '|E1|', 4}; % Absolute E1 (V)
    data_columns(242) = {'single', '|E2|', 4}; % Absolute E2 (V)
    data_columns(271) = {'single', 'Phase(Z1)', 4}; % Phase of Z1 (deg)
    data_columns(272) = {'single', 'Analog IN 1', 4}; % Analog Input 1 (V)
    data_columns(295) = {'uint16', 'I Range', 2}; % Current Range
    data_columns(301) = {'single', '|Z1|', 4}; % Impedance of Z1 (Ω)
    data_columns(302) = {'single', '|Z2|', 4}; % Impedance of Z2 (Ω)
    data_columns(326) = {'single', 'P', 4}; % Power (W)
    data_columns(331) = {'single', 'Re(Z1)', 4}; % Real Part of Z1 (Ω)
    data_columns(332) = {'single', 'Re(Z2)', 4}; % Real Part of Z2 (Ω)
    data_columns(361) = {'single', '-Im(Z1)', 4}; % Imaginary Part of Z1 (Ω)
    data_columns(362) = {'single', '-Im(Z2)', 4}; % Imaginary Part of Z2 (Ω)
    data_columns(379) = {'double', 'Energy charge', 8}; % Energy Charge (W·h)
    data_columns(391) = {'single', '<E1>', 4}; % Mean E1 (V)
    data_columns(392) = {'single', '<E2>', 4}; % Mean E2 (V)
    data_columns(422) = {'single', 'Phase(Zstack)', 4}; % Phase of Zstack (deg)
    data_columns(423) = {'single', '|Zstack|', 4}; % Impedance of Zstack (Ω)
    data_columns(424) = {'single', 'Re(Zstack)', 4}; % Real Part of Zstack (Ω)
    data_columns(425) = {'single', '-Im(Zstack)', 4}; % Imaginary Part of Zstack (Ω)
    data_columns(434) = {'single', '(Q-Qo)', 4}; % Charge Difference (C)
    data_columns(435) = {'single', 'dQ', 4}; % Charge Increment (C)
    data_columns(441) = {'single', '<Ece>', 4}; % Mean Counter Electrode Voltage (V)
    data_columns(462) = {'single', 'Temperature', 4}; % Temperature (°C)
    data_columns(467) = {'double', 'Q charge or discharge', 8}; % Charge or Discharge (mA·h)
    data_columns(468) = {'uint32', 'half cycle', 4}; % Half Cycle

    % Define flag bitmask columns (these should be combined into a single Flags byte)
    flag_columns(1)  = {3, 'mode', 2};  
    flag_columns(2)  = {4, 'ox_or_red', 1};  
    flag_columns(3)  = {8, 'error', 1};  
    flag_columns(21) = {16, 'control_changes', 1};  
    flag_columns(31) = {32, 'Ns_changes', 1};  
    flag_columns(65) = {128, 'counter_inc', 1};  

    columnNames = {}; 
    dataTypes = {}; 
    bytesPerCol = [];
    flagInfo = struct();
    haveFlagsColumn = false;

    debug_fprintf("\n=== Debugging lookupColumnTypes ===\n");
    debug_fprintf("Detected columnIDs: %s\n", mat2str(columnIDs));

    for i = 1:numel(columnIDs)
        colID = columnIDs(i);
        
        if ismember(colID, FLAGS_GROUP)
            debug_fprintf("[Flag Detected] ID: %d (Adding to Flags Byte)\n", colID);
            
            if ~haveFlagsColumn
                haveFlagsColumn = true;
                columnNames{end+1} = 'Flags';
                dataTypes{end+1} = 'uint8';
                bytesPerCol(end+1) = 1;
                flagInfo.Flags = struct();
                debug_fprintf("  -> Flags column initialized (1 byte)\n");
            end
    
            flagData = flag_columns(colID);
            bitmask = flagData{1};
            flagName = flagData{2};
            flagInfo.Flags.(flagName) = struct('bitmask', bitmask);
            debug_fprintf("  -> Flag '%s' (bitmask: %d) added to Flags Byte.\n", flagName, bitmask);
    
        elseif isKey(data_columns, colID)
            val = data_columns(colID);
            columnNames{end+1} = val{2};
            dataTypes{end+1} = val{1};
            bytesPerCol(end+1) = val{3};
    
            debug_fprintf("[Column Added] ID: %d -> Name: %s, Type: %s, Bytes: %d\n", ...
                    colID, val{2}, val{1}, val{3});
        else
            debug_fprintf("[Warning] Unknown column ID: %d\n", colID);
            %error('Unknown column ID', colID)
        end
    end

    % Final Debugging Output
    %fprintf("\n=== Final Byte Calculation ===\n");
    %fprintf("Columns Detected: %s\n", strjoin(columnNames, ", "));
    %fprintf("Data Types: %s\n", strjoin(dataTypes, ", "));
    %fprintf("Bytes per Column: %s\n", mat2str(bytesPerCol));
    
    totalBytesPerRow = sum(bytesPerCol);
    %fprintf("Total Bytes Per Row: %d\n", totalBytesPerRow);
end




function out = parseVMPLog(moduleData, modHeader)
%PARSEVMPLOG  Parse the "VMP LOG" module
%   e.g., channel_number at offset 0x0009, etc.

    out = struct();
    % ...
end

function out = parseVMPLoop(moduleData, modHeader)
%PARSEVMPLOOP  Parse the "VMP loop" module
%   e.g., n_indexes at offset 0x0000, etc.

    out = struct();
    % ...
end


function debug_fprintf(varargin)
    % Helper function to print debug logs if flag is true
    % Usage: debug_fprintf(DEBUG_FLAG, 'Your message %d', value);
    DEBUG_TRUE = false; % set to true to enable debug logs!

    if DEBUG_TRUE
        fprintf(varargin{:});
    end
end

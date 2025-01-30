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
        % Convert back to bytes for binary reading
        rawModule = uint8(rawModuleText);

                % Ensure we're reading 'short_name' correctly
        %rawStr = char(rawModule(1:48));  % First 16 bytes should be ASCII
        %disp("Raw extracted short_name:");
        %disp(rawStr);



        % === Parse the module header ===
        % We'll demonstrate reading a few fields as ASCII or numeric at fixed offsets.
        % (Adjust these offsets to match your actual .mpr layout.)
        
        % For demonstration, let's do something like the "module_header_dtypes" approach:
        %   short_name:  at offset 0,  10 bytes
        %   long_name:   at offset 10, 25 bytes
        %   length:      at offset 35, 4 bytes (little-endian uint32)
        %   oldver:      at offset 39, 4 bytes
        %   newver:      at offset 43, 4 bytes
        %   date:        at offset 47, 8 bytes
        %
        % We'll read them carefully, but the real offsets might differ.
        %
        % We also need to ensure we don't exceed rawModule length:
        
        if numel(rawModule) < 55
            % Not enough data to parse all fields
            warning('Module chunk %d is too short to parse header fully.', iModule-1);
            continue;
        end

        shortName = char(rawModule(1:10));
        %disp(shortName)
        longName  = char(rawModule(11:35));
        %disp(longName)

        % lengthBytes = typecast(rawModule(36:39), 'uint32'); % little-endian assumed
        % oldverBytes = typecast(rawModule(40:43), 'uint32');
        % newverBytes = typecast(rawModule(44:47), 'uint32');
        % dateBytes   = char(rawModule(48:55));
        % (5-field version, total 51 bytes in the header)
        lengthBytes = typecast(rawModule(36:39), 'uint32');
        oldverBytes = typecast(rawModule(40:43), 'uint32');
        dateBytes   = char(rawModule(44:51));
        moduleData  = rawModule(52:end);

        lengthVal   = double(lengthBytes);
        %disp("lengthVal: ")

        %disp(lengthVal)

        modHeader = struct();
        modHeader.short_name = strtrim(shortName);
        modHeader.long_name  = strtrim(longName);
        modHeader.length     = double(lengthBytes);
        modHeader.oldver     = double(oldverBytes);
        %modHeader.newver     = double(newverBytes);
        modHeader.date       = strtrim(dateBytes);

        % The actual module data portion is everything after these 55 bytes
        %moduleData = rawModule(56:end);

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
    raw_ns = moduleData(ns_offset : ns_offset + 1);
    ns = typecast(raw_ns, 'uint16');
    swapped_ns = swapbytes(ns);
    ns = min(ns, swapped_ns); % Take the more reasonable value

    % Validate `ns`
    % if ns > 0 && ns < 100
    %     found_ns = true;
    % else 
    %     warning('invalid ns')
    % end



    fprintf('Chosen ns_offset=0x%X: ns=%d\n', ns_offset, ns);

    % `n_params` is always **2 bytes after** `ns`
    np_offset = ns_offset + 2;
    raw_np = moduleData(np_offset : np_offset + 1);
    n_params = typecast(raw_np, 'uint16');
    swapped_np = swapbytes(n_params);
    
    % If swapped version makes more sense, use it
    if swapped_np > 0 && swapped_np < 100
        n_params = swapped_np;
    end

    fprintf('Final chosen n_params offset = 0x%X: n_params=%d\n', np_offset, n_params);
    
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
%PARSEVMPDATA - Parses "VMP data" from Biologic MPR files.

    out = struct();
    out.DataPoints  = [];
    out.Columns     = {};
    out.ColumnIDs   = [];

    offset = 1;
    nDataPoints = typecast(moduleData(offset + (0:3)), 'uint32');
    nColumns    = double(moduleData(offset + 4));

    if nColumns == 0 || nColumns > 100
        warning('[parseVMPData] Failed to detect valid nColumns.');
        return;
    end

    colIDStart = offset + 5;
    colIDEnd = colIDStart + (2 * nColumns) - 1;

    if colIDEnd > length(moduleData)
        warning('[parseVMPData] Not enough bytes for column IDs.');
        return;
    end

    columnIDBytes = moduleData(colIDStart:colIDEnd);
    columnIDs = typecast(columnIDBytes, 'uint16');
    out.ColumnIDs = columnIDs;

    [columnNames, dataTypes, bytesPerCol, flagInfo] = lookupColumnTypes(columnIDs);
    bytesPerRow = sum(bytesPerCol);

    if modHeader.oldver == 3
        dataStart = hex2dec('196');
    elseif modHeader.oldver == 2
        dataStart = hex2dec('195');
    elseif modHeader.oldver == 10 || modHeader.oldver == 11
        dataStart = hex2dec('3EF');
    else
        dataStart = hex2dec('195');
    end

    bytesNeeded = nDataPoints * bytesPerRow;
    dataEnd = dataStart + bytesNeeded - 1;

    if dataEnd > length(moduleData)
        nDataPoints = floor((length(moduleData) - dataStart) / bytesPerRow);
    end

    if nDataPoints <= 0
        return;
    end

    rawData = moduleData(dataStart : dataStart + nDataPoints * bytesPerRow - 1);
    dataMatrix = zeros(nDataPoints, numel(columnNames));
    offset = 1;
    
    for colIdx = 1:numel(columnNames)
        colType = dataTypes{colIdx};
        numBytes = bytesPerCol(colIdx);
        rawColData = rawData(offset:offset + numBytes * nDataPoints - 1);
        
        switch colType
            case 'single'
                colData = typecast(rawColData, 'single');
            case 'uint8'
                colData = typecast(rawColData, 'uint8');
            case 'uint16'
                colData = typecast(rawColData, 'uint16');
            case 'int16'
                colData = typecast(rawColData, 'int16');
            case 'uint32'
                colData = typecast(rawColData, 'uint32');
            case 'int32'
                colData = typecast(rawColData, 'int32');
            case 'double'
                colData = typecast(rawColData, 'double');
            otherwise
                colData = nan(nDataPoints, 1);
        end
        
        dataMatrix(:, colIdx) = colData;
        offset = offset + numBytes * nDataPoints;
    end

    % If there is a flag column, unpack it properly
    if isfield(flagInfo, 'Flags')
        flagIdx = find(strcmp(columnNames, 'Flags'));
        flagData = dataMatrix(:, flagIdx);
        unpackedFlags = zeros(nDataPoints, numel(fieldnames(flagInfo.Flags))); 

        flagNames = fieldnames(flagInfo.Flags);
        for i = 1:numel(flagNames)
            bitmask = flagInfo.Flags.(flagNames{i}).bitmask;
            unpackedFlags(:, i) = bitand(flagData, bitmask) > 0;
        end
        
        % Insert unpacked flag columns in place of "Flags"
        dataMatrix = [dataMatrix(:, 1:flagIdx-1), unpackedFlags, dataMatrix(:, flagIdx+1:end)];
        columnNames = [columnNames(1:flagIdx-1), flagNames', columnNames(flagIdx+1:end)];
    end

    out.DataPoints = dataMatrix;
    out.Columns = columnNames;
end





function [columnNames, dataTypes, bytesPerCol, flagInfo] = lookupColumnTypes(columnIDs)
% LOOKUPCOLUMNTYPES - Maps Biologic column IDs to data type, name, and bitmask details.

    FLAGS_GROUP = [1,2,3,21,31,65];

    data_columns = containers.Map('KeyType', 'double', 'ValueType', 'any');
    flag_columns = containers.Map('KeyType', 'double', 'ValueType', 'any');

    % Define normal columns
    data_columns(4)   = {'double', 'time', 64};  
    data_columns(5)   = {'single', 'control_V_I', 32};  
    data_columns(6)   = {'single', 'Ewe', 32}; 
    data_columns(7)   = {'double', 'dq', 64};
    data_columns(13)  = {'double', '(Q-Qo)', 64};
    data_columns(19)  = {'single', 'control_V', 32};
    data_columns(39)  = {'uint16', 'I_Range', 16};  
    data_columns(467) = {'double', 'Q_charge_or_discharge', 64};

    % Define flag bitmask columns
    flag_columns(1)  = {3, 'mode', 2};  
    flag_columns(2)  = {4, 'ox_or_red', 1};  
    flag_columns(3)  = {8, 'error', 1};  
    flag_columns(21) = {16, 'control_changes', 1};  
    flag_columns(31) = {32, 'Ns_changes', 1};  
    flag_columns(65) = {128, 'counter_inc', 1};  

    columnNames = {}; dataTypes = {}; bytesPerCol = [];
    flagInfo = struct();

    haveFlagsColumn = false;

    for i = 1:numel(columnIDs)
        colID = columnIDs(i);
        if ismember(colID, FLAGS_GROUP)
            if ~haveFlagsColumn
                haveFlagsColumn = true;
                columnNames{end+1} = 'Flags';
                dataTypes{end+1} = 'uint8';
                bytesPerCol(end+1) = 1;  % 1 byte for the entire Flags column
                flagInfo.Flags = struct();
            end
            % Extract the flag column information first
            flagData = flag_columns(colID);
            bitmask = flagData{1};
            flagName = flagData{2};
            numBits = flagData{3};

            % Assign to flagInfo struct safely
            flagInfo.Flags.(flagName) = struct('bitmask', bitmask, 'numBits', numBits);
            
        elseif isKey(data_columns, colID)
            val = data_columns(colID);
            columnNames{end+1} = val{2};
            dataTypes{end+1} = val{1};
            bytesPerCol(end+1) = val{3} / 8;
        end
    end
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

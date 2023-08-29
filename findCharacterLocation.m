function indices = findCharacterLocation(charArray, targetStr)
    % Input:
    % - charArray: The input character array (column vector)
    % - targetStr: The target string to search for (without spaces)


    % Initialize an empty array to store indices
    indices = [];
    
    % Convert target string to a character array and remove spaces
    targetCharArray = targetStr(targetStr ~= ' ');
    %targetCharArray = targetStr;
    % Iterate through the rows of the character array
    for i = 1:size(charArray, 1)
        % Convert the row to a character array and remove spaces
        rowCharArray = charArray(i, :);
        
        % Find the non-space indices
        nonSpaceIndices = find(rowCharArray ~= ' ');
        
        % Remove spaces from the row character array using non-space indices
        rowCharArrayNoSpaces = rowCharArray(nonSpaceIndices);
##        rowCharArrayNoSpaces = rowCharArray;
        % Compare character arrays
        if isequal(rowCharArrayNoSpaces, targetCharArray)
            indices = [indices, i]; % Append the index to the array
        end
    end
    
    % Display the indices if desired
%     if isempty(indices)
%         fprintf('Target string not found in the array.\n');
%     else
%         fprintf('Target string found at row indices: ');
%         fprintf('%d ', indices);
%         fprintf('\n');
%     end
end

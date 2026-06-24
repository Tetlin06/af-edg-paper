function T = paper_append_row(T, row)
    if isempty(T)
        T = row;
    else
        T = [T; row]; %#ok<AGROW>
    end
end

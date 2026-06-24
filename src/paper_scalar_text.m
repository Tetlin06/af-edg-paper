function x = paper_scalar_text(col, r)
%PAPER_SCALAR_TEXT Read one table cell as a trimmed string.
    if iscell(col)
        x = string(col{r});
    else
        x = string(col(r));
    end
    x = strtrim(x);
end

function selected_rows = selectRows(row, n_alt)
fr = floor(row/n_alt);
if fr*n_alt == row 
   rows = row-n_alt+1:row;
else
   rows = fr*n_alt+1:fr*n_alt+n_alt;
end
   selected_rows = rows;
end
function [col,row]=convertGeoCoorMatCoor(rx,ry,x,y,sz)
    % Convert georeferenced coordinates of a point into its matrix
    % coordinate, based on image coordinates
    
    n_rows=sz(1);
    n_cols=sz(2);
    col=interp1([min(rx) max(rx)],[1 n_cols],x);
    row=interp1([max(ry) min(ry)],[1 n_rows],y);
end
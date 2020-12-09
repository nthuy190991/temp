function [x,y]=convertMatCoorGeoCoor(rx,ry,col,row,sz)
    % Convert matrix coordinates of a point into its georeferenced
    % coordinate, based on image coordinates
    
    n_rows=sz(1);
    n_cols=sz(2);
    x=zeros(size(col));
    y=zeros(size(row));
    for i=1:length(col)
        x(i)=interp1([1 n_cols],[min(rx) max(rx)],col(i));
        y(i)=interp1([1 n_rows],[max(ry) min(ry)],row(i));
    end
end
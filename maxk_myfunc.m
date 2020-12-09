function [B,I2]=maxk_myfunc(A,k)
    A2=sort(A,'descend');
    B=A2(1:k);myfunc
    I=zeros(size(B));
    %for i=1:kmyfunc
    i=1;
    while (i<=k) 
        i
        idx=find(B(i)==A);
        if (numel(idx)>1)
            I(i:i+numel(idx)-1)=idx;
            i=i+numel(idx);
        else
            I(i)=idx;
            i=i+1;
        end
        
    end
    I2=I(1:k);
end
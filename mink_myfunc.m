function [B,I2]=mink_myfunc(A,k)
    A2=sort(A);
    B=A2(1:k);
    I=zeros(size(B));
    for i=1:k       
        idx=find(B(i)==A);
        if (numel(idx)>1)
            I(i:i+numel(idx)-1)=idx;
        else
            I(i)=idx;
        end
    end
    I2=I(1:k);
end
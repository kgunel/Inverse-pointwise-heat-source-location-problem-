function y = colorCodeConverterHex2Double(colorArray)
n = length(colorArray);
y =zeros(n,3);
colorArray = colorArray(:,2:end);
for i=1:n

    for j=1:3
        hexCode = colorArray(i,2*j-1:2*j);
        y(i,j) = hex2dec(hexCode)/255;
    end
end

    
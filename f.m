function z = f(source_x, source_y, x, y, epss)
    % Source term function
    % z = -10.*exp(-epss^2./(epss^2-(x-source_x).^2-(y-source_y).^2)); 
    % z = 100.*dirac(x-source_x).*dirac(y-source_y);
    z = -exp(-((x-source_x).^2+(y-source_y).^2)./(4.*epss))./(4.*pi.*epss);
end
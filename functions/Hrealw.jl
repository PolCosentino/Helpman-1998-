function Hrealw(fund, L, w, tradesh)    
    a = fund[:, 1]
    H = fund[:, 2]
    
    #domestic trade share;
    dtradesh = diag(tradesh)
    
    #real wage
    realwage = ((L / (sigma .* F .* dtradesh)).^(alpha / (sigma - 1))) .* (a.^alpha) .* ((L ./ H).^(-(1 - alpha)))
    realwage = realwage / (alpha .* ((sigma / (sigma - 1)).^alpha) .* (((1 - alpha) / alpha).^(1 - alpha)))
    
    return realwage
end

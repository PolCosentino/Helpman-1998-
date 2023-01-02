function Hpindex(fund, L, w, dtradesh)
    
    a = fund[:, 1]
    H = fund[:, 2]
    
    P= (sigma ./ (sigma-1)) .* (w ./a) .* ((L ./ (sigma .* F .* dtradesh)).^(1 ./(1-sigma)))

    
    return P
end

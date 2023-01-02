function Hwelfaregains(ctradesh, tradesh, cL, L)
    
    dtradesh = diag(tradesh)
    cdtradesh = diag(ctradesh)
    
    welfgain = (dtradesh ./ cdtradesh).^(alpha / (sigma - 1)) .* (L ./ cL).^(((sigma .* (1 - alpha)) - 1) / (sigma - 1))
    
    return welfgain
end

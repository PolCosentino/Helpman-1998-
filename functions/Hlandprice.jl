
function Hlandprice(fund, L, w)

    a = fund[:, 1]
    H = fund[:, 2]
    
    r = ((1 - alpha) / alpha) .* ((w .* L) ./ H)
    
    return r
end

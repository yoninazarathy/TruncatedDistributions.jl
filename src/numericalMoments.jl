function LL(d::BoxTruncatedMvNormal)
    @info "doing base numerical integral on dimension $(d.n)."
    hcubature((x)->pdf_nontruncated(d,x),d.a,d.b,maxevals = 10^6)[1]
end
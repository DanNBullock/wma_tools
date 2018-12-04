function fg = bsc_reorientFiber(fg)



for iStreamlines=1:length(fg.fibers)
    fg.fibers{iStreamlines}=bsc_reorientStreamline(fg.fibers{iStreamlines});
end
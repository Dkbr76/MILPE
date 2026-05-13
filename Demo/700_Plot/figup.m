








function popfig()
    figHandles = findall(0, 'Type', 'figure');
    for k = 1:length(figHandles)
        figure(figHandles(k));
    end
end

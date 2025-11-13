function Info_mats = geninfomats(xmat, pts, wts, cset)
    nsamples = length(wts);
Info_mats=zeros(size(xmat,2),size(xmat,2),nsamples);
    for idx = 1:length(wts)
        b = pts(:, idx);
        info = InfoMNL(xmat, b, cset);
        Info_mats(:,:,idx)=info;

        
    end
end
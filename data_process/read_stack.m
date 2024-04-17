function data_stack = read_stack(filename, nt, ntr)

% filename corresponds to a txt (ascii) file

data_from_txt = load(filename);

data_stack = zeros(nt,ntr);

for itr=1:ntr
    
    data_stack(:,itr) = data_from_txt((itr-1)*nt+1:itr*nt,3);

end

% for itr=1:ntr
%     for it=1:nt
%         if isnan(data_Z(it,itr))
%             data_Z(it,itr) = 0.0;
%         end
%     end
% end





end
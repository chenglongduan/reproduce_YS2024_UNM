
% find the adjacent source and receiver spacing
N = 11218;
tmp1 = zeros(1,N);
tmp2 = zeros(1,N);
src_space = zeros(1,N);
rec_space = zeros(1,N);

for i=1:N
    for j=1:N
        tmp1(j) = sqrt((geom_tab(i,1)-geom_tab(j,1))^2+(geom_tab(i,2)-geom_tab(j,2))^2);
        tmp2(j) = sqrt((geom_tab(i,3)-geom_tab(j,3))^2+(geom_tab(i,4)-geom_tab(j,4))^2);
    end
    tmp1(tmp1==0) = NaN;
    tmp2(tmp2==0) = NaN;
    src_space(i) = min(tmp1);
    rec_space(i) = min(tmp2);
end

a = mean(src_space);
b = mean(rec_space);
c = max(src_space);
d = max(rec_space);

fprintf('mean source spacing = %f\n',a);
fprintf('mean receiver spacing = %f\n',b);
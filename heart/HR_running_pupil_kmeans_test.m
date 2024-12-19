


thresh = 0.02;

pre_bout_thresh = 4;

window = 10;

hz = 1000;


binarized = x > thresh;

run_start = find(diff(binarized)==1);

good_bouts = [];

for i = 1:length(run_start)

    ii = sum(binarized(run_start(i)-pre_bout_thresh*1000:run_start(i)));

    if ii == 0
        good_bouts = [good_bouts,run_start(i)];
    end
end

X = zeros(length(good_bouts),(2*window)*hz+1);


for i = 1:length(good_bouts)

    ii = good_bouts(i);

    X(i,:) = x(ii-window*hz:ii+window*hz);

end


X_input = X(:,window*hz-4*hz:end);
X_input = zscore(X_input,0);

km = kmeans(X_input,2);



plot(X(km==1,:)')
figure;
plot(X(km==2,:)')

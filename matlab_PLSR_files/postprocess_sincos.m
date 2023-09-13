function out = postprocess_sincos(predicted)

%= Change Cartesian coordinates back to angles 

magnitude = sqrt(predicted(:,1) .^ 2 + predicted(:,2) .^2);

out = zeros(size(magnitude));
% normalize the magnitude
predicted = predicted ./ repmat(magnitude, 1, 2);

angle = atand(predicted(:,2) ./ predicted(:,1));

%process quadrant 1
selected = find(predicted(:,1) > 0 &  predicted(:,2) >0);
out(selected) = angle(selected);

%process quadrant 2
selected = find(predicted(:,1) <= 0 &  predicted(:,2) >0);
out(selected) = 180 + angle(selected);

%process quadrant 3
selected = find(predicted(:,1) <= 0 &  predicted(:,2) <=0);
out(selected) = 180 + angle(selected);

%process quadrant 4
selected = find(predicted(:,1) > 0 &  predicted(:,2) <=0);
out(selected) = 360 + angle(selected);

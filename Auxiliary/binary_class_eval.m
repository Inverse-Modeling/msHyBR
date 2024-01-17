function  [TPR,TNR,PPV, F1, TP, TN, FP, FN] = binary_class_eval(beta_ind, beta_true)


P = sum(beta_true == 1);
N = sum(beta_true == 0);

test = beta_ind+beta_true;
% sum component - wise
% count number of 2 == true positive
% count number of 0 == true negative
TP = sum(test == 2); 
TN = sum(test == 0);

test2 = beta_ind-beta_true;
% subtract component - wise
% count number of 2 == true positive
% count number of 0 == true negative
FP = sum(test2 == 1);
FN = sum(test2 == -1);



TPR = sum(test == 2) / P;
TNR = sum(test == 0) / N;
PPV = sum(test == 2)/(sum(test == 2)+sum(test2 == 1));
F1 = 2*(PPV*TPR)/(PPV+TPR);

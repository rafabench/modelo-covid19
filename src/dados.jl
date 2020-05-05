diagn_Rio  = [26, 38, 30, 114, 110, 60, 115, 86, 42, 141, 198, 184, 175, 97, 91, 326, 71, 126, 140, 287, 113, 68, 116]
obitos_Rio = [1,   3,  2,   4,   7,  5,   2,  4,  5,  12,  14,   9,  10,  6,  8,   9, 25,  27,  28,  24,  18,  8,  12]
Rio_acc_d = cumsum([498.; diagn_Rio])
Rio_acc_m = cumsum([15.; obitos_Rio]);

n_pts = length(diagn_Rio)
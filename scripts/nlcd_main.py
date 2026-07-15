import os
os.environ['OMP_NUM_THREADS'] = '1'  ## This will prevent KRR algo from thrashing the server
from scipy.stats import norm
import random
import numpy as np
from sklearn.svm import SVR
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import mean_squared_error
from sklearn.neural_network import MLPRegressor
from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning
@ignore_warnings(category=ConvergenceWarning)
def normalize(x_test,y_test):
    '''
    Function to standard normalize A|L for all the values of L to mean of 0 and variance 1
    '''
    unique_values=np.unique(x_test)
    for value in unique_values:
        indices = np.where(x_test == value)[0]
        sigma = np.std(y_test[indices])
        y_test[indices] = (y_test[indices])/sigma
    return y_test

_DEFAULT_PARAMS = {
    "SVR": {"kernel": "rbf"},
    "KRR": {"kernel": "rbf"},
    "ANN": {},
}

_REGRESSORS = {
    "SVR": SVR,
    "KRR": KernelRidge,
    "ANN": MLPRegressor,
}

def nlr_train_predict(xG, yG, algo, xL=None):
    '''
    Non-linear regression (NLR);
    `algo` may be either:
    - a string: "SVR" | "KRR" | "ANN"  (uses original defaults)
    - a dict: {"name": "KRR", "params": {<any sklearn kwarg>: value, ...}}
        e.g. {"name": "KRR", "params": {"kernel": "rbf", "alpha": 0.1, "gamma": 0.5}}
            {"name": "SVR", "params": {"kernel": "poly", "C": 10.0, "degree": 3}}
            {"name": "ANN", "params": {"hidden_layer_sizes": (64, 32), "max_iter": 500}}
    '''
    if xL is None:
        x = xG.reshape(-1, 1)
    else:
        x = np.column_stack((xL, xG))

    if isinstance(algo, dict):
        name = algo["name"]
        user_params = algo.get("params", {})
    else:
        name = algo
        user_params = {}

    if name not in _REGRESSORS:
        print("Invalid Algorithm")
        assert False

    # Merge: start from original defaults, let user_params override
    params = {**_DEFAULT_PARAMS[name], **user_params}
    regressor = _REGRESSORS[name](**params)

    regressor.fit(x, yG)
    y_pred = regressor.predict(x)
    return y_pred, regressor

def calculate_pvalue(original, loss_list):
    '''
    calculate the p value (with +1/+2 correction)
    '''
    pvalNr = sum(i <= original for i in loss_list)
    return (pvalNr + 1)/(len(loss_list) + 2)


def get_prob(L,A):
    '''
    Function to calculate the overlap score
    '''
    norm_const=0
    probs=[]
    unique_values = np.unique(L)
    for value in unique_values:
        indices = np.where(L == value)[0]
        mu = np.mean(A[indices])
        sigma = np.std(A[indices])
        p_L = np.count_nonzero(L == value) / len(L)
        p_AgivenL = norm.pdf(A, mu, sigma)
        p_LgivenA = p_AgivenL * p_L
        probs.append(p_LgivenA)
        norm_const += p_LgivenA
    diff=[]
    for i in range(len(probs)):
        for j in range(i+1,len(probs)):
            diff.append(np.minimum(probs[i],probs[j])/norm_const)

    return diff


def FI_score(x,y,overlap ):
    '''
    Function to calcuate the overlap score
    '''

    # Iterate over the pairs of lists and calculate element-wise absolute differences
    diff = np.abs(np.subtract(x, y))

    return np.sum(np.multiply(diff,overlap))/len(x)


def compute_Luniqs_predns(L,A,B,algo):
    '''
    Function to calculate the predictions of y for different values of L in the 4th (or 2nd) test
    '''
    y_pred = []
    _, regressor = nlr_train_predict(A, B, algo, L)


    unique_values = np.unique(L)
    assert (len(unique_values)==2 and (unique_values==[0,1]).all()) or (len(unique_values)==3 and (unique_values==[0,1,2]).all())
    for value in unique_values:
        L_1 = np.full_like(L, value)
        X_ = np.column_stack((L_1, A))
        y_predict_ = regressor.predict(X_)
        y_pred.append(y_predict_)

    return y_pred


def stratify_permute_variable(L, variable):
    '''
    Function to stratify permute a variable
    '''
    unique_values = np.unique(L)
    permuted_variable = np.empty_like(variable)

    for value in unique_values:
        indices = np.where(L == value)[0]
        permuted_indices = np.random.choice(indices, len(indices),replace=False)
        permuted_variable[indices] = variable[permuted_indices]

    return permuted_variable


def test_4(L,A,B,shuffles,algo):
    '''
    Function for the 4th test , same function is used for the second test
    '''

    overlap=get_prob(L,A)
    perm_loss=[]
    for i in range(shuffles):
            y_pred=compute_Luniqs_predns(L,A,stratify_permute_variable(L,B),algo)

            total_FI=0
            count=0
            for i in range(len(y_pred)):
                for j in range(i + 1, len(y_pred)):
                    total_FI+=FI_score(y_pred[i],y_pred[j],overlap[count])
                    count+=1
            assert count>0
            total_FI/=count
            perm_loss.append(total_FI)

    y_pred_original=compute_Luniqs_predns(L,A,B,algo)
    assert len(y_pred)==len(y_pred_original)
    assert len(y_pred)==2 or len(y_pred)==3
    original_loss=0
    count=0
    for i in range(len(y_pred_original)):
        for j in range(i + 1, len(y_pred_original)):
            original_loss+=FI_score(y_pred_original[i],y_pred_original[j],overlap[count])
            count+=1
    assert count>0
    original_loss/=count
    overlap_scores = [np.sum(overlaps) / len(overlaps) for overlaps in overlap]
    total_overlap_score = np.sum(overlap_scores)/count
    return [calculate_pvalue(original_loss,perm_loss),total_overlap_score,original_loss,perm_loss]


def compute_1_loss(x_test,y_test):
    '''
    Function to calculate the loss in the first test using nll
    '''
    unique_values = np.unique(x_test)

    y_mean = np.empty_like(y_test)
    y_std = np.empty_like(y_test)

    for value in unique_values:
        indices = np.where(x_test == value)[0]
        mu, sigma = np.mean(y_test[indices]), np.std(y_test[indices])
        y_mean[indices] = mu
        y_std[indices] = sigma
    #return NLL loss
    #return np.mean(((((y_test-y_mean)/y_std)**2)/2)+np.log(2*np.pi)+np.log((y_std)))
    return -np.mean(np.log(np.exp(-(((y_test-y_mean)/y_std)**2)/2)/np.sqrt(2*np.pi*(y_std**2))))


def test_1(L,B,shuffles):
    '''
    Function for the first test
    '''
    perm_loss=[]
    for i in range(0,shuffles):
        perm_loss.append(compute_1_loss(L,np.random.permutation(B)))

    original_loss=compute_1_loss(L,B)
    return calculate_pvalue(original_loss,perm_loss),original_loss,perm_loss


def test_3_loss(A,B,algo):
    '''
    function to return the loss in the third function
    '''
    Bpred, _ = nlr_train_predict(A, B, algo)
    mse = mean_squared_error(B, Bpred)
    return mse


def test_3(L,A,B,shuffles,algo):
    '''
    function for the third test
    '''
    unique_values = np.unique(L)
    p_values = []
    original_losses=[]
    perm_losslist=[]
    for value in unique_values:
        indices = np.where(L == value)[0]
        A_value = A[indices]
        B_value = B[indices]

        original_loss = test_3_loss(A_value, B_value, algo)
        original_losses.append(original_loss)
        perm_losses = [test_3_loss(A_value, np.random.permutation(B_value), algo) for _ in range(shuffles)]
        perm_losslist.append(perm_losses)
        p_value = calculate_pvalue(original_loss, perm_losses)
        p_values.append(p_value)

    #min_original_loss=original_losses[p_values.index(min(p_values))]
    min_p_value = min(p_values)
    return min_p_value,original_losses,perm_losslist


def test_2(L,A,B,shuffles,algo):
    Apred, _ = nlr_train_predict(B, A, algo)
    Aresid = A - Apred
    out,original_loss,perm_loss = test_1(L, Aresid, shuffles)
    return out,original_loss,perm_loss


def combine_tests(L,A,B,shuffles,algo,normal):
    '''
    Function to combine all the tests
    '''
    if( len(np.unique(L))==1): # fix the number 10
        print("Zero variance of L, skipping trio")
        return [None]
    if( 1 in [sum(L==x) for x in np.unique(L)] ):
        print(" Only single value for a genotype value, can't do the statistics ")
        return [None]
    if(normal==True):
        A=normalize(L,A)
        B=normalize(L,B)

    T1_p,T1_OL,T1_PL=test_1(L,B,shuffles)
    T2_p,T2_OL,T2_PL=test_2(L,A,B,shuffles,algo)
    T3_p,T3_OL,T3_PL=test_3(L,A,B,shuffles,algo)
    T4_p,T4_OS,T4_OL,T4_PL=test_4(L,A,B,shuffles,algo)
    p_final=np.max([T1_p,T2_p,T3_p,T4_p])
    return [p_final,T1_p,T2_p,T3_p,T4_p,T4_OS,T1_OL,T2_OL,T3_OL,T4_OL,T1_PL,T2_PL,T3_PL,T4_PL]


# ======================================================================
# Train/test split variants (added for R2C4: addresses overfitting concern
# by fitting models/MLEs on a train split and evaluating loss on a held-out
# test split, for both original and permuted data. A fresh random split is
# drawn for the original run and for each permuted replicate.
# ======================================================================

def _split_indices(n, train_frac):
    '''
    Plain random train/test split. Uses np.random (already seeded upstream
    in nlcd_single) so split draws are reproducible per trio seed.
    '''
    idx = np.random.permutation(n)
    n_train = int(np.floor(n * train_frac))
    if n_train < 1 or n_train >= n:
        raise ValueError(f"Bad train_frac={train_frac} for n={n}")
    return idx[:n_train], idx[n_train:]


def _fit_per_L_gaussian(L_train, y_train):
    '''
    Per-L MLE of (mu, sigma) on train data. Returns dict {L_value: (mu, sigma)}.
    Strata with <2 train samples get sigma=nan and are excluded from scoring.
    '''
    params = {}
    for value in np.unique(L_train):
        idx = np.where(L_train == value)[0]
        if len(idx) < 2:
            params[value] = (np.nan, np.nan)
        else:
            params[value] = (np.mean(y_train[idx]), np.std(y_train[idx]))
    return params


def _nll_with_train_params(L_test, y_test, train_params):
    '''
    NLL of test samples scored under per-L (mu, sigma) fit on train.
    Test samples whose L was not seen (or had <2 train samples) are dropped.
    Returns np.inf if nothing remains scorable.
    '''
    mu_arr = np.empty(len(y_test), dtype=float)
    sigma_arr = np.empty(len(y_test), dtype=float)
    mask = np.zeros(len(y_test), dtype=bool)
    for value, (mu, sigma) in train_params.items():
        if np.isnan(sigma) or sigma == 0:
            continue
        idx = np.where(L_test == value)[0]
        mu_arr[idx] = mu
        sigma_arr[idx] = sigma
        mask[idx] = True
    if not mask.any():
        return np.inf
    y = y_test[mask]
    mu = mu_arr[mask]
    sigma = sigma_arr[mask]
    return -np.mean(np.log(np.exp(-(((y - mu) / sigma) ** 2) / 2) / np.sqrt(2 * np.pi * (sigma ** 2))))


def test_1_split(L, B, shuffles, train_frac=0.8):
    '''
    Test 1 with train/test split:
      - fit per-L (mu, sigma) on train (L, B)
      - score NLL on held-out test
      - permute B then redraw a fresh split for each replicate
    '''
    n = len(L)

    def _one_run(L_eff, B_eff):
        train_idx, test_idx = _split_indices(n, train_frac)
        params = _fit_per_L_gaussian(L_eff[train_idx], B_eff[train_idx])
        return _nll_with_train_params(L_eff[test_idx], B_eff[test_idx], params)

    original_loss = _one_run(L, B)
    perm_loss = [_one_run(L, np.random.permutation(B)) for _ in range(shuffles)]
    return calculate_pvalue(original_loss, perm_loss), original_loss, perm_loss


def test_2_split(L, A, B, shuffles, algo, train_frac=0.8):
    '''
    Test 2 with train/test split:
      - on train: fit NLR B->A; compute train A residuals; fit per-L (mu,sigma)
        on (L_train, A_resid_train).
      - on test: predict A using the same trained regressor; compute test A
        residuals; score NLL on (L_test, A_resid_test) under train MLE params.
      - permute L (full) then redraw a fresh split for each replicate.
    '''
    n = len(L)

    def _one_run(L_eff):
        train_idx, test_idx = _split_indices(n, train_frac)
        A_tr, B_tr = A[train_idx], B[train_idx]
        A_te, B_te = A[test_idx], B[test_idx]
        L_tr = L_eff[train_idx]
        L_te = L_eff[test_idx]

        # fit B -> A on train (no L)
        _, regressor = nlr_train_predict(B_tr, A_tr, algo)

        A_tr_pred = regressor.predict(B_tr.reshape(-1, 1))
        A_tr_resid = A_tr - A_tr_pred
        A_te_pred = regressor.predict(B_te.reshape(-1, 1))
        A_te_resid = A_te - A_te_pred

        params = _fit_per_L_gaussian(L_tr, A_tr_resid)
        return _nll_with_train_params(L_te, A_te_resid, params)

    original_loss = _one_run(L)
    perm_loss = [_one_run(np.random.permutation(L)) for _ in range(shuffles)]
    return calculate_pvalue(original_loss, perm_loss), original_loss, perm_loss


def _test_3_split_loss(A, B, algo, train_frac):
    '''
    Single stratum: split (A, B); fit NLR A->B on train; return test MSE.
    '''
    n = len(A)
    train_idx, test_idx = _split_indices(n, train_frac)
    A_tr, B_tr = A[train_idx], B[train_idx]
    A_te, B_te = A[test_idx], B[test_idx]
    _, regressor = nlr_train_predict(A_tr, B_tr, algo)
    B_te_pred = regressor.predict(A_te.reshape(-1, 1))
    return mean_squared_error(B_te, B_te_pred)


def test_3_split(L, A, B, shuffles, algo, train_frac=0.8):
    '''
    Test 3 with train/test split, per L stratum:
      - within each L stratum: split (A, B); fit NLR A->B on train; MSE on test.
      - permute B within stratum; redraw fresh split per replicate.
      - p-value per stratum; take min across strata.
    '''
    unique_values = np.unique(L)
    p_values = []
    original_losses = []
    perm_losslist = []
    for value in unique_values:
        indices = np.where(L == value)[0]
        A_value = A[indices]
        B_value = B[indices]

        original_loss = _test_3_split_loss(A_value, B_value, algo, train_frac)
        original_losses.append(original_loss)
        perm_losses = [_test_3_split_loss(A_value, np.random.permutation(B_value), algo, train_frac)
                       for _ in range(shuffles)]
        perm_losslist.append(perm_losses)
        p_values.append(calculate_pvalue(original_loss, perm_losses))

    return min(p_values), original_losses, perm_losslist


def _test_4_split_score(L_full_unique, L_eff, A_eff, B_eff, algo, train_frac):
    '''
    Single run of test 4 with split:
      - split (L, A, B); fit NLR (L_tr, A_tr) -> B_tr.
      - compute overlap on test (L_te, A_te).
      - for each L value in L_full_unique, set L_te to that value, predict B
        on test, accumulate pairwise FI scores weighted by overlap.
    Returns (avg_FI, avg_overlap_score). If train missing any L value,
    returns (np.inf, 0.0) so that run signals an unscorable replicate.
    '''
    n = len(L_eff)
    train_idx, test_idx = _split_indices(n, train_frac)
    L_tr, A_tr, B_tr = L_eff[train_idx], A_eff[train_idx], B_eff[train_idx]
    L_te, A_te = L_eff[test_idx], A_eff[test_idx]

    # need every L value present in train so the (L,A)->B model can be queried
    # at every L level
    if not np.all(np.isin(L_full_unique, np.unique(L_tr))):
        return np.inf, 0.0

    overlap = get_prob(L_te, A_te)
    _, regressor = nlr_train_predict(A_tr, B_tr, algo, L_tr)

    y_pred = []
    for value in L_full_unique:
        L_te_set = np.full_like(L_te, value)
        X_te = np.column_stack((L_te_set, A_te))
        y_pred.append(regressor.predict(X_te))

    total_FI = 0
    count = 0
    for i in range(len(y_pred)):
        for j in range(i + 1, len(y_pred)):
            total_FI += FI_score(y_pred[i], y_pred[j], overlap[count])
            count += 1
    assert count > 0
    total_FI /= count
    overlap_scores = [np.sum(o) / len(o) for o in overlap]
    total_overlap_score = np.sum(overlap_scores) / count
    return total_FI, total_overlap_score


def test_4_split(L, A, B, shuffles, algo, train_frac=0.8):
    '''
    Test 4 with train/test split:
      - on train: fit NLR (L, A) -> B.
      - on test: for each L level, set L_test to that level, predict B,
        compute pairwise FI weighted by test-set overlap.
      - permuted via stratify_permute_variable(L, B); fresh split per replicate.
    '''
    L_full_unique = np.unique(L)

    original_loss, total_overlap_score = _test_4_split_score(L_full_unique, L, A, B, algo, train_frac)
    perm_loss = []
    for _ in range(shuffles):
        B_perm = stratify_permute_variable(L, B)
        loss_i, _ = _test_4_split_score(L_full_unique, L, A, B_perm, algo, train_frac)
        perm_loss.append(loss_i)

    return [calculate_pvalue(original_loss, perm_loss), total_overlap_score, original_loss, perm_loss]


def combine_tests_split(L, A, B, shuffles, algo, train_frac=0.8):
    '''
    Train/test-split variant of combine_tests. Same return layout as
    combine_tests so nlcd_batch packaging is unchanged. Note: --normal
    preprocessing is intentionally not supported here (the global per-L
    normalization would leak test-set statistics into preprocessing).
    '''
    if len(np.unique(L)) == 1:
        print("Zero variance of L, skipping trio")
        return [None]
    if 1 in [sum(L == x) for x in np.unique(L)]:
        print(" Only single value for a genotype value, can't do the statistics ")
        return [None]
    # ensure each stratum has enough samples that a random 80/20 split
    # is unlikely to empty out either side
    if min(sum(L == x) for x in np.unique(L)) < 5:
        print(" Insufficient samples per L stratum for train/test split (need >=5 per stratum), skipping trio")
        return [None]

    T1_p, T1_OL, T1_PL = test_1_split(L, B, shuffles, train_frac)
    T2_p, T2_OL, T2_PL = test_2_split(L, A, B, shuffles, algo, train_frac)
    T3_p, T3_OL, T3_PL = test_3_split(L, A, B, shuffles, algo, train_frac)
    T4_p, T4_OS, T4_OL, T4_PL = test_4_split(L, A, B, shuffles, algo, train_frac)
    p_final = np.max([T1_p, T2_p, T3_p, T4_p])
    return [p_final, T1_p, T2_p, T3_p, T4_p, T4_OS, T1_OL, T2_OL, T3_OL, T4_OL, T1_PL, T2_PL, T3_PL, T4_PL]

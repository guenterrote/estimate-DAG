from itertools import permutations
from random import choice, choices, random, randint, seed #,sample
from collections import defaultdict
from math import sqrt, factorial
import sys
from estimating_DAG_klee_minty import run_method
from estimating_DAG_klee_minty import random_selection_from_class
from estimating_DAG_klee_minty import commit_to_file

n = 8
n=5
if len(sys.argv)>1:
    n = int(sys.argv[1])

num_samples = 10**7
crowd_budget_B = 10 # budget for CS
problem = 'PERM'
print(f"PERMUTATIONS == {n = }, {num_samples} samples")
true_value = factorial(n)
print(f"EXACT ANSWER: {n}! = {true_value}")

fname = f"results-{problem}-{n}.py"
resultfile = open(fname,"w")

source_vertex = tuple(reversed(range(n)))

def outnbrs(p):
    return list(i for i,(p,q) in enumerate(zip(p,p[1:])) if p>q)
def innbrs(p):
    return list(i for i,(p,q) in enumerate(zip(p,p[1:])) if p<q)
def swap(p,i):
    "works equally for outneighbor or in_neighbor"
    p[i:i+2] = [p[i+1],p[i]]
def indegree(p):    
    return sum(1 for (p,q) in (zip(p,p[1:])) if p<q)

mean = dict()
Var = dict()

def compute_variance():
    """compute mean and variance of the estimate,
    proceeding from the "leaves"/target towards the source vertex"""
    target = tuple(range(n))
    old_layer = [target]
    mean[target] = 1
    Var[target] = 0
    level = 0
    while 1:
        print(f"level {level}, size {len(old_layer)}")
        level += 1
        new_layer = set()
        for p in old_layer:
            p = list(p)
            for i in innbrs(p):
                swap(p,i)
                new_layer.add(tuple(p))
                swap(p,i) # undo the change
        if not new_layer:
            break
        for p0 in new_layer:
            p = list(p0)
            mu = 0
            out_n = outnbrs(p)
            d = len(out_n)
            for i in out_n:
                swap(p,i)
                mu += mean[tuple(p)]/indegree(p)
                swap(p,i) # undo the change
            mu_d = mu/d
            V = 0
            for i in out_n:
                swap(p,i)
                q_in = indegree(p)
                V += Var[tuple(p)]/q_in**2 + (mean[tuple(p)]/q_in - mu_d)**2
                swap(p,i) # undo the change
            mean[p0] = 1+mu
            Var[p0] = d*V
        old_layer = new_layer
    return mean[source_vertex], Var[source_vertex]

Me,Va = compute_variance()
print(f"P mean={Me}, variance={Va}={Va:g}, S.D.={sqrt(Va)}={sqrt(Va):g}")
commit_to_file(Me,Va, 'Algorithm P:', true_value, resultfile, problem)

def compute_variance_importance_sampling():
    """compute mean and variance of the estimate,
    proceeding from the "leaves"/target towards the source vertex"""
    target = tuple(range(n))
    old_layer = [target]
    mean[target] = 1
    Var[target] = 0
    level = 0
    while 1:
        #print(f"level {level}, size {len(old_layer)}")
        level += 1
        new_layer = set()
        for p in old_layer:
            p = list(p)
            for i in innbrs(p):
                swap(p,i)
                new_layer.add(tuple(p))
                swap(p,i) # undo the change
        if not new_layer:
            break
        for p0 in new_layer:
            p = list(p0)
            successors = outnbrs(p)
            d = len(successors)
            mu = 0
            prob_out = []
            for i in successors:
                # fast incremental computation of new outdegree without actually carrying out the swap
                outdeg_new = (d-1 +
                          (1 if i>0   and p[i+1]<p[i-1]<p[i] else 0) +
                          (1 if i<n-2 and p[i+1]<p[i+2]<p[i] else 0))
                prob_out.append(outdeg_new)              
                swap(p,i)
                mu += mean[tuple(p)]/(n-1-outdeg_new)
                swap(p,i) # undo the change
            tot_prob = sum(prob_out)
            V = 0
            if tot_prob == 0:
                assert d==1 # last level
                i = successors[0]
                swap(p,i)
                assert Var[tuple(p)]==0
            else:
                mu_t = mu/tot_prob
                for i,prob_i in zip(successors,prob_out):
                    outdeg_new = prob_i # outdeg is the same as prob
                    q_in = n-1-outdeg_new
                    swap(p,i)
                    V += ( Var[tuple(p)]/q_in**2 + (mean[tuple(p)]/q_in - mu_t*prob_i)**2 ) / prob_i
                    swap(p,i) # undo the change
            mean[p0] = 1+mu
            Var[p0] = tot_prob*V
        old_layer = new_layer
    return mean[source_vertex], Var[source_vertex]

Me_IS,Va_IS = compute_variance_importance_sampling()
print(f"IS mean={Me_IS}, variance={Va_IS}={Va_IS:g}, S.D.={sqrt(Va_IS)}={sqrt(Va_IS):g}")
commit_to_file(Me_IS,Va_IS, "Importance sampling (IS) by outdegree:", true_value, resultfile, problem)


def insert(p,W,stratum,Store):
    z = Store[stratum]
    if z is None:
        W_total = W
    else:
        p_old,W_old = z
        W_total = W_old+W
        if random() <= W_old/W_total:
            p = p_old
    Store[stratum] = (p,W_total)

def bucket_sampling(report_levels=False):
    """heuristic sampling according to Pang Chen.
    stratum = outdegree
    """
    Store = [None]*(n+1)
    X = 0
    num_visited = 0
    num_inserted = 0
    Store[n] = source_vertex,1 # start with W=1
    active = True
    while active:
        for outdeg in range(n,-1,-1):
            z = Store[outdeg] # expand highest outdegree
            #print(level,outdeg,z,Store)
            if z is not None:
                Store[outdeg] = None
                num_visited += 1
                # expand children
                p0,W = z                
                X += W
                q = list(p0)
                successors = outnbrs(q)
                d = len(successors)
                num_inserted += d
                for i in successors:
                    swap(q,i)
                    # fast incremental computation of new outdegree
                    outdeg_new = (d-1 +
                                  (1 if i>0   and q[i]<q[i-1]<q[i+1] else 0) +
                                  (1 if i<n-2 and q[i]<q[i+2]<q[i+1] else 0))
                    assert outdeg_new == n-1-indegree(q)
                    indeg_new = n-1-outdeg_new
                    insert(tuple(q), W/indeg_new, outdeg_new, Store)
                    swap(q,i) # undo the change               
                    insert(q, W/indeg_new, outdeg_new, Store)
                break
        else:
            active = False
    return X,num_visited,num_inserted                

def path_sampling():
    "Algorithm P"
    num_visited = 0
    X = 0
    u = list(source_vertex) # u will be changed during the walk
    W = 1
    while 1:
        num_visited += 1
        X += W
        successors = outnbrs(u)
        d_out = len(successors)
        #print(u,successors,indegree(u),f"{W=}, {X=}")
        if d_out==0:
            return X,num_visited,0
        W *= d_out
        i = choice(successors)
        swap(u,i)
        W /= indegree(u)

def crowd_sampling():
    "Algorithm CS" # crowd sampling
    num_visited = 0
    num_inserted = 0
    num_condensed = 0 # after eliminating multiple equal elements
    X = 0
    A = [(source_vertex,1)]
    while A:
        new_A = defaultdict(float)
        for u,W in A:
            u = list(u)
            num_visited += 1
            X += W
            #successors = outnbrs(u)
            # d = len(successors)
            for i in outnbrs(u):
                num_inserted += 1
                swap(u,i)
                # indegree computation could be speeded up
                new_A[tuple(u)] += W/indegree(u)
                swap(u,i) # undo
        num_condensed += len(new_A)
        new_A = list(new_A.items())
        if len(new_A) <= crowd_budget_B:
            A = new_A
        else:
            A = [random_selection_from_class(new_A,k,crowd_budget_B)
                 for k in range(crowd_budget_B)]
    return X,num_visited,num_inserted, num_inserted-num_condensed
        

def importance_sampling():
    num_visited = 0
    X = 0
    u = list(source_vertex) # u will be changed during the walk
    W = 1
    while 1:
        num_visited += 1
        X += W
        succ = outnbrs(u)
        d = len(succ)
        if d==0:
            return X,num_visited,0
        prob_out = []
        for i in succ:
            # fast incremental computation of new outdegree without actually carrying out the swap
            outdeg_new = (d-1 +
                          (1 if i>0   and u[i+1]<u[i-1]<u[i] else 0) +
                          (1 if i<n-2 and u[i+1]<u[i+2]<u[i] else 0))
            prob_out.append(outdeg_new)
        tot_prob = sum(prob_out)
        if tot_prob == 0:
            assert d==1
            i = succ[0]
            swap(u,i)
            W *= 1/indegree(u)
        else:
            ind = choices(range(d), prob_out)[0] # weighted selection.
            i = succ[ind]
            swap(u,i)
            W *= tot_prob/(prob_out[ind]*indegree(u))

print("first run P: ", path_sampling())
print("first run BS:", bucket_sampling())
print("first run IS:", importance_sampling())


powers_of_ten = [10**i for i in range(1,20)]

for name,method,rand_seed,num_s in [
        ("Algorithm P:", path_sampling, 327069990,num_samples),
        ("Importance sampling (IS) by outdegree:", importance_sampling, None,num_samples),
        (f"Crowd sampling (CS) ({crowd_budget_B=}):", crowd_sampling, None,
           num_samples//10),
        ("Bucket sampling (BS) by outdegree:", bucket_sampling, None,
           num_samples//10), # bucket sampling does more work per sample.
    ]:
    run_method(resultfile, problem,n,name,method,rand_seed,num_s, true_value)
resultfile.close()

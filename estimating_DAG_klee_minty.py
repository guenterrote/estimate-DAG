from itertools import permutations
from random import choice, choices, random, randint, seed #,sample
from collections import defaultdict
from math import sqrt
import sys

n = 12
n=6

if len(sys.argv)>1:
    n = int(sys.argv[1])

num_samples = 10**7
horde_limit = 10 # budget for SE        
problem = 'KM'


source_vertex = 0

def outnbrs(p):
    return list(i for i in range(n) if p&(1<<i)==0)
def innbrs(p):
    return list(i for i in range(n) if p&(1<<i))
def swap(p,i):
    "works equally for outneighbor or in_neighbor"
    return p ^ ((1<<(i+1))-1)
def indegree(p):    
    return sum(1 for i in range(n) if p&(1<<i))
def outdegree(p):    
    return sum(1 for i in range(n) if p&(1<<i)==0)

if 0: # explore various methods for assigning vertex labels
    # check out-degrees
    n=3
    for p in range(2**n):
        o = outnbrs(p)
        print(p, [swap(p,i) for i in o], "       ", o,)
    exit()
              
mean = dict()
Var = dict()

def compute_variance():
    """compute mean and variance of the estimate,
    proceeding from the "leaves"/target towards the source vertex"""
#    target = 2**n-1
 #   mean[target] = 1
  #  Var[target] = 0
    for p in range(2**n-1,-1,-1):
        mu = 0
        out_n = outnbrs(p)
        d = len(out_n)
        for i in out_n:
            q = swap(p,i)
            mu += mean[q]/indegree(q)
        if d:
            mu_d = mu/d
        V = 0
        for i in out_n:
            q = swap(p,i)
            q_in = indegree(q)
            V += Var[q]/q_in**2 + (mean[q]/q_in - mu_d)**2
        mean[p] = 1+mu
        Var[p] = d*V
    return mean[source_vertex], Var[source_vertex]

def compute_variance_importance_sampling(use_size=0):
    """compute mean and variance of the estimate,
    proceeding from the "leaves"/target towards the source vertex"""
    target = 2**n-1
    mean[target] = 1
    Var[target] = 0
    for u in range(2**n-2,-1,-1):
        successors = outnbrs(u)
        d = len(successors)
        assert d>0
        prob_out = []
        mu = 0
        for i in successors:
            v = swap(u,i)
            outdeg_new = outdegree(v)
            if use_size==1:
                weight = 2**n-v # size of reachable subgraph
            elif use_size==2:
                weight = (2**n-v)/indegree(v)
            else:
                weight = max(outdeg_new,1)
                # at least 1: makes variance slightly worse
            prob_out.append(weight) 
            mu += mean[v]/(n-outdeg_new)
        tot_prob = sum(prob_out)
        V = 0
        mu_t = mu/tot_prob
        for i,prob_i in zip(successors,prob_out):
            if prob_i>0:
                v = swap(u,i)
                q_in = indegree(v)
                V += ( Var[v]/q_in**2 + (mean[v]/q_in - mu_t*prob_i)**2 ) / prob_i
        mean[u] = 1+mu
        Var[u] = tot_prob*V
    return mean[source_vertex], Var[source_vertex]

def commit_to_file(mean,var, name, true_value, resultfile):
    assert(abs(mean-true_value)<1e-10*mean)
    print (f"true_SD[{repr(problem)},{n},{repr(name)}] = {sqrt(var)}",
           file = resultfile)
    resultfile.flush()

if __name__ == "__main__":
    print(f"KLEE-MINTY-BOX ##### {n = }, {num_samples} samples")
    true_value = 2**n # truth
    print(f"EXACT ANSWER: {2**n=}")

    fname = f"results-{problem}-{n}.py"
    resultfile = open(fname,"w")

    Me,Va = compute_variance()
    print(f"P    mean={Me}, variance={Va}={Va:g}, S.D.={sqrt(Va)}={sqrt(Va):g}")
    commit_to_file(Me,Va, 'Algorithm P:', true_value, resultfile)

    Me_IS,Va_IS = compute_variance_importance_sampling()
    print(f"IS   mean={Me_IS}, variance={Va_IS}={Va_IS:g}, S.D.={sqrt(Va_IS)}={sqrt(Va_IS):g}")
    commit_to_file(Me_IS,Va_IS, "Importance sampling (IS) by outdegree:", true_value, resultfile)

    Me_IS,Va_IS = compute_variance_importance_sampling(1)
    print(f"IS*  mean={Me_IS}, variance={Va_IS}={Va_IS:g}, S.D.={sqrt(Va_IS)}={sqrt(Va_IS):g}")
    commit_to_file(Me_IS,Va_IS, "Importance sampling (IS*) by reachable size:", true_value, resultfile)
    
    Me_IS,Va_IS = compute_variance_importance_sampling(2)
    print(f"IS** mean={Me_IS}, variance={Va_IS}={Va_IS:g}, S.D.={sqrt(Va_IS)}={sqrt(Va_IS):g}")
    commit_to_file(Me_IS,Va_IS, "Importance sampling (IS**) with size/indegree:", true_value, resultfile)

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

def heuristic_sampling(report_levels=False):
    """according to Pang Chen.
    stratum = outdegree
    """
    Store = [None]*(n+1)
    X = 0
    num_visited = 0
    num_inserted = 0
    Store[n] = source_vertex,1
    active = True
    while active:
        for outdeg in range(n,-1,-1):
            z = Store[outdeg] # expand highest outdegree
            #print(level,outdeg,z,Store)
            if z is not None:
                Store[outdeg] = None
                num_visited += 1
                # expand children
                p,W = z                
                X += W
                successors = outnbrs(p)
                num_inserted += len(successors)
                for i in successors:
                    q = swap(p,i)
                    outdeg_new = outdegree(q)
                    indeg_new = n-outdeg_new
                    insert(q, W/indeg_new, outdeg_new, Store)
                break
        else:
            active = False
    return X,num_visited,num_inserted                

def path_sampling():
    "Algorithm P"
    X = 0
    num_visited = 0
    u = source_vertex
    W = 1
    while 1:
        num_visited += 1
        X += W
        out = outnbrs(u)
        d_out = len(out)
        #print(u,out,indegree(u),f"{W=}, {X=}")
        if d_out==0:
            return X,num_visited
        W *= d_out
        i = choice(out)
        u = swap(u,i)
        W /= indegree(u)

def random_selection_from_class(S,k,n):
    # from the k-th class out of n almost equal classes
    bucket = S[k::n]
    vertices = [x for x,w in bucket]
    weights =  [w for x,w in bucket]
    return (choices(vertices, weights)[0], # weighted selection.
            sum(weights))       

def stochastic_estimation():
    "Algorithm SE" # horde sampling
    num_visited = 0
    num_inserted = 0
    num_condensed = 0 # after eliminating multiple equal elements
    X = 0
    A = [(source_vertex,1)]
    while A:
        new_A = defaultdict(float)
        for u,W in A:
            num_visited += 1
            X += W
            for i in outnbrs(u):
                num_inserted += 1
                v = swap(u,i)
                new_A[v] += W/indegree(v)
        num_condensed += len(new_A)
        new_A = list(new_A.items())
        if len(new_A) <= horde_limit:
            A = new_A
        else:
            A = [random_selection_from_class(new_A,k,horde_limit)
                 for k in range(horde_limit)]
    return X,num_visited,num_inserted, num_inserted-num_condensed

def importance_sampling(use_size=0):
    num_visited = 0
    X = 0
    u = source_vertex
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
            v = swap(u,i)
            if use_size==1:
                weight = 2**n-v # size of reachable subgraph
            if use_size==2:
                # also take into account the indegree
                weight = (2**n-v)/indegree(v)
            else:
                weight = max(outdegree(v), 1) # at least 1
                # increasing to 1 makes variance slightly worse
            prob_out.append(weight) 
        tot_prob = sum(prob_out)
        ind = choices(range(d), prob_out)[0] # weighted selection.
        i = succ[ind]
        u = swap(u,i)
        W *= tot_prob/(prob_out[ind]*indegree(u))

def importance_sampling_1():
    return importance_sampling(use_size=1)
def importance_sampling_2():
    return importance_sampling(use_size=2)

powers_of_ten = [10**i for i in range(1,20)]


def run_method(resultfile, problem,n,name,method,rand_seed,num_samples,true_value):
    if rand_seed is None:
        rand_seed = randint(0,10**9)
    seed(rand_seed) # can set a particular seed for random number generation

    print(name, f"random seed = {rand_seed}")        
    SX = SMax = 0
    Var0 = 0
    tot_visited = max_visited = tot_inserted = tot_saved = 0
    for m in range(1,1+num_samples):
        result = method()
        if len(result)==2:
            X,visited = method()
        elif len(result)==3:
            X,visited,inserted = result
            tot_inserted += inserted
        elif len(result)==4:
            X,visited,inserted,duplicates = result
            tot_inserted += inserted
            tot_saved += duplicates
        tot_visited += visited
        max_visited = max(max_visited,visited)
        SX += X
        SMax = max(SMax,X)
        Var0 += (X-true_value)**2
        #print(m,SX)
        if m in powers_of_ten:
            mean = SX/m
            Var_M = Var0/m
            Var_X = Var_M - (mean-true_value)**2
            CD = SMax/SX # Chatterjee-Diaconis estimated, between 0 and 1
            vis_statistics = (f"avg.#visited={tot_visited/m:0.2f}, "
                              f"max.#vis={max_visited}, ")
            if tot_inserted:
                vis_statistics += f"avg.#inserted={tot_inserted/m:0.2f}, "
            if tot_saved:
                vis_statistics += f"avg.#saved={tot_saved/m:0.2f}, "
            std_dev = sqrt(m/(m-1)*Var_X)
            print (f"{m:6} samples, {vis_statistics}"
                   f"{mean=:0.2f}±{std_dev:0.7g}, max={SMax:0.2f}, {CD=:0.5f}"
                   )
            sys.stdout.flush()
            print (f"""result[{repr(problem)},{n},{repr(name)},{m}] = {{
              'mean' : {mean}, 's.d.': {std_dev}, 'max': {SMax}, 'CD': {CD},
              'avg.#visited': {tot_visited/m}, 'max.#vis': {max_visited},
              "avg.#inserted": {tot_inserted/m:0.2f}, "avg.#saved": {tot_saved/m:0.2f} }}""",
                   file = resultfile)
            resultfile.flush()



            
    #print("true Var = variance with respect to the true mean of the distribution")

if __name__ == "__main__":
    print("first run P: ", path_sampling())
    print("first run HS:", heuristic_sampling())
    print("first run IS:", importance_sampling())
    print("first run SE:", stochastic_estimation())
    sys.stdout.flush()

    for name,method,rand_seed,num_s in [
        ("Algorithm P:", path_sampling, 755223781,num_samples),
        ("Importance sampling (IS) by outdegree:", importance_sampling, None,num_samples),
        ("Importance sampling (IS*) by reachable size:",
         importance_sampling_1, None,num_samples),
        ("Importance sampling (IS**) with size/indegree:",
         importance_sampling_2, None,num_samples),
        ("Heuristic sampling (HS) by outdegree:", heuristic_sampling, None,
           num_samples//10), # heuristic sampling does more work per sample.
        (f"Stochastic estimation (SE) ({horde_limit=}):", stochastic_estimation, None,num_samples//10),
          ]:
        run_method(resultfile, problem,n,name,method,rand_seed,num_s, true_value)
    resultfile.close()

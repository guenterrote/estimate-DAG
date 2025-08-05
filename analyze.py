from math import log,sqrt,exp, factorial
from collections import defaultdict

true_SD = {}
results = result = {}

filenames = """
results-KM-8.py  results-KM-10.py results-KM-12.py
results-KM-14.py results-KM-16.py results-KM-18.py
results-PERM-7.py  results-PERM-8.py  results-PERM-9.py
results-PERM-10.py results-PERM-11.py results-PERM-12.py 
"""
filenames = filenames.split()
#filenames = ["results-KM-10.py"]
for f in filenames:
    print("file",f)
    exec(open(f).read())
print(f"{len(results)=}")
#print(result.keys())

for x in true_SD.items():
    print("true_SD:", *x)

def pow10(m):
    k = int(log(m,10)+0.1)
    return f"10^{k}"

method = 'Algorithm P:'
tables = [ (method, "KM", "Klee--Minty cubes", "-P",
            [10,12]),
           (method,"PERM","Permutations", "-P",
            [8,9]),
           (method,"PERM","Permutations", "-P-all",
            range(100)),
           (method,"KM", "Klee--Minty cubes", "-P-all",
            range(100)),          
           ]
method,key = 'Importance sampling (IS) by outdegree:',"IS"
tables += [
    (method,"PERM","Permutations", f"-{key}-all",
     range(100)),
    (method,"KM", "Klee--Minty cubes",  f"-{key}-all",
     range(100)),
    (method, "KM", "Klee--Minty cubes", f"-{key}",
     [10,12]),
    (method,"PERM","Permutations", f"-{key}",
     [8,9]),
    ]
for method,key in [
        ('Stochastic estimation (SE) (horde_limit=10):', "CS"), # new name: Crowd Sampling 
        ('Heuristic sampling (HS) by outdegree:',"BS"), # new name: Bucket Sampling
        ]:
    tables += [
        (method, "KM", "Klee--Minty cubes", f"-{key}",
         [10,12]),
        (method,"PERM","Permutations", f"-{key}",
         [8,9]),
        (method,"PERM","Permutations", f"-{key}-all",
         range(100)),
        (method,"KM", "Klee--Minty cubes",  f"-{key}-all",
         range(100)),
    ]

for method,key in [
        ('Importance sampling (IS*) by reachable size:',"IS1"),
        ('Importance sampling (IS**) with size/indegree:',"IS2"),
        ]:
    tables += [
        (method,"KM", "Klee--Minty cubes",  f"-{key}-all",
         range(100)),
        (method, "KM", "Klee--Minty cubes", f"-{key}",
         [10,12]),
    ]

for meth0,prob0,PROB,name_suffix,list_n in tables:
    tables = defaultdict(list)
    tableslatex = defaultdict(list)
    headers = defaultdict(str)
    avg_visited = dict() # later entries overwrite earlier entries
    max_visited = dict()
    avg_inserted = dict()
    avg_saved = dict()
    latexfile = f"../tableP-{prob0}.tex"
    latexfile = f"../table-{prob0}{name_suffix}.tex"
    latexfile = f"tables/table-{prob0}{name_suffix}.tex"
    out = open(latexfile,"w")
    
    for (prob,n,meth,m),v in results.items():
#        if n not in list_n: #(prob,n) in (("KM",18),("KM",16),("PERM",10),("PERM",11),("PERM",12),("KM",14)) or prob!=prob0:
#            continue
        if meth==meth0 and prob==prob0 and n in list_n:
            headers[prob,n] = f"{prob} {n=}"
            line = (f" | {v['mean']:0,.1f}±{v['s.d.']:0.2e}, max={v['max']:0.2e}, CD={v['CD']:0.5f} |"
                       )
            stddev = v['s.d.']
            avg_visited[n] = v['avg.#visited']
            max_visited[n] = v['max.#vis']
            avg_inserted[n] = v['avg.#inserted']
            avg_saved[n] = v['avg.#saved']
            
            linelatex = (f"&${v['mean']:0,.1f}±" +
                         (f"{stddev:0,.1f}" if stddev <10**5 else f"{stddev:0.2e}") +
                         f"$& ${v['max']:0.2e}$ & ${v['CD']:0.5f}$"
                       )
            tables[prob,n].append(pow10(m) + line)
            tableslatex[prob,n].append((m,linelatex))
    
    numcolumns = len(headers)
    
    def prepareline(l):
        o = ""
        i = 0
        while i < len(l):
            c = l[i]
            if c=="±":
                o += r"$&$\,\pm\,"
            elif c==",":
                o += "{,}"
            elif c=="e" and l[i+1] in "+-":
#                o += fr"{{\times}}\!10^{{{int(l[i+1:i+4])}}}"
                if l[i+4] in "012345789":
                    exp_length = 3
                else:
                    exp_length = 2
                o += fr"{{\times}}10^{{{int(l[i+1:i+2+exp_length])}}}"
                i += exp_length+1
            else:
                o += c
            i += 1
        return o
                
    longlines = defaultdict(str)
    groupheaders = ""
    truth = ""
    visitedline = ""
    generatedline = ""
    savedline = ""
    
    for k,head in sorted(headers.items()):
        print(head)
        prob,n = k
        groupheaders += r"&\multicolumn4{c|}{$n=%d$}"%n
        if prob=="KM":
            truevalue = 2**n
        else:
            truevalue = factorial(n)
        #print("S",avg_visited[n],max_visited[n])
        if avg_visited[n]==max_visited[n]:
            vis = rf'\ \ \ \#visited = {max_visited[n]:d}'
        else:
            vis = ""
            visitedline += rf'&\multicolumn4{{l|}}{{\#visited: avg\ = {avg_visited[n]:0.1f}, max\ = {max_visited[n]:d}}}'
        if avg_inserted[n]:
            generatedline += rf'&\multicolumn4{{l|}}{{\#generated: avg\ = {avg_inserted[n]:0.1f}}}'
        if avg_saved[n]:
            percent = avg_saved[n]/avg_inserted[n]*100
            savedline += rf'&\multicolumn4{{l|}}{{\#saved: avg\ = {avg_saved[n]:0.1f} = {percent:0.1f}\%}}'
        SD_value = true_SD.get((prob,n,meth0))
        if SD_value is not None:
            t_SD = fr"±\rlap{{${SD_value:1.3e}$"
        else:
            t_SD = r"$&$\rlap{"
        truth += fr"&${truevalue:0,d}\phantom{{.0}}{t_SD}{vis}}}$&&"
        print ("Truth",truth)
        print("\n".join(tables[k]))
        for m,line in tableslatex[k]:
            longlines[m] += prepareline(line)
    
    entry = r"\par\noindent\begin{tabular}{|c" + "|r@{}lll" *numcolumns + "|}" + "\n\\hline"
    print (entry, file=out)
    print (rf"\multicolumn{{{1+4*numcolumns}}}{{|c|}}{{{PROB}}}\\"+ "\n\\hline", file=out)
    print (groupheaders + r'\\', file=out)
    header = ("$m$ "+ r"& mean &\,$\pm$\,std.\,dev.&\multicolumn 1 c{max.}&\multicolumn 1 {c|}{$\hat\chi$}"*numcolumns
              + r"\\" + "\n\\hline" + r"\vbox to 11pt{}%")
    print (header, file=out)
    for m,lo in sorted(longlines.items()):
        m2 = pow10(m) # only single-digit powers so far
        print(f"${m2}${lo}  \\\\", file=out)
    print("\\hline\n" + r"\vbox to 11pt{}" + prepareline(truth) + r'\\', file=out)
    if visitedline:
        print(visitedline + r"\\", file=out)        
    if generatedline:
        print(generatedline + r"\\", file=out)        
    if savedline:
        print(savedline + r"\\", file=out)        
    print("\\hline\n" + r"\end{tabular}", file=out)
    out.close()
        
print(sorted(set((prob,meth) for (prob,n,meth,m) in results)))
print(sorted(set(meth for (prob,n,meth,m) in results)))
    

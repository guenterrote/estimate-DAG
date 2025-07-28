from math import log,sqrt,exp, factorial
from collections import defaultdict

true_SD = {}
results = result = {}
filenames = """
results-KM-10.py  results-KM-14.py  results-KM-18.py	results-PERM-11.py  results-PERM-8.py
results-KM-12.py  results-KM-16.py  results-PERM-10.py	results-PERM-12.py  results-PERM-9.py""".split()
filenames = ["results-KM-10.py"]
for f in filenames:
    print("file",f)
    exec(open(f).read())
print(len(results))
#print(result.keys())

for x in true_SD.items():
    print("true_SD:", *x)

def pow10(m):
    k = int(log(m,10)+0.1)
    return f"10^{k}"

meth0 = 'Algorithm P:'
for prob0,PROB in (("KM", "Klee--Minty cubes"), ("PERM","Permutations")):
    tables = defaultdict(list)
    tableslatex = defaultdict(list)
    headers = defaultdict(str)
    avg_visited = dict() # later entries overwrite earlier entries
    max_visited = dict()
    latexfile = f"../tableP-{prob0}.tex"
    latexfile = f"tableP-{prob0}.tex"
    out = open(latexfile,"w")
    
    
    for (prob,n,meth,m),v in results.items():
        if (prob,n) in (("KM",18),("KM",16),("PERM",10),("PERM",11),("PERM",12),("KM",14)) or prob!=prob0:
            continue
        if meth==meth0:
            headers[prob,n] = f"{prob} {n=}"
            line = (f" | {v['mean']:0,.1f}±{v['s.d.']:0.2e}, max={v['max']:0.2e}, CD={v['CD']:0.4f} |"
                       )
            stddev = v['s.d.']
            avg_visited[n] = v['avg.#visited']
            max_visited[n] = v['max.#vis']
            
            linelatex = (f"&${v['mean']:0,.1f}±" +
                         (f"{stddev:0,.1f}" if stddev <10**5 else f"{stddev:0.2e}") +
                         f"$& ${v['max']:0.2e}$ & ${v['CD']:0.4f}$"
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
                o += fr"{{\times}}\!10^{{{int(l[i+1:i+4])}}}"
                i += 3
            else:
                o += c
            i += 1
        return o
                
    longlines = defaultdict(str)
    groupheaders = ""
    truth = ""
    visitedline = ""
    
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
            visitedline += rf'&\multicolumn4{{l|}}{{\#visited: avg={avg_visited[n]:0.1f}, max={max_visited[n]:d}}}'
        truth += fr"&${truevalue:0,d}\phantom{{.0}}±\rlap{{${true_SD[prob,n]:1.3e}${vis}}}$&&"
        print ("T",truth)
        print("\n".join(tables[k]))
        for m,line in tableslatex[k]:
            longlines[m] += prepareline(line)
    
    entry = r"\par\noindent\begin{tabular}{|r" + "|r@{}lll" *numcolumns + "|}" + "\n\\hline"
    print (entry, file=out)
    print (rf"\multicolumn{{{1+4*numcolumns}}}{{|c|}}{{{PROB}}}\\"+ "\n\\hline", file=out)
    print (groupheaders + r'\\', file=out)
    header = ("$m$ "+ r"& mean &\,$\pm$\,std.dev.&\multicolumn 1 c{max.}&\multicolumn 1 {c|}{$\hat\chi$}"*numcolumns
              + r"\\" + "\n\\hline" + r"\vbox to 11pt{}%")
    print (header, file=out)
    for m,lo in sorted(longlines.items()):
        m2 = pow10(m) # only single-digit powers so far
        print(f"${m2}${lo}  \\\\", file=out)
    print("\\hline\n" + prepareline(truth) + r'\\', file=out)
    if visitedline:
        print(visitedline  + r'\\', file=out)        
    print("\\hline\n" + r"\end{tabular}", file=out)
    out.close()
        
        
             
    

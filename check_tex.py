import re

with open("I:/Scratch/tough-math/erdos-straus/paper/erdos_straus_gateway.tex") as f:
    txt = f.read()

pat_env = re.compile(r"\\(begin|end)\{(\w+)\}")
stack = []
for m in pat_env.finditer(txt):
    action, name = m.group(1), m.group(2)
    if action == "begin":
        stack.append(name)
    else:
        if stack and stack[-1] == name:
            stack.pop()
        else:
            print("MISMATCH:", name)
if stack:
    print("UNCLOSED:", stack)
else:
    print("All environments balanced")

pat_cite = re.compile(r"\\cite\{([^}]+)\}")
pat_bib = re.compile(r"\\bibitem\{([^}]+)\}")
cites = set()
for m in pat_cite.finditer(txt):
    for k in m.group(1).split(","):
        cites.add(k.strip())
bibs = set(b.group(1) for b in pat_bib.finditer(txt))
miss = cites - bibs
if miss:
    print("MISSING bibitems:", miss)
else:
    print("All cited keys have bibitems")
unused = bibs - cites
if unused:
    print("UNUSED bibitems:", unused)

pat_ref = re.compile(r"\\ref\{([^}]+)\}")
pat_lab = re.compile(r"\\label\{([^}]+)\}")
refs = set(r.group(1) for r in pat_ref.finditer(txt))
labels = set(l.group(1) for l in pat_lab.finditer(txt))
miss2 = refs - labels
if miss2:
    print("MISSING labels:", miss2)
else:
    print("All refs have labels")
print("Lines:", len(txt.splitlines()))

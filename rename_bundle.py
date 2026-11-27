"""Rename the bundle everywhere.
Usage: python rename_bundle.py NewName "New Display Name" [new_command_prefix]
Example: python rename_bundle.py PilotDock "Lattice PilotDock" pdock
Update CURRENT below if you edit names by hand.
"""
import os, re, sys
CURRENT = {"name": "VinaPilotX", "display": "Lattice VinaPilotX", "cmd": "vinax"}

if len(sys.argv) not in (3, 4):
    sys.exit(__doc__)
new, new_display = sys.argv[1], sys.argv[2]
new_cmd = sys.argv[3] if len(sys.argv) == 4 else CURRENT["cmd"]
cur = CURRENT
subs = [(cur["display"], new_display), (cur["name"], new), (cur["name"].lower(), new.lower())]
cmd_re = re.compile(r"(?<![\w-])%s (gui|vina|box|run|cancel|pdbqt)\b" % re.escape(cur["cmd"]))
root = os.path.dirname(os.path.abspath(__file__))
for d, _, files in os.walk(root):
    if any(p in d for p in (".git", "__pycache__", ".pytest_cache")):
        continue
    for f in files:
        if f == "rename_bundle.py" or f.startswith("LICENSE") or not f.endswith((".py", ".xml", ".md", ".html", ".txt")):
            continue
        p = os.path.join(d, f)
        t = open(p, encoding="utf-8").read(); o = t
        for a, b in subs:
            t = t.replace(a, b)
        t = cmd_re.sub(new_cmd + r" \1", t)
        if t != o:
            open(p, "w", encoding="utf-8").write(t)
old_html = os.path.join(root, "src", "docs", "user", "tools", cur["name"].lower() + ".html")
if os.path.exists(old_html) and new.lower() != cur["name"].lower():
    os.rename(old_html, os.path.join(os.path.dirname(old_html), new.lower() + ".html"))
# record the new current names inside this script
s = open(__file__, encoding="utf-8").read()
s = re.sub(r'CURRENT = \{.*\}', 'CURRENT = {"name": "%s", "display": "%s", "cmd": "%s"}' % (new, new_display, new_cmd), s, count=1)
open(__file__, "w", encoding="utf-8").write(s)
print("Renamed to", new, "/", new_display, "/ command:", new_cmd)

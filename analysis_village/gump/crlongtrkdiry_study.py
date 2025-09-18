from gump_util import *
from gump_cuts import *

df_nd, _, _ = load_data("/exp/sbnd/app/users/nrowe/cafpyana/sbnd_no_cuts.df")
df_nd = df_nd[df_nd.crlongtrkdiry > -1000]
vals = breakdown_top(df_nd.crlongtrkdiry, df_nd)
n, b, _ = plt.hist(vals, bins=20, stacked=True, label=top_labels)
og_sig = sum((vals[0]))
best_pur = 0.0

for bi in b[:-1]:
    sig = sum((vals[0] > bi))
    tot = sum((df_nd.crlongtrkdiry > bi))
    pur = sig/tot
    eff = sig/og_sig
    print(f"cut on {bi}")
    print(f"purity: {pur}")
    print(f"efficiency: {eff}")
    if pur > best_pur:
        best_pur = pur
        best_bin = bi

plt.xlabel("crlongtrkdiry")
plt.ylabel("Events")
plt.vlines([best_bin], 0, max(n[-1])*1.1, color='red')
#plt.text(0.5, max(n[-1])/2, str(best_pur))
plt.legend()
plt.savefig('crlongtrkdiry.png')
plt.clf()

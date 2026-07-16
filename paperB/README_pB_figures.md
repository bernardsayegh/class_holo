# Paper B figure scripts

| Figure | File produced | Script | Inputs |
|---|---|---|---|
| 1 | pB_fig1_capacity_history.pdf | xh_posterior_freeze.py | frozen non-ladder chain (modelB_A2_noprior_v5_acc.*) |
| 2 | pB_fig2_pushforward.pdf | xh_posterior_freeze.py | same run (both figures from one invocation) |
| 3 | pB_fig3_native_plane.pdf | pB_fig3_native_plane.py | none (analytic) |
| 4 | pB_fig4_closure_loop.pdf | pB_fig4_closure_loop.py | none (schematic) |

Run:
```
python3 xh_posterior_freeze.py chains_test/modelB_A2_noprior_v5_acc.1.txt
python3 pB_fig3_native_plane.py
python3 pB_fig4_closure_loop.py
```
xh_posterior_freeze.py merges all parallel chains (30% burn-in each), thins
to ~400 samples, re-integrates the interacting background per sample, and
self-certifies against the chain's stored H0_local (achieved: max relative
deviation 7.7e-5). Reported freeze: <X_H> = 0.03592, sigma(X_H) < 1e-5
(flow invariant), z_entry = 0.640(8), z_exit = 0.014(5).

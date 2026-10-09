import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

DIR = "/Users/xzhao/Downloads/github/lr-annotation/lrGS_coverage_blacklist_analyses"

p90_bp = 172409732
p70_bp = 178705532
p50_bp = 185630332

# perfect nesting confirmed via bedtools subtract (p90 subset of p70 subset of p50)
only_p50 = p50_bp - p70_bp
only_p70 = p70_bp - p90_bp
only_p90 = p90_bp

fig, ax = plt.subplots(figsize=(7, 7))
# venn3 subset order: (Abc, aBc, ABc, abC, AbC, aBC, ABC); A=p50, B=p70, C=p90
subsets = (only_p50, 0, only_p70, 0, 0, 0, only_p90)
v = venn3(subsets=subsets, set_labels=("blacklist >=50%", "blacklist >=70%", "blacklist >=90%"), ax=ax)

label_map = {"100": ("only >=50%", only_p50), "110": ("only >=70% (& >=50%)", only_p70), "111": (">=90% (& >=70%, >=50%)", only_p90)}
# rings are too thin (only ~4% difference) for default label placement to avoid
# overlap; pin each label to a distinct spot with a leader line instead.
positions = {"100": (-1.55, 0.55), "110": (-1.55, 0.15), "111": (-1.55, -0.25)}
for label_id, (desc, bp) in label_map.items():
    lbl = v.get_label_by_id(label_id)
    if lbl:
        orig_xy = lbl.get_position()
        lbl.set_visible(False)
        new_x, new_y = positions[label_id]
        ax.annotate(f"{desc}:\n{bp / 1e6:.2f} Mb", xy=orig_xy, xytext=(new_x, new_y),
                    fontsize=9, ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.8))

ax.set_xlim(-2.0, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.set_title("hgsvc_hprc long-read low-coverage blacklist:\noverlap between sample-fraction cutoffs\n(perfectly nested: p90 ⊆ p70 ⊆ p50)")
fig.tight_layout()
fig.savefig(f"{DIR}/hgsvc_hprc.blacklist.cutoff_venn.png", dpi=150)
print("saved", f"{DIR}/hgsvc_hprc.blacklist.cutoff_venn.png")

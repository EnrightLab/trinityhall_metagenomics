import matplotlib.pyplot as plt
import pandas as pd

# List of files and sample names
files = [
    ('barcode13.bracken.microb.ncbi.txt', 'Sample 1'),
    ('barcode14.bracken.microb.ncbi.txt', 'Sample 2'),
    ('barcode15.bracken.microb.ncbi.txt', 'Sample 3'),
    ('barcode16.bracken.microb.ncbi.txt', 'Sample 4')
]

# High-contrast custom color list (20 colors)
custom_colors = [
    "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
    "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000",
    "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080"
]

# Create plots
fig, axes = plt.subplots(1, 4, figsize=(20, 10), sharey=True)

for i, (filepath, sample_label) in enumerate(files):
    ax = axes[i]

    # Read and get top 20 species
    df = pd.read_csv(filepath, sep='\t')
    df = df.sort_values(by='new_est_reads', ascending=False).head(20)
    sample_data = df[['name', 'new_est_reads']].set_index('name')
    species = sample_data.index.tolist()
    values = sample_data['new_est_reads'].tolist()

    # Assign fixed high-contrast colors
    colors = custom_colors[:len(species)]

    bottoms = 0
    for sp, val, color in zip(species[::-1], values[::-1], colors[::-1]):  # reverse for top-down order
        bar = ax.bar(0, val, bottom=bottoms, color=color, label=sp)
        ax.text(
            x=bar[0].get_x() + bar[0].get_width() / 2,
            y=bottoms + val / 2,
            s=f'{val:.0f}',
            ha='center',
            va='center',
            fontsize=6,
            color='black'
        )
        bottoms += val

    ax.set_title(sample_label)
    ax.set_xticks([0])
    ax.set_xticklabels([sample_label])
    if i == 0:
        ax.set_ylabel('Abundance')

    ax.legend(
        loc='center left',
        bbox_to_anchor=(1.0, 0.5),
        title='Species',
        fontsize=6,
        title_fontsize=8
    )

# Adjust layout
plt.suptitle('Top 20 Species per Sample (Highest on Top)', fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.88, right=0.95)
plt.savefig("top20_species.png", dpi=300, bbox_inches='tight')
plt.show()

import matplotlib.pyplot as plt
import pandas as pd

# List of files and sample names
files = [
    ('barcode13.G.bracken.microb.txt', 'Sample 1'),
    ('barcode14.G.bracken.microb.txt', 'Sample 2'),
    ('barcode15.G.bracken.microb.txt', 'Sample 3'),
    ('barcode16.G.bracken.microb.txt', 'Sample 4')
]

# High-contrast custom color list (20 colors)
custom_colors = [
    "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
    "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000",
    "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080"
]

# Create plots
fig, axes = plt.subplots(1, 4, figsize=(20, 10), sharey=True)

for i, (filepath, sample_label) in enumerate(files):
    ax = axes[i]

    # Read and get top 20 species
    df = pd.read_csv(filepath, sep='\t')
    df = df.sort_values(by='new_est_reads', ascending=False).head(20)
    sample_data = df[['name', 'new_est_reads']].set_index('name')
    species = sample_data.index.tolist()
    values = sample_data['new_est_reads'].tolist()

    # Assign fixed high-contrast colors
    colors = custom_colors[:len(species)]

    bottoms = 0
    for sp, val, color in zip(species[::-1], values[::-1], colors[::-1]):  # reverse for top-down order
        bar = ax.bar(0, val, bottom=bottoms, color=color, label=sp)
        ax.text(
            x=bar[0].get_x() + bar[0].get_width() / 2,
            y=bottoms + val / 2,
            s=f'{val:.0f}',
            ha='center',
            va='center',
            fontsize=6,
            color='black'
        )
        bottoms += val

    ax.set_title(sample_label)
    ax.set_xticks([0])
    ax.set_xticklabels([sample_label])
    if i == 0:
        ax.set_ylabel('Abundance')

    ax.legend(
        loc='center left',
        bbox_to_anchor=(1.0, 0.5),
        title='Genus',
        fontsize=6,
        title_fontsize=8
    )

# Adjust layout
plt.suptitle('Top 20 Genus per Sample (Highest on Top)', fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.88, right=0.95)
plt.savefig("top20_genus.png", dpi=300, bbox_inches='tight')
plt.show()
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.patches import Patch

refseqV205 = pd.read_csv('kraken2_refseqV205_compare_with_truth.csv', index_col=0, header=0)
test_results = pd.read_csv('test_results_compare_with_truth.csv', index_col=0, header=0)

samples, confidence = [], [str(i/100).ljust(4, '0') for i in range(0,101) if i % 5 == 0]
for row in refseqV205.index:
  samples.append(row.split('-')[0])
  
samples = list(set(samples))
metrics = refseqV205.columns
metric_names = {'proportion_classified':'Proportion of reads classified', 'precision_taxa':'Precision taxa', 'recall_taxa':'Recall taxa', 'f1_taxa':'F1 score taxa', 'precision_reads':'Precision reads', 'recall_reads':'Recall reads', 'f1_reads':'F1 score reads', 'dissimilarity_species':'Dissimilarity in number of species identified', 'shannon_diversity':'Dissimilarity in Shannon diversity', 'aitchisons_distance':"Robust Aitchison's distance", 'braycurtis_distance':'Bray-Curtis distance', 'L1_distance':'L1 distance'}

fig = plt.figure(figsize=(20,20))
count = 0
for metric in metrics:
  ax = plt.subplot(4,3,count+1)
  plot_count = 0
  for df in [refseqV205, test_results]:
    vals_conf = [[] for c in confidence]
    for c in range(len(confidence)):
      conf = confidence[c]
      for sample in samples:
        for row in df.index.values:
          if conf in row and sample in row:
            vals_conf[c].append(df.loc[row, metric])
    mean_conf, upper_conf, lower_conf, plot_conf = [], [], [], []
    for v in range(len(vals_conf)):
      if len(vals_conf[v]) > 0:
        mean_conf.append(np.median(vals_conf[v]))
        upper_conf.append(np.percentile(vals_conf[v], 75))
        lower_conf.append(np.percentile(vals_conf[v], 25))
        plot_conf.append(float(confidence[v]))
    if plot_count == 0: color = 'k'
    else: color = 'firebrick'
    line = ax.plot(plot_conf, mean_conf, color=color)
    line = ax.fill_between(plot_conf, upper_conf, lower_conf, color=color, alpha=0.2)
    plot_count += 1
  xt = plt.xticks([float(conf) for conf in confidence], confidence, rotation=90)
  xl = plt.xlabel('Confidence threshold')
  xl = plt.xlim([-0.025, 1.025])
  ti = plt.title(metric_names[metric], fontweight='bold')
  if count == 2:
    handles = [Patch(facecolor='k', edgecolor='k', label='RefSeq Complete V205'), Patch(facecolor='firebrick', edgecolor='k', label='Test database')]
    leg = ax.legend(handles=handles, bbox_to_anchor=(0.5,1.1), loc='lower center')
  count += 1

plt.tight_layout()
plt.savefig('Comparison_of_databases.png', dpi=600, bbox_inches='tight')

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.patches import Patch

refseqV205 = pd.read_csv('kraken2_refseqV205_compare_with_truth.csv', index_col=0, header=0)
metaphlan3 = pd.read_csv('metaphlan3_compare_with_truth.csv', index_col=0, header=0)
test_results = pd.read_csv('test_results_compare_with_truth.csv', index_col=0, header=0)

samples, confidence = [], [str(i/100).ljust(4, '0') for i in range(0,101) if i % 5 == 0]
for row in refseqV205.index:
  samples.append(row.split('-')[0])

samples = list(set(samples))
metrics = refseqV205.columns
metric_names = {'proportion_classified':'Proportion of reads classified', 'precision_taxa':'Precision taxa', 'recall_taxa':'Recall taxa', 'f1_taxa':'F1 score taxa', 'precision_reads':'Precision reads', 'recall_reads':'Recall reads', 'f1_reads':'F1 score reads', 'dissimilarity_species':'Dissimilarity in number of species identified', 'shannon_diversity':'Dissimilarity in Shannon diversity', 'aitchisons_distance':"Robust Aitchison's distance", 'braycurtis_distance':'Bray-Curtis distance', 'L1_distance':'L1 distance'}

fig = plt.figure(figsize=(30,20))
count1, count2 = 0, 0
for metric in metrics:
  if count1 == 12: 
    count1 = 0
    count2 += 1
  ax_krak = plt.subplot2grid((4,12), (count2,count1), colspan=3)
  ax_met = plt.subplot2grid((4,12), (count2,count1+3), sharey=ax_krak)
  count1 += 4
  plot_count = 0
  krak_plot = [refseqV205]
  if '-' in test_results.index[0]: krak_plot.append(test_results)
  for df in krak_plot:
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
    line = ax_krak.plot(plot_conf, mean_conf, color=color)
    line = ax_krak.fill_between(plot_conf, upper_conf, lower_conf, color=color, alpha=0.2)
    plot_count += 1
  plt.sca(ax_krak)
  xt = plt.xticks([float(conf) for conf in confidence], confidence, rotation=90)
  xl = plt.xlabel('Confidence threshold')
  xl = plt.xlim([-0.025, 1.025])
  ti = plt.title(metric_names[metric], fontweight='bold')
  if count1 == 12 and count2 == 0:
    handles = [Patch(facecolor='k', edgecolor='k', label='RefSeq Complete V205'), Patch(facecolor='gray', edgecolor='k', label='MetaPhlAn 3'), Patch(facecolor='firebrick', edgecolor='k', label='Test database/tool')]
    leg = ax_krak.legend(handles=handles, bbox_to_anchor=(0.5,1.1), loc='lower center')
  met_plot = [metaphlan3]
  if '-' not in test_results.index[0]: met_plot.append(test_results)
  count = 0
  colors = ['gray', 'firebrick']
  plt.sca(ax_met)
  for df in met_plot:
    df = df.loc[samples, :]
    vals = df.loc[:, metric].values
    sc = ax_met.scatter(np.random.normal(count, scale=0.075, size=len(vals)), vals, color=colors[count], alpha=0.5)
    box = ax_met.boxplot(vals, positions=[count], showfliers=False, widths=0.5, vert=True)
    for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']: plt.setp(box[item], color='k')
    count += 1
  if 'dissimilarity' in metric or 'shannon' in metric:
    plt.yscale('symlog')

plt.tight_layout()
plt.savefig('Comparison_of_databases.png', dpi=600, bbox_inches='tight')

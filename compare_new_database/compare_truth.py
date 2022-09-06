import os
import pandas as pd
import numpy as np
from skbio.diversity import alpha_diversity, beta_diversity
from deicode.preprocessing import rclr
from scipy.spatial import distance

truth = pd.read_csv('truth.csv', index_col=0, header=0)
# results = pd.read_csv('subset_kraken2_refseqV205_combined.csv', index_col=0, header=0)
results = pd.read_csv('bracken_out.csv', index_col=0, header=0)
col_rename = {}
for col in results.columns:
  if len(col.split('.',1)[1]) == 4:
    col_rename[col] = col.split('.',1)[0]+'-'+col.split('.',1)[1]

results = results.rename(columns=col_rename)

output = 'test_results_compare_with_truth.csv'

truth.index = truth.index.map(int)
results.index = results.index.map(int)

metrics = ['sample', 'proportion_classified', 'precision_taxa', 'recall_taxa', 'f1_taxa', 'precision_reads', 'recall_reads', 'f1_reads', 'dissimilarity_species', 'shannon_diversity', 'aitchisons_distance', 'braycurtis_distance', 'L1_distance']

def prec_rec_f1(sample_df):
  y_true, y_pred = sample_df.iloc[:, 0].values, sample_df.iloc[:, 1].values
  correct_reads, correct_taxa = 0, 0
  for v in range(len(y_true)):
    if y_true[v] > 0 and y_pred[v] > 0:
      if y_pred[v] >= y_true[v]: correct_reads += y_true[v]
      else: correct_reads += y_pred[v]
      if y_pred[v] > 0 and y_true[v] > 0: correct_taxa += 1
  precision_reads = correct_reads/sum(y_pred) 
  recall_reads = correct_reads/sum(y_true)
  if precision_reads == 0 or recall_reads == 0: f1_score_reads = 0
  else: f1_score_reads = 2*((precision_reads*recall_reads)/(precision_reads+recall_reads))
  precision_taxa = correct_taxa/sum([1 for y in y_pred if y > 0])
  recall_taxa = correct_taxa/sum([1 for y in y_true if y > 0])
  if precision_taxa == 0 or recall_taxa == 0: f1_score_taxa = 0
  else: f1_score_taxa = 2*((precision_taxa*recall_taxa)/(precision_taxa+recall_taxa))
  return precision_taxa, recall_taxa, f1_score_taxa, precision_reads, recall_reads, f1_score_reads

metric_results = []
for col in results.columns:
  #if 'ani95_cLOW_stFalse_r0' not in col: continue
  this_sample = [col]
  ts = col.split('-')[0]
  sam_df = results.copy(deep=True).loc[:, [col]]
  tru_df = truth.copy(deep=True).loc[:, [ts]]
  sam_df = sam_df[sam_df.max(axis=1) > 0]
  tru_df = tru_df[tru_df.max(axis=1) > 0]
  comb_df = pd.concat([sam_df, tru_df]).fillna(value=0)
  comb_df = comb_df.groupby(by=comb_df.index, axis=0).sum()
  #proportion classified
  this_sample.append(sam_df.sum(axis=0).values[0]/tru_df.sum(axis=0).values[0])
  #precision, recall, f1
  if this_sample[0] == 0: pt, rt, ft, pr, rr, fr = 0, 0, 0, 0, 0, 0
  else: pt, rt, ft, pr, rr, fr = prec_rec_f1(comb_df)
  this_sample = this_sample+[pt, rt, ft, pr, rr, fr]
  #dissimilarity species
  sam_sp, tru_sp = sam_df.shape[0], tru_df.shape[0]
  this_sample.append(sam_sp-tru_sp)
  #dissimilarity shannon
  sam_shan = list(alpha_diversity('shannon', list(sam_df.loc[:, col].values)))[0]
  tru_shan = list(alpha_diversity('shannon', list(tru_df.loc[:, ts].values)))[0]
  this_sample.append(sam_shan-tru_shan)
  #aitchisons distance
  X = comb_df.iloc[0:].values
  rclr_sample = rclr(X)
  rclr_sample = pd.DataFrame(rclr_sample, columns=comb_df.columns, index=comb_df.index.values).fillna(value=0)
  X = rclr_sample.transpose().iloc[0:].values
  similarities = np.nan_to_num(distance.cdist(X, X, 'euclidean'))
  this_sample.append(similarities[0][1])
  #braycurtis distance
  X = comb_df.transpose().iloc[0:].values
  similarities = np.nan_to_num(distance.cdist(X, X, 'braycurtis')) 
  this_sample.append(similarities[0][1])
  #l1 distance
  X = comb_df.transpose().iloc[0:].values
  similarities = np.nan_to_num(distance.cdist(X, X, 'cityblock')) 
  this_sample.append(similarities[0][1])
  #end of calculations
  metric_results.append(this_sample)

all_results = pd.DataFrame(metric_results, columns=metrics).fillna(value=0).set_index('sample')
all_results.to_csv(output)

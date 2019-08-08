import matplotlib
matplotlib.use('agg')
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats


def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2

# CheckM stdout and sag-plus-err-0.1.py output
checkm_stdout_path = '/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_1/2/checkM/checkM_stdout.tsv'
err_stats_path = '/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_1/2/error_analysis/All_stats.tsv'

checkm_df = pd.read_csv(checkm_stdout_path, sep='\t', header=0)
# extract sag ID from bin id
checkm_df['sag_id'] = [x.rsplit('.', 1)[0] for x in checkm_df['Bin Id']]

# reshape sag+ err output df
err_stats_df = pd.read_csv(err_stats_path, sep='\t', header=0)
algorithm_list = list(set(err_stats_df['algorithm']))
level_list = list(set(err_stats_df['level']))
df_list = []
for algorithm in algorithm_list:
	algorithm_df = err_stats_df.loc[err_stats_df['algorithm'] == algorithm]
	for level in level_list:
		level_df = algorithm_df.loc[algorithm_df['level'] == level]
		pivot_df = level_df[['sag_id', 'statistic', 'score']
							].pivot(index='sag_id', columns='statistic', values='score'
								).reset_index()
		pivot_df['algorithm'] = algorithm
		pivot_df['level'] = level
		df_list.append(pivot_df)
concat_df = pd.concat(df_list)

# merge checkm and sag+ err
merge_df = concat_df.merge(checkm_df, on='sag_id', how='left')
merge_df.to_csv('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_1/2/error_analysis/' +
				'/Merged_checkM_errstat.tsv', index=False, sep='\t')
# create 1-precision col
merge_df['1-precision'] = 1 - merge_df['precision']
combined_df = merge_df.loc[(merge_df['algorithm'] == 'combined') &
							(merge_df['level'] == 'species')
							]
sns.set_context("paper")
ax = sns.jointplot(x='Completeness', y='sensitivity', data=combined_df,
					kind="reg", stat_func=r2
					)
plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_1/2/error_analysis/' +
			'Compl_v_sensi_scatter.png',
			bbox_inches='tight'
			)
plt.clf()

sns.set_context("paper")
ax = sns.jointplot(x='Contamination', y='1-precision', data=combined_df,
					kind="reg", stat_func=r2
					)
plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_1/2/error_analysis/' +
			'Conta_v_invprec_scatter.png',
			bbox_inches='tight'
			)
plt.clf()


contam_df = combined_df[['sag_id', 'algorithm', 'level', 'Contamination']]
#contam_df['assessment_type'] = 'checkM'

precision_df = combined_df[['sag_id', 'algorithm', 'level',  '1-precision']]
precision_df['1-precision'] = precision_df['1-precision']*100
precision_df.columns = ['sag_id', 'algorithm', 'level',  '1-precision']

#precision_df['assessment_type'] = 'true'
#contam_prec_df = pd.concat([contam_df, precision_df])

comple_df = combined_df[['sag_id', 'algorithm', 'level',  'Completeness']]
#comple_df['assessment_type'] = 'checkM'

sensitivity_df = combined_df[['sag_id', 'algorithm', 'level',  'sensitivity']]
sensitivity_df['sensitivity'] = sensitivity_df['sensitivity']*100
sensitivity_df.columns = ['sag_id', 'algorithm', 'level',  'sensitivity']

#sensitivity_df['assessment_type'] = 'true'
#compl_sens_df = pd.concat([comple_df, sensitivity_df])

merge_concomp_df = contam_df.merge(comple_df, on=['sag_id', 'algorithm', 'level'],
									how='left'
									)
merge_presen_df = precision_df.merge(sensitivity_df, on=['sag_id', 'algorithm', 'level'],
									how='left'
									)

# checkm
sns.set_context("paper")
sns.set(font_scale=1.5)

#ax = sns.scatterplot(x='Completeness', y='Contamination', hue='assessment_type',
#						data=combined_contam_compl_df
#						)
#ax = sns.jointplot(x='Completeness', y='Contamination', data=merge_concomp_df)
#ax.ax_marg_x.set_xlim(-5, 105)
#ax.ax_marg_y.set_ylim(-5, 105)
#ax.ax_joint.plot([50.0, 50.0], [0.0, 100.0], ':k')
#ax.ax_joint.plot([0.0, 100.0], [50.0, 50.0], ':k')
ax = sns.kdeplot(merge_concomp_df['Completeness'].dropna(), color='blue', label='CheckM',
					bw=2, shade=True, legend=False
					)
ax = sns.kdeplot(merge_presen_df['sensitivity'].dropna(), color='orange', label='True',
					bw=2, shade=True, legend=False
					)
ax.set(xlabel='Completeness/Sensitivity', ylabel='')
plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_1/2/error_analysis/' +
			'Compl_v_Sensi_kde.pdf',
			bbox_inches='tight'
			)
plt.clf()

# true
sns.set_context("paper")
sns.set(font_scale=1.5)

#ax = sns.jointplot(x='sensitivity', y='1-precision', data=merge_presen_df)
#ax.ax_marg_x.set_xlim(-5, 105)
#ax.ax_marg_y.set_ylim(-5, 105)
#x0, x1 = ax.ax_joint.get_xlim()
#y0, y1 = ax.ax_joint.get_ylim()
#lims = [max(x0, y0), min(x1, y1)]
#ax.ax_joint.plot([50.0, 50.0], [0.0, 100.0], ':k')
#ax.ax_joint.plot([0.0, 100.0], [50.0, 50.0], ':k')
ax = sns.kdeplot(merge_concomp_df['Contamination'].dropna(), color='blue', label='CheckM',
					bw=2, shade=True, legend=False
					)
ax = sns.kdeplot(merge_presen_df['1-precision'].dropna(), color='orange', label='True',
					bw=2, shade=True
					)
ax.set(xlabel='Contamination/1-Precision', ylabel='')
plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_1/2/error_analysis/' +
			'Contam_v_Precis_kde.pdf',
			bbox_inches='tight'
			)
plt.clf()


# Box plot of 51_1 checkm results
checkm_all_path = '/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_1/checkM_stdout/ALL_checkM_stdout.tsv'
checkm_all_df = pd.read_csv(checkm_all_path, sep='\t', header=0)

sub_checkm_df = checkm_all_df[['Bin Id', 'Completeness', 'Contamination',
								'Strain heterogeneity','data_type'
								]].set_index(['Bin Id', 'data_type'])
piv_sub_df = sub_checkm_df.stack().reset_index()
piv_sub_df.columns = ['sag_id', 'data_type', 'metric', 'score']

ax = sns.catplot(x="metric", y="score", hue='data_type', kind='box',
						data=piv_sub_df, aspect=2
						)
plt.plot([-1, 6], [25, 25], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 6], [50, 50], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 6], [75, 75], linestyle='--', alpha=0.3, color='k')

plt.ylim(-5, 105)
plt.title('SAG-plus')
ax._legend.set_title('')
plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_1/checkM_stdout/All_checkm_boxplox.svg',
			bbox_inches='tight'
			)
plt.clf()

# Box plot for 51_1 true error
true_all_path = '/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_1/trueerror_stdout/All_error_stats.tsv'
true_all_df = pd.read_csv(true_all_path, sep='\t', header=0)
#comb_true_df = true_all_df.loc[((true_all_df['algorithm'] == 'combined') |
#									(true_all_df['algorithm'] == 'MinHash')) &
#								(true_all_df['level'] == 'genus') &
#								(true_all_df['statistic'] != 'F1_score')]
fifty_true_df = true_all_df.loc[(true_all_df['data_type'] == 0.5) &
								(true_all_df['level'] == 'genus') &
								(true_all_df['statistic'] != 'F1_score')]
twenty_true_df = true_all_df.loc[(true_all_df['data_type'] == 0.2) &
								(true_all_df['level'] == 'genus') &
								(true_all_df['statistic'] != 'F1_score')]
ten_true_df = true_all_df.loc[(true_all_df['data_type'] == 0.1) &
								(true_all_df['level'] == 'genus') &
								(true_all_df['statistic'] != 'F1_score')]

sns.set_context("paper")
sns.set(font_scale=1.5)
ax = sns.catplot(x="statistic", y="score", hue='algorithm', kind='box',
						data=fifty_true_df, aspect=2
						)
plt.plot([-1, 6], [0.25, 0.25], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 6], [0.50, 0.50], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 6], [0.75, 0.75], linestyle='--', alpha=0.3, color='k')

plt.ylim(-0.05, 1.05)
plt.title('')
ax._legend.set_title('')
plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_1/trueerror_stdout/' +
			'50_trueerror_boxplox.pdf',
			bbox_inches='tight'
			)
plt.clf()

sns.set_context("paper")
sns.set(font_scale=1.5)
ax = sns.catplot(x="statistic", y="score", hue='algorithm', kind='box',
						data=twenty_true_df, aspect=2
						)
plt.plot([-1, 6], [0.25, 0.25], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 6], [0.50, 0.50], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 6], [0.75, 0.75], linestyle='--', alpha=0.3, color='k')

plt.ylim(-0.05, 1.05)
plt.title('')
ax._legend.set_title('')
plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_1/trueerror_stdout/' +
			'20_trueerror_boxplox.pdf',
			bbox_inches='tight'
			)
plt.clf()

sns.set_context("paper")
sns.set(font_scale=1.5)
ax = sns.catplot(x="statistic", y="score", hue='algorithm', kind='box',
						data=ten_true_df, aspect=2
						)
plt.plot([-1, 6], [0.25, 0.25], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 6], [0.50, 0.50], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 6], [0.75, 0.75], linestyle='--', alpha=0.3, color='k')

plt.ylim(-0.05, 1.05)
plt.title('')
ax._legend.set_title('')
plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_1/trueerror_stdout/' +
			'10_trueerror_boxplox.pdf',
			bbox_inches='tight'
			)
plt.clf()





'''
# Box plot of all checkm results
checkm_all_path = '/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_51/checkM_stdout/ALL_checkM_stdout.tsv'
checkm_all_df = pd.read_csv(checkm_all_path, sep='\t', header=0)

sub_checkm_df = checkm_all_df[['Bin Id', 'Completeness', 'Contamination',
								'Strain heterogeneity','data_type'
								]].set_index(['Bin Id', 'data_type'])
piv_sub_df = sub_checkm_df.stack().reset_index()
piv_sub_df.columns = ['sag_id', 'data_type', 'metric', 'score']

ax = sns.catplot(x="metric", y="score", hue='data_type', kind='box',
						data=piv_sub_df, aspect=2
						)
plt.plot([-1, 6], [25, 25], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 6], [50, 50], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 6], [75, 75], linestyle='--', alpha=0.3, color='k')

plt.ylim(-5, 105)
plt.title('SAG-plus')
ax._legend.set_title('')
plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_51/checkM_stdout/All_checkm_boxplox.svg',
			bbox_inches='tight'
			)
plt.clf()

# Box plot for all true error
true_all_path = '/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_51/trueerror_stdout/All_error_stats.tsv'
true_all_df = pd.read_csv(true_all_path, sep='\t', header=0)
comb_true_df = true_all_df.loc[(true_all_df['algorithm'] == 'combined') &
								(true_all_df['level'] == 'genus')]

ax = sns.catplot(x="statistic", y="score", hue='data_type', kind='box',
						data=comb_true_df, aspect=2
						)
plt.plot([-1, 6], [0.25, 0.25], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 6], [0.50, 0.50], linestyle='--', alpha=0.3, color='k')
plt.plot([-1, 6], [0.75, 0.75], linestyle='--', alpha=0.3, color='k')

plt.ylim(-0.05, 1.05)
plt.title('SAG-plus')
ax._legend.set_title('')
plt.savefig('/home/rmclaughlin/Ryan/SAG-plus/CAMI_I_HIGH/sag_redux/51_51/trueerror_stdout/All_trueerror_boxplox.svg',
			bbox_inches='tight'
			)
plt.clf()
'''
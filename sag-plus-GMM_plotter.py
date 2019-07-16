

def draw_ellipse(position, covariance, ax=None, **kwargs):
	"""Draw an ellipse with a given position and covariance"""
	ax = ax or plt.gca()
	
	# Convert covariance to principal axes
	if covariance.shape == (2, 2):
		U, s, Vt = np.linalg.svd(covariance)
		angle = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
		width, height = 2 * np.sqrt(s)
	else:
		angle = 0
		width, height = 2 * np.sqrt(covariance)
	
	# Draw the Ellipse
	for nsig in range(1, 3):
		nsig_width = nsig * width
		nsig_height = nsig * height
		ax.add_patch(Ellipse(position, nsig_width, nsig_height,
							 angle, **kwargs))

	return position, nsig_width, nsig_height, angle


def plot_ellispe_membership(df, plot_save_path, mean, covar):
	# Draw ellispe that colors by membership
	position, width, height, angle = draw_ellipse(mean, covar,
													edgecolor='green',
													linewidth=1, alpha=0.1
													)
	cos_angle = np.cos(np.radians(180.-angle))
	sin_angle = np.sin(np.radians(180.-angle))
	xc = df[df.columns[0]] - position[0]
	yc = df[df.columns[1]] - position[1]
	xct = xc * cos_angle - yc * sin_angle
	yct = xc * sin_angle + yc * cos_angle 
	rad_cc = (xct**2/(width/2.)**2) + (yct**2/(height/2.)**2)
	membr_category = []
	for r in rad_cc:
		if r <= 1.:
			# point in ellipse
			membr_category.append('SAG+')
		else:
			# point not in ellipse
			membr_category.append('Not SAG')
	if membr_category[0] == 'SAG+':
		pal = ['green', 'gray']
	else:
		pal = ['gray', 'green']
	ax = sns.scatterplot(x=df[df.columns[0]], y=df[df.columns[1]],
							hue=membr_category,	palette=pal
							)
	plt.gca().set_aspect('equal', 'datalim')
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.savefig(plot_save_path, bbox_inches="tight")
	plt.clf()
	# add membership to df
	#isSAG_col = '_'.join(['isSAG', df.columns[0], df.columns[1]])
	#df[isSAG_col] = [1 if x == 'SAG+' else 0 for x in membr_category]
	#df.drop(df.columns[:2], axis=1, inplace=True)

	#return df, isSAG_col


# Build plots for visualizing Tetra-Hz in 2-D space






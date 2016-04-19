# Discriminant Functions

Two R scripts described below.

# genData.R

Generates a dataset that could be troublesome for discriminant functions and/or linear probability models to correctly identify binary classified data due to the presence of high-leverage outliers.  The dataset is
created in such a way as to be easy for humans to visually identify a 100% correct boundary despite the model's difficulty in doing so.

### Public functions:
	generate.dataset()
	generate.troublesome.dataset()

Running the script will automatically call generate.troublesome.dataset() which will write the data to 'dataset.csv' and render a plot to 'dataPlot.pdf'.

To customize the dataset, call generate.dataset() directly.  See the code comments for parameters.


# loanData3C.R

Generates a dataset containing three classifications of loan applications: Approved, Denied, and Undecided.  It separates these classifications with two separate discriminant functions and generates a graph showing the data
points and boundary lines.  The first boundary line splits denied loans from all others because, presumably, it is better to deny a loan that should be approved than it is to approve a loan that should be denied.
Second, the remaining points are split with another boundary that separates approved from undecided loans.  Note that this boundary line does not extend beyond the first boundary line so as to avoid creating ambiguous
regions.

### Public function:
	run.analysis()

Running the script will automatically call run.analysis(), write the results to a file 'predictions.csv', and render a plot to 'discFunctions3C.pdf'.

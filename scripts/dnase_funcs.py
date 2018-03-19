"""
EXAMPLE DNASE LOOKUP FUNCTIONS TO TRY
"""

# Basic: One value per position in the motif.
def lookup_DNase1(seq, data, chrom, start, end, *args, **kwargs):
	if start > end:
		trace_start = end
		trace_end = start - 1
	else:
		trace_start = start
		trace_end = end - 1
	# Load bigwig file
	bwh = bx.bbi.bigwig_file.BigWigFile(open(data))
	# Extract region of interest
	trace = bwh.get_as_array(chrom, trace_start,trace_end)
	# If nothing is returned, return all zeros
	if trace is None:
		trace = np.zeros(trace_end - trace_start)
	# Replace nan's with zeros
	trace[np.isnan(trace)] = 0.0
	# Reverse if interested in the other strand
	if start > end: trace = trace[::-1]
	return trace.reshape([1,-1])

# Basic + Normalized
def lookup_DNase2(seq, data, chrom, start, end, *args, **kwargs):
	if start > end:
		trace_start = end
		trace_end = start - 1
	else:
		trace_start = start
		trace_end = end - 1
	# Load bigwig file
	bwh = bx.bbi.bigwig_file.BigWigFile(open(data))
	# Extract region of interest
	trace = bwh.get_as_array(chrom, trace_start,trace_end)
	# If nothing is returned, return all zeros
	if trace is None:
		trace = np.zeros(trace_end - trace_start)
	# Replace nan's with zeros
	trace[np.isnan(trace)] = 0.0
	# Reverse if interested in the other strand
	if start > end: trace = trace[::-1]
	# Normalize the trace
	trace -= np.min(trace)
	if not all(trace==0.0):
		trace /= np.max(trace)
	return trace.reshape([1,-1])

# Nucleotide resolution, including adjacent bases -- extended window
def lookup_DNase3(seq, data, chrom, start, end, *args, **kwargs):
	n_windows = 1
	if start > end:
		trace_start = end
		trace_end = start - 1
	else:
		trace_start = start
		trace_end = end - 1
	width = trace_end - trace_start
	# Load bigwig file
	bwh = bx.bbi.bigwig_file.BigWigFile(open(data))
	# Extract region of interest
	trace = bwh.get_as_array(chrom, trace_start-width*n_windows,trace_end+width*n_windows)
	# If nothing is returned, return all zeros
	if trace is None:
		trace = np.zeros(trace_end - trace_start + 2*width*n_windows)
	# Replace nan's with zeros
	trace[np.isnan(trace)] = 0.0
	# Reverse if interested in the other strand
	if start > end: trace = trace[::-1]
	return trace.reshape([2*n_windows+1,-1])

# Extended window, normalized
def lookup_DNase4(seq, data, chrom, start, end, *args, **kwargs):
	n_windows = 1
	if start > end:
		trace_start = end
		trace_end = start - 1
	else:
		trace_start = start
		trace_end = end - 1
	width = trace_end - trace_start
	# Load bigwig file
	bwh = bx.bbi.bigwig_file.BigWigFile(open(data))
	# Extract region of interest
	trace = bwh.get_as_array(chrom, trace_start-width*n_windows,trace_end+width*n_windows)
	# If nothing is returned, return all zeros
	if trace is None:
		trace = np.zeros(trace_end - trace_start + 2*width*n_windows)
	# Replace nan's with zeros
	trace[np.isnan(trace)] = 0.0
	# Reverse if interested in the other strand
	if start > end: trace = trace[::-1]
	# Normalize the trace
	trace -= np.min(trace)
	if not all(trace==0.0):
		trace /= np.max(trace)
	return trace.reshape([2*n_windows+1,-1])

# Bin local region
def lookup_DNase5(seq, data, chrom, start, end, *args, **kwargs):
	n_windows = 1
	if start > end:
		trace_start = end
		trace_end = start - 1
	else:
		trace_start = start
		trace_end = end - 1
	width = trace_end - trace_start
	# Load bigwig file
	bwh = bx.bbi.bigwig_file.BigWigFile(open(data))
	# Extract region of interest
	trace = bwh.get_as_array(chrom, trace_start-width*n_windows,trace_end+width*n_windows)
	# If nothing is returned, return all zeros
	if trace is None:
		trace = np.zeros(trace_end - trace_start + 2*width*n_windows)
	# Replace nan's with zeros
	trace[np.isnan(trace)] = 0.0
	# Reverse if interested in the other strand
	if start > end: trace = trace[::-1]

	# Bin the values
	binsize = n_windows*2+1
	trace = [np.sum(trace[i:i+binsize]) for i in range(0, len(trace), binsize)]
	trace = np.asarray(trace)

	return trace.reshape([1,-1])

# Bin local region, normalized
def lookup_DNase6(seq, data, chrom, start, end, *args, **kwargs):
	n_windows = 2
	if start > end:
		trace_start = end
		trace_end = start - 1
	else:
		trace_start = start
		trace_end = end - 1
	width = trace_end - trace_start
	# Load bigwig file
	bwh = bx.bbi.bigwig_file.BigWigFile(open(data))
	# Extract region of interest
	trace = bwh.get_as_array(chrom, trace_start-width*n_windows,trace_end+width*n_windows)
	# If nothing is returned, return all zeros
	if trace is None:
		trace = np.zeros(trace_end - trace_start + 2*width*n_windows)
	# Replace nan's with zeros
	trace[np.isnan(trace)] = 0.0
	# Reverse if interested in the other strand
	if start > end: trace = trace[::-1]

	# Bin the values
	binsize = n_windows*2+1
	trace = [np.sum(trace[i:i+binsize]) for i in range(0, len(trace), binsize)]
	trace = np.asarray(trace)

	# Normalize the trace
	trace -= np.min(trace)
	if not all(trace==0.0):
		trace /= np.max(trace)
	return trace.reshape([1,-1])
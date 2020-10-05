
import os
import sys
#sys.path.append('../')
import numpy as np
import scipy.signal as sp_signal
import scipy.stats as sp_stats
import pyphi
import itertools
from fly_phi_pyphi1_2_0 import *

def phi_compute(tpm, state_counters, nValues, out_dir, out_file):
	# Assumes tpm is state-by-node
	# Intputs:
	#	tpm = state-by-node TPM
	#	state_counters = vector of state counts
	#	nValues = number of states each element can take
	#	out_dir = string; directory to output results file
	#	out_name = string; name of output files
	
	# Build pyphi network
	network = pyphi.Network(tpm)
	# Determine number of system elements
	nChannels = np.shape(tpm)[1]
	# Determine number of system states
	n_states = nValues ** nChannels

	# Results structure
	phi = dict()

	# Initialise results storage structures
	phi_value = 0; # Doesn't need initialisation
	mips = np.empty((n_states), dtype=tuple)
	big_mips = np.empty((n_states), dtype=object)
	state_phis = np.zeros((n_states))

	# sys.exit()

	# Calculate all possible phi values (number of phi values is limited by the number of possible states)
	for state_index in range(0, n_states):
		#print('State ' + str(state_index))
		# Figure out the state
		#state = pyphi.convert.loli_index2state(state_index, nChannels)
		state = pyphi.convert.le_index2state(state_index, nChannels)
		# As the network is already limited to the channel set, the subsystem would have the same nodes as the full network
		subsystem = pyphi.Subsystem(network, state, network.node_indices)
		
		#sys.exit()
		
		# Compute phi values for all partitions
		big_mip = pyphi.compute.sia(subsystem)
		
		# Store phi and associated MIP
		state_phis[state_index] = big_mip.phi
		mips[state_index] = big_mip.cut
		
		# MATLAB friendly storage format (python saves json as nested dict)
		big_mip_json = big_mip.to_json()
		# Sort each constellation by their mechanisms
		big_mip_json['partitioned_'] = sort_constellation_by_mechanism(big_mip_json['ces'])
		big_mip_json['unpartitioned_constellation'] = sort_constellation_by_mechanism(big_mip_json['partitioned_ces'])
		# Store big_mip
		big_mips[state_index] = big_mip_json
		
		print('State ' + str(state_index) + ' Phi=' + str(big_mip.phi))
		#print(big_mip)

	phi_total = 0
	for state_index in range(0, n_states):
		if state_counters[state_index] > 0:
			# Add phi to total, weighted by the number of times the state occurred
			phi_total += state_phis[state_index] * state_counters[state_index]
	phi_value = phi_total / np.sum(state_counters)

	phi['phi'] = phi_value
	"""
	phi['state_counters'] = state_counters
	phi['big_mips'] = big_mips
	phi['state_phis'] = state_phis
	phi['tpm'] = tpm
	phi['mips'] = mips
	"""
	# Save
	save_mat(out_dir+out_file, {'phi': phi})
	print('saved ' + out_dir + out_file)

# Setup ############################################################################

# Parameters for loading data and calculating
tpm_type = sys.argv[1]

# TPM locations
tpm_dir = "tpms/" + tpm_type + "/"

# Output location
results_directory = "results/split/" + tpm_type + "/"
if not os.path.exists(results_directory):
	os.makedirs(results_directory)

# Loop through TPMs ############################################################################

for tpm_file in os.listdir(tpm_dir):
	
	if tpm_file != "params.mat":
		# Load TPM
		loaded = load_mat(tpm_dir + tpm_file)
		tpm = loaded['tpm']
		state_counters = loaded['state_counters']
		nValues = loaded['nValues'][0][0]
		tpm_formatted = pyphi.convert.to_2dimensional(pyphi.convert.state_by_state2state_by_node(tpm)) # assumes loaded TPMs are state-by-state
		#tpm_formatted = np.reshape(tpm_formatted, (16,4))
		#print(tpm_formatted)
		phi_compute(tpm_formatted, state_counters, nValues, results_directory, tpm_file)

"""
# Single PHI calculation to see all steps ###########################################################
# ~.~.~.~.~.~.~.~
tpm_file = 'tau1_thresh1_run1.mat'
# ~.~.~.~.~.~.~.~

loaded = load_mat(tpm_dir + tpm_file)
tpm = loaded['tpm']
state_counters = loaded['state_counters']
nValues = loaded['nValues'][0][0]
tpm_formatted = pyphi.convert.state_by_state2state_by_node(tpm) # assumes loaded TPMs are state-by-state
twod_format_tpm = pyphi.convert.to_2dimensional((tpm_formatted))
print(twod_format_tpm)



phi_compute(twod_format_tpm, state_counters, nValues, results_directory, tpm_file)

"""
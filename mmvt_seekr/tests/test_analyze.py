"""
tests for mmvt_seekr analysis module.
"""

import pickle
import mmvt_seekr
import pytest
import os


def test_analyze_kinetics():

	picklefile = open('test_data/test_model.pkl', 'rb')
	model = pickle.load(picklefile)
	
	p_equil, N, R, T, T_tot, Q, n_conv, r_conv, k_cell= mmvt_seekr.analyze.analyze_kinetics(
		model, [0], verbose=True,)
	MFPT = T[0]
	k_off = 1/MFPT
	assert T_tot == pytest.approx(7.856675637067191e-08)
	assert N[2,1] == pytest.approx(1877485711.4764836)
	assert k_off == pytest.approx(6101107.593969266)

def test_k_on_from_bd():

	picklefile = open('test_data/test_model.pkl', 'rb')
	model = pickle.load(picklefile)

	picklefile = open('test_data/test_Q.pkl', 'rb')
	Q = pickle.load(picklefile)

	os.chdir('test_data/test_filetree')
	k_on = mmvt_seekr.analyze.calc_kon_from_bd(model, [0], Q)
	assert k_on == pytest.approx(4289609646.4275894)
	os.chdir('../../')



def test_big_check_milestone_convergence():

	picklefile = open('test_data/test_model.pkl', 'rb')
	model = pickle.load(picklefile)

	N_conv, R_conv, k_cell_conv, p_equil_conv, k_conv, k_on_conv, conv_intervals = mmvt_seekr.analyze.check_milestone_convergence(
		model, [0], 500000, 5000000, 46288850,) 

	assert len(conv_intervals) == 92
	assert k_conv[12] == pytest.approx(5210133.55768566)

	min_anchor_times = mmvt_seekr.analyze.calc_RMSD_conv(model, 
		N_conv, R_conv, conv_intervals, 30, 0.05, 20)

	assert min_anchor_times[4] == 45000000

	p_equil, N, R, T, T_tot, Q, n_conv, r_conv, k_cell= mmvt_seekr.analyze.analyze_kinetics(
		model, [0], max_steps = min_anchor_times)

	MFPT = T[0]
	k_off = 1/MFPT

	assert k_off == pytest.approx(5925305.798841562)




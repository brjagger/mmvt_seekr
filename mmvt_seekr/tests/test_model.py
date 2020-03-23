"""
tests for mmvt_seekr model() module.
"""

# Import package, test suite, and other packages as needed
import mmvt_seekr
import pytest
import sys

def test_make_model():
    """Sample test, will always pass so long as import statement worked"""
    model, max_steps = mmvt_seekr.model.make_model(milestone_filename="test_data/milestones.xml", verbose=False)
    assert model.md_time_factor == 456 #check reading of initial xml info
    assert model.sites[0].num_anchors == 2 #check all available anchor info is read from XML
    assert max_steps == 7286640 #make sure transition data was parsed correctly and max steps determined
    assert model.sites[0].anchors[1].transitions[1].time == 6520.0 #check individual transitions are populated

def test_get_bd_transition_statistics():
	model, max_steps = mmvt_seekr.model.make_model(milestone_filename="test_data/milestones.xml", verbose=False)
	model.bd_milestone.directory = "test_data/test_filetree/bd_milestone"
	bd_counts, bd_total_counts, bd_total_times, bd_avg_times = model.bd_milestone.get_bd_transition_statistics(results_filename="results.xml", bd_time=0.0)
	assert bd_counts[model.bd_milestone.index]['6'] == 811673
#!/usr/bin/python

"""
analyze.py
Simulation Enabled Estimation of Kinetic Rates (SEEKR) is a tool that facilitates the preparation, running and analysis of multiscale MD/BD/Milestoning  simulations for the calculation of protein-ligand binding kinetics.

Performs kinetic analysis, including calculation of rate matrix, MFPT (on and off rates), and Milestone free energy

Parameters
                ----------
                calc_type: string default= "off" 
                                type of calculation, options "on" or "off"
                model: object Required
                        the SEEKR milestoning model
                bound_dict: dictionary Required
                        dictionary of bound states in the milestoning model



                Returns
                -------
                model : class 
                                contains all required information for milestoning analysis
"""

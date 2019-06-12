%
%	EPG - Extended Phase Graph Analysis Functions
%
%	These functions allow pulse sequence analysis
%	using the Extended Phase Graph (EPG) techniques - 
%	see Hennig, Weigel, and other articles.
%
%	EPG is based on dephased states, often called
%	F+, F- and Z, where the "dephasing order" indicates
%	how many cycles of dephasing the spins have undergone.
%
%
%	FUNCTIONS:	(help epg_xx to see syntax)
%
%	epg_rf	- propagate states through an RF pulse.
%	epg_grad - propagate through positive gradient only.
%	epg_mgrad - propagate through negative gradient only.
%	epg_grelax - propagate through relaxation and diffusion 
%			with or without a constant gradient.
%	epg_plot - simple plot of F,Z states.
%	epg_plotcomp - plot of Mz,My,Mz components across a voxel
%	epg_plot3Df - 3D plot of F states as spins in Mx,My,Mz space
%	epg_trim - remove states no longer needed.
%
%	epg_FZ2spins - convert EPG F,Z states to a set of M vectors
%	epg_spins2FZ - convert M vectors to EPG F,Z states
%
%	epg_gradanim - animate the transition due to a positive gradient
%
%	** - Requires arrow3D package!
%	epg_show** - make a multiple-plot to display states graphically
%	epg_showstate** - make a single plot of the magnetization
%	epg_showorder** - make a 2x2 diagram of nth F/Z state
%
%	EXAMPLES: 
%
%	epg_cpmg - CPMG sequence simulations
%	epg_spgr - SPGR (RF-spoiled) sequence simulation
%	epg_gre	- Gradient spoiled sequence simulation
%	epg_sediff - Spin Echo Diffusion
%	epg_stim - Simple stimulated echo example
%	epg_bssfp - bSSFP example





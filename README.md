# Agent-Based and PDE Models for Foraging Locusts

This repo contains codes used in creating the work of
- A J Bernoff, M Culshaw-Maurer, R A Everett, M E Hohn, W C Strickland, J Weinburd. *Agent-based and continuous models of hopper bands for the Australian plague locust: How resource consumption mediates pulse formation and geometry*. [(Preprint)](https://arxiv.org/abs/1910.14553)

The repo contains one Matlab script and three directories:
- `video_ABM.m`
- ABM_in_R
- Publication2020_MATLAB
- Sobol

### `video_ABM.m`
This self-contained script can be run in Matlab to simulate our agent-based model and save a video of the simulation. Adjust the variable `plot_frames` to control the frequency of plotting and saving frames to the video.

### ABM_in_R
Another version of the ABM, written in R. Includes two subdirectories:
- misc_scripts (includes code used to generate Fig 1)
- model_scripts (includes a full version of the ABM)

### Publication2020_MATLAB
Everything needed to reproduce Figs 2-8. Each Matlab script `fig#.m` generates a figure from the paper. By setting the value of the flag `new_data`, the script will either load precisely the data used in the paper (`=0`) or generate new data (`=1`) and overwrite existing. Includes two subdirectories:
- data (contains numerically generated data for figures, `.mat` format)
- functions (includes Matlab scripts called by the figure scripts)

### Sobol
All the code used to conduct the Sobol sensitivity analysis and generate Figs 9-11.

## Acknowledgements
This material grew out of a collaboration that began at a Mathematics Research Communities of the American Mathematical Society on Agent-based Modeling in Biological and Social Systems (2018), under National Science Foundation grant DMS-1321794. Additional support came from the Institute for Advanced Study Summer Collaborators Program, two Simons Foundation Collaboration Grants (Bernoff and Strickland), and an NSF Mathematical Sciences Postdoctoral Research Fellowship grant DMS-1902818 (Weinburd).
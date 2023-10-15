This repository includes custom codes for the results and figures in the article: Indirect reciprocity with assessments of group reputation.

The basic folder contains the color blendent folder and four basic matrics files describing different situations of population compositions and game compositions.

The advanced folder contains the custom codes for the article. 
Default.m contains default settings of the experiment and should be adjusted and run before any test.  rep_evol_sym_cr.m is the main program for the reputation dynamics, including inner_game2.m and outer_eval.m which are the functions describing the reputation updating process within and out of the game group. 
After collecting the cooperation rates, use cr2pofv.m to derive the expectations of payoff under different population compositions. 
This step involves payoff.m as the basic payoff calculation in one game. 
Then, use slc_mut_equ.m to calculate the selection-mutation equilibrium, as well as fixprob.m and fixprob_calculate.m for  fixation probabilities, based on the payoff results. 
mugg0.m corresponds to the situations with significant mutation rates. 
Finally, data_vis.m contains the data visualization code of each experiment, appearing as a series of code blocks. This file should be applied within in each block, instead of running through all the lines.

The lambda_evolve folder includes codes for the evolutionary process of group assessment criterion.
main.m is the main body code, including lambda_evol_cr.m as the evolutionary function. 
lainner_game.m and laouter_eval.m which are the functions describing the reputation updating process within and out of the game group.
lacr2pofv use cooperation rates to compute expectations of payoff.
laAvFr.m calculates the selection-mutation equilibrium when mutations are rare.
laSMEq.m calculates the selection-mutation equilibrium with higher mutation rates.
Finally, data_vis_laev.m contains the data visualization code of each experiment, appearing as a series of code blocks. This file should be applied within in each block, instead of running through all the lines.  

Notice that, each position of file path in the code file mentioned above is an example of data path with respect to the  corresponding experiment. Please replace them with legal format of file paths before running the test.

The results fold contains the formal results for each part of the main text. Please contact the author for any question about detailed information.

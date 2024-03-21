import argparse
import os

import runBMCA

# ant, data_file, output_dir, n_iter,

# all data 
runBMCA.runBMCA('data/interim/Antimony/JSexample22.ant',
                'data/interim/generated_data/JSexample22-noReg/JSexample22',
                'data/interim/generated_data/JSexample22-noReg', 5)
runBMCA.runBMCA('data/interim/Antimony/JSexample22_reg1.ant',
                'data/interim/generated_data/JSexample22-reg1/JSexample22_reg1',
                'data/interim/generated_data/JSexample22-reg1', 5)
runBMCA.runBMCA('data/interim/Antimony/JSexample22_reg2.ant',
                'data/interim/generated_data/JSexample22-reg2/JSexample22_reg2',
                'data/interim/generated_data/JSexample22-reg2', 5)

# omitted fluxes
runBMCA.runBMCA('data/interim/Antimony/JSexample22.ant',
                'data/interim/generated_data/JSexample22-noReg/JSexample22',
                'data/interim/generated_data/JSexample22-noReg', 5,
                omit='fluxes')
runBMCA.runBMCA('data/interim/Antimony/JSexample22_reg1.ant',
                'data/interim/generated_data/JSexample22-reg1/JSexample22_reg1',
                'data/interim/generated_data/JSexample22-reg1', 5,
                omit='fluxes')
runBMCA.runBMCA('data/interim/Antimony/JSexample22_reg2.ant',
                'data/interim/generated_data/JSexample22-reg2/JSexample22_reg2',
                'data/interim/generated_data/JSexample22-reg2', 5,
                omit='fluxes')

# omitted enzymes
runBMCA.runBMCA('data/interim/Antimony/JSexample22.ant',
                'data/interim/generated_data/JSexample22-noReg/JSexample22',
                'data/interim/generated_data/JSexample22-noReg', 5,
                omit='enzymes')
runBMCA.runBMCA('data/interim/Antimony/JSexample22_reg1.ant',
                'data/interim/generated_data/JSexample22-reg1/JSexample22_reg1',
                'data/interim/generated_data/JSexample22-reg1', 5,
                omit='enzymes')
runBMCA.runBMCA('data/interim/Antimony/JSexample22_reg2.ant',
                'data/interim/generated_data/JSexample22-reg2/JSexample22_reg2',
                'data/interim/generated_data/JSexample22-reg2', 5,
                omit='enzymes')

# omitted internal
runBMCA.runBMCA('data/interim/Antimony/JSexample22.ant',
                'data/interim/generated_data/JSexample22-noReg/JSexample22',
                'data/interim/generated_data/JSexample22-noReg', 5,
                omit='internal')
runBMCA.runBMCA('data/interim/Antimony/JSexample22_reg1.ant',
                'data/interim/generated_data/JSexample22-reg1/JSexample22_reg1',
                'data/interim/generated_data/JSexample22-reg1', 5,
                omit='internal')
runBMCA.runBMCA('data/interim/Antimony/JSexample22_reg2.ant',
                'data/interim/generated_data/JSexample22-reg2/JSexample22_reg2',
                'data/interim/generated_data/JSexample22-reg2', 5,
                omit='internal')    

# omitted external
runBMCA.runBMCA('data/interim/Antimony/JSexample22.ant',
                'data/interim/generated_data/JSexample22-noReg/JSexample22',
                'data/interim/generated_data/JSexample22-noReg', 5,
                omit='external')
runBMCA.runBMCA('data/interim/Antimony/JSexample22_reg1.ant',
                'data/interim/generated_data/JSexample22-reg1/JSexample22_reg1',
                'data/interim/generated_data/JSexample22-reg1', 5,
                omit='external')
runBMCA.runBMCA('data/interim/Antimony/JSexample22_reg2.ant',
                'data/interim/generated_data/JSexample22-reg2/JSexample22_reg2',
                'data/interim/generated_data/JSexample22-reg2', 5,
                omit='external')    
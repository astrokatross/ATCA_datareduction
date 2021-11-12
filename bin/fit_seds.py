#!/usr/bin/python3
# This script fits the entire MWA and ATCA seds
# By K.Ross 29/09/21

import analysis_functs

# Source/run information
save_dir = "/data/ATCA/analysis/"
data_dir = "/data/ATCA/ATCA_datareduction/"
gleam_targets = ["GLEAM J001513-472706", "GLEAM J015445-232950", "GLEAM J020507-110922", "GLEAM J021246-305454", "GLEAM J022744-062106", "GLEAM J024838-321336", "GLEAM J032213-462646", "GLEAM J032836-202138", "GLEAM J033023-074052", "GLEAM J042502-245129", "GLEAM J044033-422918", "GLEAM J044737-220335", "GLEAM J052824-331104", "GLEAM J223933-451414", "GLEAM J224408-202719"]


for i in range(len(gleam_targets)):
    gleam_tar = gleam_targets[i]
    analysis_functs.run_everything(save_dir, data_dir, gleam_tar)
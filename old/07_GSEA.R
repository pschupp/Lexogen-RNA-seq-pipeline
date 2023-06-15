source('~/code/git/GSEA/GSEAfxsV3.r')
MyGSHGloop(kmecut1='fdr',exclude="none",pvalcut1=NULL, moduleDir='/home/patrick/opentrons_testing_run/07_sample_networks/ruvg_protein_Modules', geneSets='/home/patrick/code/GSEA/genesets_slim', outputName='ruvg_protein')

source('~/code/git/GSEA/GSEAfxsV3.r')
MyGSHGloop(kmecut1='fdr',exclude="none",pvalcut1=NULL, moduleDir='/home/patrick/opentrons_testing_run/07_sample_networks/ruvg_prot_qn_Modules', geneSets='/home/patrick/code/GSEA/genesets_slim', outputName='ruvg_prot_qn')

source('~/code/git/GSEA/GSEAfxsV3.r')
MyGSHGloop(kmecut1='fdr',exclude="none",pvalcut1=NULL, moduleDir='/home/patrick/opentrons_testing_run/07_sample_networks/ruvg_all_pc1_Modules', geneSets='/home/patrick/code/GSEA/genesets_slim', outputName='ruvg_prot_all_pc1')

source('~/code/git/GSEA/GSEAfxsV3.r')
MyGSHGloop(kmecut1='fdr',exclude="none",pvalcut1=NULL, moduleDir='/home/patrick/opentrons_testing_run/07_sample_networks/ruvg_all_Modules', geneSets='/home/patrick/code/GSEA/genesets_slim', outputName='ruvg_all')

# best outome is from ruvg_prot_qn_Modules and ruvg_all_pc1_Modules:/mnt/fidelity/opentrons_testing_run/07_sample_networks/ruvg_all_pc1_Modules/Bicor-None_signum0.497_minSize20_merge_ME_0.8_31558
/mnt/fidelity/opentrons_testing_run/07_sample_networks/ruvg_prot_qn_Modules/Bicor-None_signum0.625_minSize8_merge_ME_0.8_17671


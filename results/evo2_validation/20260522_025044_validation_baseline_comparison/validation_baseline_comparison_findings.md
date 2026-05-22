# Does Evo2 Generalise? Validating Functional Classification AUROC on Glycolysis vs T-cell Signaling

## Hypothesis

Run the baseline MLP classifier (frozen Evo2-7B layer blocks.28.mlp.l3, 256->128->1 MLP) on glycolysis vs T-cell signaling genes. Report AUROC, accuracy, F1. Compare directly to the DNA repair vs TSG result (0.78 AUROC). If AUROC is similarly high, it suggests Evo2 broadly encodes functional class information. If AUROC is lower, the DNA repair/TSG result may reflect something specific to those classes (e.g. cancer-related sequence features). Also compute raw cosine similarity between class centroids as a baseline — does the unsupervised geometry show more separation here than in the DNA repair/TSG case?

## Results

- **final_train_loss**: 3.0945
- **final_val_accuracy**: 0.7895
- **final_val_auroc**: 0.6179
- **final_val_f1**: 0.8667
- **best_val_auroc**: 0.6464
- **best_val_accuracy**: 0.7895
- **n_train**: 151.0000
- **n_test**: 38.0000
- **n_glycolysis**: 67.0000
- **n_tcell**: 122.0000
- **input_dim**: 4096.0000
- **centroid_cosine_similarity**: 0.9986
- **centroid_cosine_distance**: 0.0014
- **within_glycolysis_cosine_sim**: 0.9554
- **within_tcell_cosine_sim**: 0.9519

## Scores

- Interestingness: 9/10
- Feasibility: 9/10
- Novelty: 8/10

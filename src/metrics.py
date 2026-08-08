import numpy as np

def calculate_mcc(true_matrix: np.ndarray, pred_matrix: np.ndarray) -> float:
    """
    Calculates the Matthews Correlation Coefficient (MCC) between ground truth
    and predicted base-pair binary adjacency matrices.
    """
    TP = np.sum((true_matrix == 1) & (pred_matrix == 1))
    TN = np.sum((true_matrix == 0) & (pred_matrix == 0))
    FP = np.sum((true_matrix == 0) & (pred_matrix == 1))
    FN = np.sum((true_matrix == 1) & (pred_matrix == 0))

    numerator = (TP * TN) - (FP * FN)
    denominator = np.sqrt(float((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))

    if denominator == 0:
        return 0.0

    return numerator / denominator

def evaluate(
    true_matrix: np.ndarray,
    pred_matrix: np.ndarray,
) -> dict:
    """
    Evaluate a predicted RNA secondary structure.
    """

    TP = np.sum((true_matrix == 1) & (pred_matrix == 1))
    TN = np.sum((true_matrix == 0) & (pred_matrix == 0))
    FP = np.sum((true_matrix == 0) & (pred_matrix == 1))
    FN = np.sum((true_matrix == 1) & (pred_matrix == 0))

    mcc = calculate_mcc(
        true_matrix,
        pred_matrix,
    )

    return {

        "TP": int(TP),

        "TN": int(TN),

        "FP": int(FP),

        "FN": int(FN),

        "MCC": float(mcc),

    }
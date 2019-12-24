package DifferentialMethylation.PseudoPrior;

import java.util.ArrayList;

/**
 * estimate parameters of Beta distribution,
 * reference to https://pdfs.semanticscholar.org/ddc7/cb7bc07d44dc4b3395def445016ecc8c00cc.pdf
 */
public class PseudoBetaDistribution {
    public static double[] estimate(ArrayList<Double> dataList) {
        double mean = CommonMethod.calcMean(dataList);
        double variance = CommonMethod.calcVariance(dataList);

        double param1 = mean * (mean * (1 - mean) / Math.pow(variance, 2) - 1);
        double param2 = (1 - mean) * (mean * (1 - mean) / Math.pow(variance, 2) - 1);

        // avoid org.apache.commons.math3.exception.NotStrictlyPositiveException which caused by param1 and param2 equals 0
        param1 = (Math.abs(param1) < 0.00001)? 0.001: param1;
        param2 = (Math.abs(param2) < 0.00001)? 0.001: param2;

        return new double[] {param1, param2};
    }
}

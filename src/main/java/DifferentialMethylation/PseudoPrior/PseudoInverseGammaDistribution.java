package DifferentialMethylation.PseudoPrior;


import java.util.ArrayList;

/**
 * inverse Gamma distribution parameter estimation with moment method
 * reference to https://arxiv.org/pdf/1605.01019.pdf
 */
public class PseudoInverseGammaDistribution {

    public static double[] estimate(ArrayList<Double> dataList) {
        if (dataList.size() == 0)
            return null;
        double mean = CommonMethod.calcMean(dataList);
        double variance = CommonMethod.calcVariance(dataList);
        double shape = Math.pow(mean, 2) / variance + 2;
        double scale = mean * (Math.pow(mean, 2) / variance +1);

        return new double[] {shape, scale};
    }

}

package DifferentialMethylation.PseudoPrior;

import java.util.ArrayList;
import java.util.NoSuchElementException;

/**
 * estimate gamma distribution parameters with maximum likelihood method,
 * reference to https://en.wikipedia.org/wiki/Gamma_distribution#Maximum_likelihood_estimation
 */
@Deprecated
public class PseudoGammaDistribution {
    public static double[] estimate(ArrayList<Double> dataList) {
        double sum, logSum, sampleNum;
        try {
            sum = dataList.stream().reduce((x, y) -> x + y).get();
            logSum = dataList.stream().map(x -> Math.log(x)).reduce((x, y) -> x + y).get();
            sampleNum = dataList.size();
            double s = Math.log(sum / sampleNum) - logSum / sampleNum;
            double shape = ((3 - s) + Math.sqrt(Math.pow(s - 3, 2) + 24 * s)) / 12 * s;
            double scale = sum / (shape * sampleNum);

            // avoid org.apache.commons.math3.exception.NotStrictlyPositiveException which caused by shape and scale equals 0
            shape = (Math.abs(shape) < 0.00001)? 0.001: shape;
            scale = (Math.abs(scale) < 0.00001)? 0.001: scale;

            return new double[] {shape, scale};
        } catch (NoSuchElementException nee) {
            return new double[] {0, 0};
        }
    }
}

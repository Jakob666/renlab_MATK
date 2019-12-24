package DifferentialMethylation.PseudoPrior;

import java.util.ArrayList;

/**
 * estimate parameters of log normal distribution,
 * reference to https://scholarsarchive.byu.edu/cgi/viewcontent.cgi?article=2927&context=etd
 */
public class PseudoLogNormalDistribution {
    public static double[] estimate(ArrayList<Double> dataList) {
        double dataSum = dataList.stream().reduce((x, y) -> x + y).get();
        double squareSum = dataList.stream().map(x -> Math.pow(x, 2)).reduce((x, y) -> x + y).get();

        double miu = -1 * Math.log(squareSum) / 2 + 2 * Math.log(dataSum) - 3 * Math.log(dataList.size()) / 2;
        double sigma = Math.sqrt(Math.log(squareSum) - 2 * Math.log(dataSum) + Math.log(dataList.size()));

        // avoid org.apache.commons.math3.exception.NotStrictlyPositiveException which caused by miu and sigma equals 0
        miu = (Math.abs(miu) < 0.00001)? 0.001: miu;
        sigma = (Math.abs(sigma) < 0.00001)? 0.001: sigma;

        return new double[] {miu, sigma};
    }
}

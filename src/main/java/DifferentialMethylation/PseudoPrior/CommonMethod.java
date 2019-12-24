package DifferentialMethylation.PseudoPrior;

import java.util.ArrayList;

public class CommonMethod {

    public static double calcMean(ArrayList<Double> dataList) {
        double sum = dataList.stream().reduce((x, y) -> x+y).get();

        return sum / dataList.size();
    }

    public static double calcVariance(ArrayList<Double> dataList) {
        double mean = calcMean(dataList);
        double variance = dataList.stream().map(x -> Math.pow(x-mean, 2)).reduce((x, y) -> x+y).get() / dataList.size();
        // avoid zero
        if (Math.abs(variance) < 0.00001)
            variance = 0.0001;
        return variance;
    }
}

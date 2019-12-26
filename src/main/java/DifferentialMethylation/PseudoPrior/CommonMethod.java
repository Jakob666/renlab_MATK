package DifferentialMethylation.PseudoPrior;

import java.util.ArrayList;
import java.util.NoSuchElementException;

public class CommonMethod {

    public static double calcMean(ArrayList<Double> dataList) {
        try {
            double sum = dataList.stream().reduce((x, y) -> x+y).get();
            return sum / dataList.size();
        } catch (NoSuchElementException nee) {
            return 0;
        }
    }

    public static double calcVariance(ArrayList<Double> dataList) {
        try {
            double mean = calcMean(dataList);
            double variance = dataList.stream().map(x -> Math.pow(x-mean, 2)).reduce((x, y) -> x+y).get() / dataList.size();
            // avoid zero
            if (Math.abs(variance) < 0.00001)
                variance = 0.0001;
            return variance;
        } catch (NoSuchElementException nee) {
            return 0;
        }

    }
}

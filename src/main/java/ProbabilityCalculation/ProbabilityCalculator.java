package ProbabilityCalculation;

import DifferentialMethylation.ModelSelection;
import SeqDataModel.NegativeBinomialDistribution;

public class ProbabilityCalculator {

    /**
     * negative binomial distribution probability
     * @param readsCount reads count of each individual cover on a gene, shape 1 × individualNumber
     * @param readsExpectation reads count expectation of each individual cover on a gene, shape 1 × individualNumber
     * @param overdispersion reads count overdispersion of a gene
     * @return probability cumPro(p(rc_i|expect_i, overdispersion)), i=1,2,..., individualNumber
     */
    public static double logNegativeProbability(double[] readsCount, double[] readsExpectation, double overdispersion) {
        assert readsCount.length == readsExpectation.length;
        NegativeBinomialDistribution nbd;
        double proba = 0;
        for (int i=0; i<readsCount.length; i++) {
            nbd = new NegativeBinomialDistribution(readsExpectation[i], overdispersion);
            proba += nbd.logDensity(readsCount[i]);
        }
        nbd = null;

        return proba;
    }

    /**
     * negative binomial distribution probability
     * @param readsCount reads count of each individual cover on a gene, shape 1 × individualNumber
     * @param readsExpectation reads count expectation of each individual cover on a gene, shape 1 × individualNumber
     * @param overdispersion reads count overdispersion of a gene
     * @return probability cumPro(p(rc_i|expect_i, overdispersion)), i=1,2,..., individualNumber
     */
    public static double logNegativeProbability(int[] readsCount, double[] readsExpectation, double overdispersion) {
        assert readsCount.length == readsExpectation.length;
        NegativeBinomialDistribution nbd;
        double proba = 0;
        for (int i=0; i<readsCount.length; i++) {
            nbd = new NegativeBinomialDistribution(readsExpectation[i], overdispersion);
            proba += nbd.logDensity(readsCount[i]);
        }
        nbd = null;

        return proba;
    }
}

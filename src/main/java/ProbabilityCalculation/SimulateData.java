package ProbabilityCalculation;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

public class SimulateData {
    private int geneNumber, individualNumber;
    private double[] methylationLevel = new double[9];
    private double[] inputBkg, ipBkg;
    private double[][] inputSizeFactor, ipSizeFactor;
    private GammaDistribution ipOverdispersionDistribution, inputOverdispersionDistribution;
    private BetaDistribution methylationDistribution;

    public SimulateData(int geneNumber, int individualNumber, double bkgExpMean, double bkgExpStd, double methParam1, double methParam2,
                        double ipOverdispersionParam1, double ipOverdispersionParam2, double inputOverdispersionParam1, double inputOverdispersionParam2) {
        this.geneNumber = geneNumber;
        this.individualNumber = individualNumber;

    }



    /**
     * methylation level 0.1~0.9
     */
    private void methLevelRange() {
        for (int i=0; i<methylationLevel.length; i++) {
            methylationLevel[i] = 0.1 * (i+1);
        }
    }
}
